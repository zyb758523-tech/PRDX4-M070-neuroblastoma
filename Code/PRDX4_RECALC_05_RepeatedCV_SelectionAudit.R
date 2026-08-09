rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20260808)

## ============================================================
## PRDX4_RECALC_05_RepeatedCV_SelectionAudit.R
##
## Strict model-selection stability audit
##   frozen 24-gene cross-platform panel from completed SC28 RECALC 04C
##   -> E-MTAB-179 only: 10 repeats x stratified 5-fold OOF comparison
##   -> all historical 101 combinations compete in every repeat
##   -> summarize C-index/AUC rank stability and feature stability
##   -> DO NOT lock a final model and DO NOT read SEQC outcomes
##
## Important safeguards
##   1. Historical ComboID is preserved; M050 = SuperPC -> Ridge.
##   2. The exact 24 genes already confirmed measurable in both cohorts are frozen.
##      PRDX4/LDHA are not appended and SEQC is not opened in this script.
##   3. Validation-fold survival outcomes are never supplied to prediction functions.
##   4. Ranking uses continuous discrimination metrics; Cox/log-rank P values do not
##      break ties or determine the winner.
##   5. Repeat/fold checkpoints allow safe resume after interruption.
## ============================================================

root_dir <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"
seed_main <- 20260808
k_folds <- 5L
repeats_n <- 10L
resume_existing_folds <- TRUE

emtab_expr_file <- file.path(root_dir, "SC13A_MANUAL_V6_EMTAB179_ADF_TabFix", "Tables", "05_EMTAB179_expr_gene_symbol_matched.rds")
emtab_clin_file <- file.path(root_dir, "SC13A_MANUAL_V5_EMTAB179_SemicolonFix", "Tables", "04_EMTAB179_clin_matched_raw.csv")
out_dir <- file.path(root_dir, "SC28_RECALC_05_RepeatedCV_SelectionAudit_20260808")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Checkpoints"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "Logs", "SC28_RECALC_05_log.txt")
sink(log_file, split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("SC28 RECALC 05: 10x5 repeated-CV selection stability audit\n")
cat("Start:", as.character(Sys.time()), "\n")
cat("Seed:", seed_main, " folds:", k_folds, " repeats:", repeats_n, "\n")
cat("============================================================\n\n")

## ---------------------------- helpers ----------------------------
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

clean_event <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(y))
  out[y %in% c("1","dead","deceased","death","event","yes","true","died","deceased due to disease")] <- 1
  out[y %in% c("0","alive","living","censored","no","false","alive without event","alive with relapse/primary tumor")] <- 0
  z <- suppressWarnings(as.numeric(y))
  out[is.na(out) & !is.na(z)] <- z[is.na(out) & !is.na(z)]
  out
}

as_expr_matrix <- function(x) {
  if (is.matrix(x)) mat <- x else if (is.data.frame(x)) {
    df <- x
    first_num <- suppressWarnings(as.numeric(as.character(df[[1]])))
    if (sum(is.na(first_num)) > length(first_num) * 0.5) {
      rownames(df) <- make.unique(as.character(df[[1]])); df <- df[, -1, drop = FALSE]
    }
    mat <- as.matrix(as.data.frame(lapply(df, safe_num), check.names = FALSE))
    rownames(mat) <- rownames(df)
  } else stop("Unsupported expression object")
  mode(mat) <- "numeric"; mat
}

collapse_duplicate_genes <- function(mat) {
  rn <- toupper(trimws(as.character(rownames(mat))))
  keep <- !is.na(rn) & rn != ""
  mat <- mat[keep, , drop = FALSE]; rn <- rn[keep]; rownames(mat) <- rn
  if (!anyDuplicated(rn)) return(mat)
  sp <- split(seq_along(rn), rn)
  z <- do.call(rbind, lapply(sp, function(i) colMeans(mat[i, , drop = FALSE], na.rm = TRUE)))
  rownames(z) <- names(sp); z
}

find_time_event_cols <- function(clin) {
  cn <- colnames(clin); low <- tolower(cn)
  tp <- c("OS_time","OS.time","os_time","os.time","OS_time_days","os_time_days","OVERALLSURVIVAL","overall_survival","survival_time")
  ep <- c("OS_event","OS.event","os_event","os.event","DEATHOFDISEASE","death_of_disease","overall_survival_event","vital_status","event","status")
  tc <- tp[tp %in% cn][1]; ec <- ep[ep %in% cn][1]
  if (is.na(tc)) for (p in c("os.*time","overall.*survival","survival.*time","follow.*up","time")) {
    h <- cn[grepl(p, low)]; if (length(h)) for (q in h) if (sum(!is.na(safe_num(clin[[q]]))) >= 30) { tc <- q; break }
    if (!is.na(tc)) break
  }
  if (is.na(ec)) for (p in c("os.*event","death.*disease","death","vital.*status","event","status")) {
    h <- cn[grepl(p, low)]; if (length(h)) for (q in h) {
      e <- clean_event(clin[[q]]); if (sum(!is.na(e)) >= 30 && length(unique(e[!is.na(e)])) >= 2) { ec <- q; break }
    }
    if (!is.na(ec)) break
  }
  list(time = tc, event = ec)
}

scale_train_apply <- function(train_mat, valid_mat) {
  mu <- rowMeans(train_mat, na.rm = TRUE)
  sdv <- apply(train_mat, 1, sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- 1
  tr <- sweep(sweep(train_mat, 1, mu, "-"), 1, sdv, "/")
  va <- sweep(sweep(valid_mat, 1, mu, "-"), 1, sdv, "/")
  list(train_z = tr, valid_z = va, center = mu, scale = sdv)
}

make_folds <- function(event, k = 5, seed = 1) {
  set.seed(seed); f <- integer(length(event))
  for (lev in sort(unique(event))) {
    ix <- which(event == lev); ix <- sample(ix)
    f[ix] <- rep(seq_len(k), length.out = length(ix))
  }
  f
}

univ_rank <- function(train_z, genes, time, event) {
  ans <- lapply(genes, function(g) {
    d <- data.frame(time=time, event=event, x=as.numeric(train_z[g,]))
    fit <- tryCatch(survival::coxph(survival::Surv(time,event) ~ x, data=d), error=function(e) NULL)
    if (is.null(fit)) return(NULL); s <- summary(fit)
    data.frame(Gene=g, HR=s$coefficients[1,"exp(coef)"], P=s$coefficients[1,"Pr(>|z|)"])
  })
  ans <- ans[!vapply(ans,is.null,logical(1))]; if (!length(ans)) return(data.frame())
  z <- do.call(rbind, ans); z[order(z$P, -abs(log(z$HR))), , drop=FALSE]
}

orient_scores <- function(risk_train, risk_valid, time, event) {
  d <- data.frame(time=time,event=event,risk=as.numeric(risk_train))
  fit <- tryCatch(survival::coxph(survival::Surv(time,event) ~ risk, data=d), error=function(e) NULL)
  if (!is.null(fit) && is.finite(coef(fit)[1]) && coef(fit)[1] < 0) {
    risk_train <- -risk_train; risk_valid <- -risk_valid
  }
  list(train=risk_train, valid=risk_valid)
}

finalize_fit_output <- function(out, x, xv, time, event) {
  if (is.null(out) || is.null(out$rt) || is.null(out$rv)) return(NULL)
  out$rt <- safe_num(out$rt); out$rv <- safe_num(out$rv)
  if (length(out$rt) != nrow(x) || length(out$rv) != nrow(xv)) return(NULL)
  if (length(unique(out$rt[is.finite(out$rt)])) < 2 ||
      length(unique(out$rv[is.finite(out$rv)])) < 2) return(NULL)
  oo <- orient_scores(out$rt, out$rv, time, event)
  out$rt <- oo$train; out$rv <- oo$valid
  ## Put OOF scores from separately trained folds onto a common scale.
  ## Only the TRAINING-score distribution is used; validation outcomes and
  ## validation-score distribution are not used for this transformation.
  risk_mu <- mean(out$rt, na.rm=TRUE); risk_sd <- sd(out$rt, na.rm=TRUE)
  if (!is.finite(risk_mu) || !is.finite(risk_sd) || risk_sd <= 0) return(NULL)
  out$rt <- (out$rt-risk_mu)/risk_sd
  out$rv <- (out$rv-risk_mu)/risk_sd
  names(out$rt) <- rownames(x); names(out$rv) <- rownames(xv)
  out
}

## ---------------------------- packages ----------------------------
required <- c("survival","glmnet","quadprog")
optional <- c("timeROC","randomForestSRC","gbm","CoxBoost","superpc","survivalsvm","plsRcox")

## IMPORTANT FOR RSTUDIO:
## Do not install through a temporary loop variable such as `p` here.
## If RStudio restarts R to update a loaded dependency, that temporary
## variable disappears and the deferred install call would fail with
## "object 'p' not found".  Literal package names make a restart safe.
## After any RStudio-requested restart/install completes, click Source again.
if (!requireNamespace("survival", quietly=TRUE)) {
  cat("Installing missing package: survival\n")
  try(install.packages("survival", repos="https://cloud.r-project.org"), silent=TRUE)
}
if (!requireNamespace("glmnet", quietly=TRUE)) {
  cat("Installing missing package: glmnet\n")
  try(install.packages("glmnet", repos="https://cloud.r-project.org"), silent=TRUE)
}
if (!requireNamespace("quadprog", quietly=TRUE)) {
  cat("Installing missing package required by SurvivalSVM: quadprog\n")
  try(install.packages("quadprog", repos="https://cloud.r-project.org"), silent=TRUE)
}
if (!requireNamespace("timeROC", quietly=TRUE)) {
  cat("Installing missing package: timeROC\n")
  try(install.packages("timeROC", repos="https://cloud.r-project.org"), silent=TRUE)
}
if (!requireNamespace("randomForestSRC", quietly=TRUE)) {
  cat("Installing missing package: randomForestSRC\n")
  try(install.packages("randomForestSRC", repos="https://cloud.r-project.org"), silent=TRUE)
}
if (!requireNamespace("gbm", quietly=TRUE)) {
  cat("Installing missing package: gbm\n")
  try(install.packages("gbm", repos="https://cloud.r-project.org"), silent=TRUE)
}
if (!requireNamespace("CoxBoost", quietly=TRUE)) {
  cat("Installing missing package: CoxBoost\n")
  try(install.packages("CoxBoost", repos="https://cloud.r-project.org"), silent=TRUE)
}
if (!requireNamespace("superpc", quietly=TRUE)) {
  cat("Installing missing package: superpc\n")
  try(install.packages("superpc", repos="https://cloud.r-project.org"), silent=TRUE)
}
if (!requireNamespace("survivalsvm", quietly=TRUE)) {
  cat("Installing missing package: survivalsvm\n")
  try(install.packages("survivalsvm", repos="https://cloud.r-project.org"), silent=TRUE)
}

## plsRcox imports two Bioconductor packages. CRAN can install the plsRcox
## binary while warning that mixOmics/survcomp are unavailable from a
## CRAN-only repository. Install those dependencies explicitly first.
if (!requireNamespace("mixOmics", quietly=TRUE) || !requireNamespace("survcomp", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) {
    cat("Installing BiocManager for plsRcox dependencies\n")
    try(install.packages("BiocManager", repos="https://cloud.r-project.org"), silent=TRUE)
  }
  if (requireNamespace("BiocManager", quietly=TRUE)) {
    if (!requireNamespace("mixOmics", quietly=TRUE)) {
      cat("Installing Bioconductor dependency: mixOmics\n")
      try(BiocManager::install("mixOmics", ask=FALSE, update=FALSE), silent=TRUE)
    }
    if (!requireNamespace("survcomp", quietly=TRUE)) {
      cat("Installing Bioconductor dependency: survcomp\n")
      try(BiocManager::install("survcomp", ask=FALSE, update=FALSE), silent=TRUE)
    }
  }
}
if (!requireNamespace("plsRcox", quietly=TRUE)) {
  cat("Installing missing package: plsRcox\n")
  try(install.packages("plsRcox", repos="https://cloud.r-project.org"), silent=TRUE)
}
if (!all(vapply(required, requireNamespace, logical(1), quietly=TRUE))) stop("Required packages survival/glmnet/quadprog are unavailable")
pkg_ok <- setNames(vapply(optional, requireNamespace, logical(1), quietly=TRUE), optional)
write.csv(data.frame(Package=c(required,optional),Available=unname(vapply(c(required,optional),requireNamespace,logical(1),quietly=TRUE))), file.path(out_dir,"Tables","00_package_status.csv"), row.names=FALSE)

## ---------------------------- selector ----------------------------
select_genes <- function(method, genes, train_z, time, event, seed=1) {
  set.seed(seed); genes <- unique(genes[genes %in% rownames(train_z)])
  if (length(genes) < 2) return(character())
  x <- t(train_z[genes,,drop=FALSE]); y <- survival::Surv(time,event)
  out <- tryCatch({
    if (method == "All") return(genes)
    if (method == "UnivCox") {
      u <- univ_rank(train_z,genes,time,event); if (!nrow(u)) return(character())
      s <- u$Gene[u$P < .05]; if (length(s)<5) s <- head(u$Gene,min(20,nrow(u))); return(unique(s))
    }
    if (method %in% c("Lasso","Enet","RidgeTop")) {
      a <- if (method=="Lasso") 1 else if (method=="Enet") .5 else 0
      cv <- glmnet::cv.glmnet(x,y,family="cox",alpha=a,nfolds=10,standardize=FALSE,maxit=100000)
      cf <- as.matrix(coef(cv,s="lambda.min")); d <- data.frame(Gene=rownames(cf),Coef=as.numeric(cf[,1]))
      if (method=="RidgeTop") return(head(d$Gene[order(-abs(d$Coef))],min(20,nrow(d))))
      return(unique(d$Gene[d$Coef != 0]))
    }
    if (method == "StepCox") {
      u <- univ_rank(train_z,genes,time,event); if (!nrow(u)) return(character())
      gs <- head(u$Gene,min(20,nrow(u))); nc <- make.names(gs,unique=TRUE)
      d <- as.data.frame(t(train_z[gs,,drop=FALSE])); colnames(d)<-nc; d$time<-time; d$event<-event
      ft <- survival::coxph(as.formula(paste0("survival::Surv(time,event)~",paste(nc,collapse="+"))),data=d)
      st <- step(ft,direction="both",trace=0); return(unique(gs[match(names(coef(st)),nc)]))
    }
    if (method == "CoxBoost") {
      if (!pkg_ok["CoxBoost"]) return(character())
      pen <- 9 * sum(event == 1, na.rm = TRUE)
      if (!is.finite(pen) || pen <= 0) return(character())
      cv <- CoxBoost::cv.CoxBoost(time=time,status=event,x=x,maxstepno=100,K=10,type="verweij",penalty=pen)
      opt_step <- as.integer(cv$optimal.step[1])
      if (!length(opt_step) || !is.finite(opt_step) || opt_step < 0) return(character())
      ft <- CoxBoost::CoxBoost(time=time,status=event,x=x,stepno=opt_step,penalty=pen)
      cf <- coef(ft,at.step=opt_step); return(names(cf)[is.finite(cf) & cf!=0])
    }
    if (method == "RSF") {
      if (!pkg_ok["randomForestSRC"]) return(character())
      d<-as.data.frame(x); d$time<-time; d$event<-event
      rsf_formula <- stats::as.formula("Surv(time,event) ~ .", env=asNamespace("survival"))
      ft<-randomForestSRC::rfsrc(rsf_formula,data=d,ntree=500,nodesize=10,importance=TRUE)
      im<-sort(ft$importance,decreasing=TRUE); return(head(names(im),min(20,length(im))))
    }
    if (method == "GBM") {
      if (!pkg_ok["gbm"]) return(character())
      d<-as.data.frame(x); d$time<-time; d$event<-event
      ft<-gbm::gbm(survival::Surv(time,event)~.,data=d,distribution="coxph",n.trees=2000,interaction.depth=2,shrinkage=.01,n.minobsinnode=10,bag.fraction=.7,train.fraction=1,verbose=FALSE)
      im<-summary(ft,plotit=FALSE); return(head(as.character(im$var),min(20,nrow(im))))
    }
    if (method == "SuperPC") {
      if (!pkg_ok["superpc"]) return(character())
      dat<-list(x=t(x),y=time,censoring.status=event,featurenames=genes)
      ft<-superpc::superpc.train(dat,type="survival"); fs<-ft$feature.scores; if(is.matrix(fs)) fs<-fs[,1]
      return(head(names(sort(abs(fs),decreasing=TRUE)),min(20,length(fs))))
    }
    character()
  }, error=function(e){cat("Selector failed",method,":",conditionMessage(e),"\n");character()})
  unique(out[!is.na(out) & out %in% rownames(train_z)])
}

## ---------------------------- learners ----------------------------
fit_learner <- function(learner, genes, train_z, valid_z, time, event,
                        valid_time=NULL, valid_event=NULL, model_name="model", seed=1) {
  set.seed(seed); genes<-unique(genes[genes %in% rownames(train_z) & genes %in% rownames(valid_z)])
  if (!length(genes)) return(NULL)
  x<-t(train_z[genes,,drop=FALSE]); xv<-t(valid_z[genes,,drop=FALSE])
  if (is.null(valid_time)) valid_time<-rep(1,nrow(xv)); if(is.null(valid_event)) valid_event<-rep(0,nrow(xv))
  out <- tryCatch({
    if (learner == "Cox") {
      nc<-make.names(genes,unique=TRUE); d<-as.data.frame(x); colnames(d)<-nc; d$time<-time; d$event<-event
      dv<-as.data.frame(xv); colnames(dv)<-nc
      ft<-survival::coxph(as.formula(paste0("survival::Surv(time,event)~",paste(nc,collapse="+"))),data=d)
      rt<-as.numeric(predict(ft,newdata=d,type="lp")); rv<-as.numeric(predict(ft,newdata=dv,type="lp")); cf<-coef(ft)
      cd<-data.frame(Model=model_name,Gene=genes[match(names(cf),nc)],Coef=as.numeric(cf)); return(finalize_fit_output(list(rt=rt,rv=rv,genes=cd$Gene,coef=cd),x,xv,time,event))
    }
    if (learner == "StepCox") {
      if(length(genes)<2) return(NULL); nc<-make.names(genes,unique=TRUE); d<-as.data.frame(x); colnames(d)<-nc; d$time<-time; d$event<-event
      dv<-as.data.frame(xv); colnames(dv)<-nc; ft<-survival::coxph(as.formula(paste0("survival::Surv(time,event)~",paste(nc,collapse="+"))),data=d)
      st<-step(ft,direction="both",trace=0); rt<-as.numeric(predict(st,newdata=d,type="lp")); rv<-as.numeric(predict(st,newdata=dv,type="lp")); cf<-coef(st)
      fg<-genes[match(names(cf),nc)]; cd<-data.frame(Model=model_name,Gene=fg,Coef=as.numeric(cf)); return(finalize_fit_output(list(rt=rt,rv=rv,genes=fg,coef=cd),x,xv,time,event))
    }
    if (learner %in% c("Lasso","Enet","Ridge")) {
      a<-if(learner=="Lasso")1 else if(learner=="Enet").5 else 0
      cv<-glmnet::cv.glmnet(x,survival::Surv(time,event),family="cox",alpha=a,nfolds=10,standardize=FALSE,maxit=100000)
      rt<-as.numeric(predict(cv,newx=x,s="lambda.min",type="link")); rv<-as.numeric(predict(cv,newx=xv,s="lambda.min",type="link"))
      cf<-as.matrix(coef(cv,s="lambda.min")); cd<-data.frame(Model=model_name,Gene=rownames(cf),Coef=as.numeric(cf[,1]))
      if(learner!="Ridge") cd<-cd[cd$Coef!=0,,drop=FALSE] else cd<-head(cd[order(-abs(cd$Coef)),,drop=FALSE],min(20,nrow(cd)))
      return(finalize_fit_output(list(rt=rt,rv=rv,genes=cd$Gene,coef=cd),x,xv,time,event))
    }
    if (learner == "CoxBoost") {
      if(!pkg_ok["CoxBoost"]) return(NULL); pen<-9*sum(event==1,na.rm=TRUE); if(!is.finite(pen)||pen<=0)return(NULL)
      cv<-CoxBoost::cv.CoxBoost(time=time,status=event,x=x,maxstepno=100,K=10,type="verweij",penalty=pen); opt_step<-as.integer(cv$optimal.step[1]); if(!length(opt_step)||!is.finite(opt_step)||opt_step<0)return(NULL)
      ft<-CoxBoost::CoxBoost(time=time,status=event,x=x,stepno=opt_step,penalty=pen)
      rt<-as.numeric(predict(ft,newdata=x,type="lp",at.step=opt_step)); rv<-as.numeric(predict(ft,newdata=xv,type="lp",at.step=opt_step)); cf<-coef(ft,at.step=opt_step)
      nz<-is.finite(cf)&cf!=0; cd<-data.frame(Model=model_name,Gene=names(cf)[nz],Coef=as.numeric(cf[nz])); return(finalize_fit_output(list(rt=rt,rv=rv,genes=cd$Gene,coef=cd),x,xv,time,event))
    }
    if (learner == "RSF") {
      if(!pkg_ok["randomForestSRC"]) return(NULL); d<-as.data.frame(x); d$time<-time; d$event<-event; dv<-as.data.frame(xv)
      rsf_formula<-stats::as.formula("Surv(time,event) ~ .",env=asNamespace("survival")); ft<-randomForestSRC::rfsrc(rsf_formula,data=d,ntree=1000,nodesize=10,importance=TRUE)
      rt<-as.numeric(predict(ft,newdata=d)$predicted); rv<-as.numeric(predict(ft,newdata=dv)$predicted); im<-sort(ft$importance,decreasing=TRUE); fg<-head(names(im),min(20,length(im)))
      cd<-data.frame(Model=model_name,Gene=fg,Coef=as.numeric(im[fg])); return(finalize_fit_output(list(rt=rt,rv=rv,genes=fg,coef=cd),x,xv,time,event))
    }
    if (learner == "GBM") {
      if(!pkg_ok["gbm"]) return(NULL); d<-as.data.frame(x); d$time<-time; d$event<-event; dv<-as.data.frame(xv)
      ft<-gbm::gbm(survival::Surv(time,event)~.,data=d,distribution="coxph",n.trees=3000,interaction.depth=2,shrinkage=.01,n.minobsinnode=10,bag.fraction=.7,train.fraction=1,verbose=FALSE)
      rt<-as.numeric(predict(ft,newdata=d,n.trees=3000,type="link")); rv<-as.numeric(predict(ft,newdata=dv,n.trees=3000,type="link")); im<-summary(ft,plotit=FALSE); fg<-head(as.character(im$var),min(20,nrow(im)))
      cd<-data.frame(Model=model_name,Gene=fg,Coef=im$rel.inf[match(fg,im$var)]); return(finalize_fit_output(list(rt=rt,rv=rv,genes=fg,coef=cd),x,xv,time,event))
    }
    if (learner == "SuperPC") {
      if(!pkg_ok["superpc"]) return(NULL); dt<-list(x=t(x),y=time,censoring.status=event,featurenames=genes); dv<-list(x=t(xv),y=valid_time,censoring.status=valid_event,featurenames=genes)
      ft<-superpc::superpc.train(dt,type="survival"); fs<-ft$feature.scores; if(is.matrix(fs)) fs<-fs[,1]; th<-as.numeric(quantile(abs(fs),.75,na.rm=TRUE))
      pt<-superpc::superpc.predict(ft,dt,dt,threshold=th,n.components=1,prediction.type="continuous"); pv<-superpc::superpc.predict(ft,dt,dv,threshold=th,n.components=1,prediction.type="continuous")
      fg<-head(names(sort(abs(fs),decreasing=TRUE)),min(20,length(fs))); cd<-data.frame(Model=model_name,Gene=fg,Coef=as.numeric(fs[fg])); return(finalize_fit_output(list(rt=as.numeric(pt$v.pred),rv=as.numeric(pv$v.pred),genes=fg,coef=cd),x,xv,time,event))
    }
    if (learner == "SurvivalSVM") {
      if(!pkg_ok["survivalsvm"]) return(NULL); d<-as.data.frame(x); d$time<-time; d$event<-event; dv<-as.data.frame(xv); dv$time<-valid_time; dv$event<-valid_event
      ft<-survivalsvm::survivalsvm(survival::Surv(time,event)~.,data=d,type="regression",gamma.mu=1,opt.meth="quadprog",kernel="lin_kernel")
      rt<-as.numeric(predict(ft,newdata=d)$predicted); rv<-as.numeric(predict(ft,newdata=dv)$predicted); cd<-data.frame(Model=model_name,Gene=genes,Coef=NA_real_); return(finalize_fit_output(list(rt=rt,rv=rv,genes=genes,coef=cd),x,xv,time,event))
    }
    NULL
  }, error=function(e){cat("Learner failed",model_name,":",conditionMessage(e),"\n");NULL})
  if (is.null(out) || length(unique(out$rt[is.finite(out$rt)]))<2 || length(unique(out$rv[is.finite(out$rv)]))<2) return(NULL)
  oo<-orient_scores(out$rt,out$rv,time,event); out$rt<-oo$train; out$rv<-oo$valid
  names(out$rt)<-rownames(x); names(out$rv)<-rownames(xv); out
}

fit_plsr <- function(genes,train_z,valid_z,time,event,model_name,seed=1) {
  if(!pkg_ok["plsRcox"]) return(NULL); set.seed(seed); genes<-genes[genes%in%rownames(train_z)&genes%in%rownames(valid_z)]
  if(length(genes)<2)return(NULL); x<-t(train_z[genes,,drop=FALSE]); xv<-t(valid_z[genes,,drop=FALSE])
  out<-tryCatch({ft<-plsRcox::plsRcox(Xplan=x,time=time,event=event,nt=min(3,ncol(x))); rt<-as.numeric(predict(ft,newdata=x,type="lp")); rv<-as.numeric(predict(ft,newdata=xv,type="lp")); list(rt=rt,rv=rv,genes=genes,coef=data.frame(Model=model_name,Gene=genes,Coef=NA_real_))},error=function(e)NULL)
  if(is.null(out))return(NULL); oo<-orient_scores(out$rt,out$rv,time,event); out$rt<-oo$train;out$rv<-oo$valid
  risk_mu<-mean(out$rt,na.rm=TRUE);risk_sd<-sd(out$rt,na.rm=TRUE);if(!is.finite(risk_mu)||!is.finite(risk_sd)||risk_sd<=0)return(NULL)
  out$rt<-(out$rt-risk_mu)/risk_sd;out$rv<-(out$rv-risk_mu)/risk_sd
  names(out$rt)<-rownames(x);names(out$rv)<-rownames(xv);out
}

## ---------------------------- metric ----------------------------
calc_metrics <- function(time,event,score) {
  d<-data.frame(time=safe_num(time),event=safe_num(event),score=safe_num(score)); d<-d[complete.cases(d)&d$time>0,,drop=FALSE]
  if(nrow(d)<30||sum(d$event==1)<5||length(unique(d$score))<2)return(data.frame(N=nrow(d),Events=sum(d$event==1),HR=NA,Lower95=NA,Upper95=NA,CoxP=NA,Cindex=NA,KMHR=NA,LogRankP=NA,AUC1=NA,AUC3=NA,AUC5=NA,MeanAUC=NA))
  d$z<-as.numeric(scale(d$score)); ft<-survival::coxph(survival::Surv(time,event)~z,data=d); s<-summary(ft)
  ci<-tryCatch(as.numeric(survival::concordance(survival::Surv(time,event)~score,data=d,reverse=TRUE)$concordance),error=function(e)NA_real_)
  med<-median(d$score); d$group<-factor(ifelse(d$score>med,"High","Low"),levels=c("Low","High")); km<-tryCatch(survival::coxph(survival::Surv(time,event)~group,data=d),error=function(e)NULL); sd<-tryCatch(survival::survdiff(survival::Surv(time,event)~group,data=d),error=function(e)NULL)
  kmhr<-if(is.null(km))NA_real_ else exp(coef(km)[1]); lrp<-if(is.null(sd))NA_real_ else 1-pchisq(sd$chisq,length(sd$n)-1)
  auc<-rep(NA_real_,3); if(pkg_ok["timeROC"]) { rr<-tryCatch(timeROC::timeROC(T=d$time,delta=d$event,marker=d$score,cause=1,weighting="marginal",times=c(365,1095,1825),iid=FALSE),error=function(e)NULL); if(!is.null(rr))auc<-as.numeric(rr$AUC) }
  data.frame(N=nrow(d),Events=sum(d$event==1),HR=s$coefficients[1,"exp(coef)"],Lower95=s$conf.int[1,"lower .95"],Upper95=s$conf.int[1,"upper .95"],CoxP=s$coefficients[1,"Pr(>|z|)"],Cindex=ci,KMHR=kmhr,LogRankP=lrp,AUC1=auc[1],AUC3=auc[2],AUC5=auc[3],MeanAUC=if(all(is.na(auc)))NA_real_ else mean(auc,na.rm=TRUE))
}

## ---------------------------- inputs ----------------------------
for(f in c(emtab_expr_file,emtab_clin_file)) if(!file.exists(f)) stop("Missing input: ",f)
emtab_expr<-collapse_duplicate_genes(as_expr_matrix(readRDS(emtab_expr_file)))
emtab_clin<-read.csv(emtab_clin_file,check.names=FALSE)
if(!"Sample"%in%colnames(emtab_clin))stop("E-MTAB clinical lacks Sample")
cm<-intersect(colnames(emtab_expr),emtab_clin$Sample); emtab_expr<-emtab_expr[,cm,drop=FALSE];emtab_clin<-emtab_clin[match(cm,emtab_clin$Sample),,drop=FALSE]
emtab_time<-safe_num(emtab_clin$OVERALLSURVIVAL); emtab_event<-clean_event(emtab_clin$DEATHOFDISEASE); keep<-is.finite(emtab_time)&emtab_time>0&!is.na(emtab_event)
emtab_expr<-emtab_expr[,keep,drop=FALSE];emtab_clin<-emtab_clin[keep,,drop=FALSE];emtab_time<-emtab_time[keep];emtab_event<-emtab_event[keep];emtab_samples<-colnames(emtab_expr)

## Frozen before this audit from SC28 RECALC 04C Table 01.
## These 24 genes were the intersection of the audited BALANCED_40 panel with
## genes measurable in both E-MTAB-179 and SEQC.  SEQC is NOT opened here.
candidate_genes<-c(
  "PRDX4","SLC16A1","LDHA","LDHB","PGK1","ENO1","PDK1","HK2",
  "GPI","ALDOA","PRDX6","PRDX2","TXN","PRDX5","TXNRD1","NQO1",
  "NUDT5","ELK1","SLC25A5","SNAPC1","NUP37","RUVBL1","SLC1A5","AHCY"
)
missing_emtab<-setdiff(candidate_genes,rownames(emtab_expr))
if(length(missing_emtab))stop("Frozen audit genes missing from E-MTAB expression: ",paste(missing_emtab,collapse=", "))
if(length(candidate_genes)!=24L||anyDuplicated(candidate_genes))stop("Frozen candidate panel integrity check failed")
write.csv(data.frame(Gene=candidate_genes,FrozenForRepeatedCV=TRUE,In_EMTAB=TRUE),file.path(out_dir,"Tables","01_frozen_24gene_panel.csv"),row.names=FALSE)
cat("Frozen candidate genes:",length(candidate_genes)," PRDX4:","PRDX4"%in%candidate_genes," LDHA:","LDHA"%in%candidate_genes,"\n")

## Historical 101 numbering.
selectors<-c("All","UnivCox","Lasso","Enet","RidgeTop","StepCox","CoxBoost","RSF","GBM","SuperPC")
learners<-c("Cox","StepCox","Lasso","Enet","Ridge","CoxBoost","RSF","GBM","SuperPC","SurvivalSVM")
combos<-expand.grid(Selector=selectors,Learner=learners,stringsAsFactors=FALSE);combos$ComboID<-seq_len(nrow(combos));combos<-rbind(combos,data.frame(Selector="All",Learner="PLSRcox_extra",ComboID=101));combos$Model<-paste0("M",sprintf("%03d",combos$ComboID),"_",combos$Selector,"_",combos$Learner)
write.csv(combos,file.path(out_dir,"Tables","02_standard_101_plan_preserved.csv"),row.names=FALSE)
if(combos$Model[50]!="M050_SuperPC_Ridge")stop("Historical combination numbering mismatch")

## ---------------------------- repeated 5-fold OOF ----------------------------
## Pre-create every fold assignment before any model is fit.  Each patient is
## validation exactly once per repeat; events are stratified across the 5 folds.
fold_assignments<-lapply(seq_len(repeats_n),function(rep_i){
  fid<-make_folds(emtab_event,k_folds,seed_main+rep_i*1000L)
  data.frame(Repeat=rep_i,Sample=emtab_samples,OS_time=emtab_time,OS_event=emtab_event,Fold=fid)
})
fold_assignment_df<-do.call(rbind,fold_assignments)
write.csv(fold_assignment_df,file.path(out_dir,"Tables","03_EMTAB_repeated5fold_assignments.csv"),row.names=FALSE)

for(rep_i in seq_len(repeats_n)) {
  fold_id<-fold_assignments[[rep_i]]$Fold
  cat("\n================ REPEAT",rep_i,"/",repeats_n,"================\n")
  for(fold in seq_len(k_folds)) {
    ck<-file.path(out_dir,"Checkpoints",sprintf("repeat_%02d_fold_%d_results.rds",rep_i,fold))
    if(resume_existing_folds && file.exists(ck)){
      cat("Repeat",rep_i,"fold",fold,"checkpoint exists; skipping.\n")
      next
    }
    cat("\n===== REPEAT",rep_i,"FOLD",fold,"/",k_folds,"=====\n")
    vi<-which(fold_id==fold);ti<-which(fold_id!=fold)
    sc<-scale_train_apply(emtab_expr[candidate_genes,ti,drop=FALSE],emtab_expr[candidate_genes,vi,drop=FALSE])
    tz<-sc$train_z;vz<-sc$valid_z

    sel_cache<-list()
    for(si in seq_along(selectors)){
      sm<-selectors[si]
      cat(" selector",sm,"\n")
      sel_cache[[sm]]<-select_genes(sm,candidate_genes,tz,emtab_time[ti],emtab_event[ti],seed_main+rep_i*100000L+fold*100L+si)
    }

    risks<-list();mans<-list();fails<-list()
    for(i in seq_len(nrow(combos))){
      mo<-combos$Model[i];sm<-combos$Selector[i];lr<-combos$Learner[i];gs<-sel_cache[[sm]]
      cat(" ",mo,"\n")
      if(!length(gs)){
        fails[[mo]]<-data.frame(Repeat=rep_i,Fold=fold,Model=mo,Reason="selector returned no genes")
        next
      }
      fit_seed<-seed_main+rep_i*1000000L+fold*10000L+i
      ## Validation outcomes are deliberately NOT passed to prediction functions.
      fit<-if(lr=="PLSRcox_extra") {
        fit_plsr(gs,tz,vz,emtab_time[ti],emtab_event[ti],mo,fit_seed)
      } else {
        fit_learner(lr,gs,tz,vz,emtab_time[ti],emtab_event[ti],NULL,NULL,mo,fit_seed)
      }
      if(is.null(fit)){
        fails[[mo]]<-data.frame(Repeat=rep_i,Fold=fold,Model=mo,Reason="fit/prediction failed")
        next
      }
      rv<-as.numeric(fit$rv);valid_ids<-colnames(vz)
      if(length(rv)!=length(valid_ids)){
        fails[[mo]]<-data.frame(Repeat=rep_i,Fold=fold,Model=mo,Reason=paste0("prediction length mismatch: ",length(rv)," vs ",length(valid_ids)))
        next
      }
      risks[[mo]]<-data.frame(Repeat=rep_i,Fold=fold,Model=mo,Sample=valid_ids,RiskScore=rv)
      mans[[mo]]<-data.frame(
        Repeat=rep_i,Fold=fold,Model=mo,Selector=sm,Learner=lr,
        SelectedGeneN=length(gs),FinalGeneN=length(unique(fit$genes)),
        SelectedContainsPRDX4="PRDX4"%in%gs,FinalContainsPRDX4="PRDX4"%in%fit$genes,
        SelectedGenes=paste(gs,collapse=";"),FinalGenes=paste(unique(fit$genes),collapse=";")
      )
    }
    bind<-function(x)if(length(x))do.call(rbind,x)else data.frame()
    saveRDS(list(risk=bind(risks),manifest=bind(mans),fail=bind(fails)),ck)
    cat("Repeat",rep_i,"fold",fold,"saved:",ck,"\n")
  }
}

## ---------------------------- collect checkpoints ----------------------------
ck_grid<-expand.grid(Repeat=seq_len(repeats_n),Fold=seq_len(k_folds))
ck_files<-mapply(function(r,f)file.path(out_dir,"Checkpoints",sprintf("repeat_%02d_fold_%d_results.rds",r,f)),ck_grid$Repeat,ck_grid$Fold,USE.NAMES=FALSE)
if(any(!file.exists(ck_files)))stop("Missing repeat/fold checkpoints: ",paste(basename(ck_files[!file.exists(ck_files)]),collapse=", "))
fold_objs<-lapply(ck_files,readRDS)
bind_nonempty<-function(x){
  x<-x[vapply(x,function(z)is.data.frame(z)&&nrow(z)>0,logical(1))]
  if(!length(x))return(data.frame())
  do.call(rbind,x)
}
oof_risk<-bind_nonempty(lapply(fold_objs,`[[`,"risk"))
oof_manifest<-bind_nonempty(lapply(fold_objs,`[[`,"manifest"))
oof_fail<-bind_nonempty(lapply(fold_objs,`[[`,"fail"))
if(!nrow(oof_risk)||!nrow(oof_manifest))stop("No usable repeated-CV predictions were produced")
write.csv(oof_risk,file.path(out_dir,"Tables","04_all_repeated_OOF_risk_scores.csv"),row.names=FALSE)
write.csv(oof_manifest,file.path(out_dir,"Tables","05_repeat_fold_model_manifest.csv"),row.names=FALSE)
write.csv(oof_fail,file.path(out_dir,"Tables","06_repeat_fold_failures.csv"),row.names=FALSE)

## ---------------------------- metrics for every repeat x model ----------------------------
repeat_metric_list<-list();zz<-1L
for(rep_i in seq_len(repeats_n)){
  for(m in combos$Model){
    z<-oof_risk[oof_risk$Repeat==rep_i & oof_risk$Model==m,,drop=FALSE]
    z<-z[!duplicated(z$Sample),,drop=FALSE]
    ix<-match(z$Sample,emtab_samples)
    me<-calc_metrics(emtab_time[ix],emtab_event[ix],z$RiskScore)
    mn<-mean(oof_manifest$FinalGeneN[oof_manifest$Repeat==rep_i & oof_manifest$Model==m],na.rm=TRUE)
    repeat_metric_list[[zz]]<-cbind(Repeat=rep_i,Model=m,CoverageN=nrow(z),CoverageFrac=nrow(z)/length(emtab_samples),MeanFinalGeneN=mn,me)
    zz<-zz+1L
  }
}
repeat_metrics<-do.call(rbind,repeat_metric_list)
numcols<-setdiff(colnames(repeat_metrics),"Model")
repeat_metrics[numcols]<-lapply(repeat_metrics[numcols],safe_num)
repeat_metrics$Eligible<-repeat_metrics$CoverageFrac>=.95 & is.finite(repeat_metrics$Cindex)
repeat_metrics$CindexRank<-NA_real_;repeat_metrics$MeanAUCRank<-NA_real_;repeat_metrics$RankSum<-NA_real_
for(rep_i in seq_len(repeats_n)){
  ii<-which(repeat_metrics$Repeat==rep_i & repeat_metrics$Eligible)
  repeat_metrics$CindexRank[ii]<-rank(-repeat_metrics$Cindex[ii],ties.method="min")
  if(any(is.finite(repeat_metrics$MeanAUC[ii]))){
    aa<-ii[is.finite(repeat_metrics$MeanAUC[ii])]
    repeat_metrics$MeanAUCRank[aa]<-rank(-repeat_metrics$MeanAUC[aa],ties.method="min")
  }
  repeat_metrics$RankSum[ii]<-repeat_metrics$CindexRank[ii]+repeat_metrics$MeanAUCRank[ii]
}
write.csv(repeat_metrics,file.path(out_dir,"Tables","07_repeat_level_101_model_metrics.csv"),row.names=FALSE)

## ---------------------------- model-level stability summary ----------------------------
safe_mean<-function(x)if(any(is.finite(x)))mean(x[is.finite(x)])else NA_real_
safe_sd<-function(x)if(sum(is.finite(x))>=2)sd(x[is.finite(x)])else NA_real_
safe_median<-function(x)if(any(is.finite(x)))median(x[is.finite(x)])else NA_real_
safe_min<-function(x)if(any(is.finite(x)))min(x[is.finite(x)])else NA_real_
safe_max<-function(x)if(any(is.finite(x)))max(x[is.finite(x)])else NA_real_

model_summary<-do.call(rbind,lapply(combos$Model,function(m){
  z<-repeat_metrics[repeat_metrics$Model==m,,drop=FALSE]
  data.frame(
    Model=m,RepeatsN=nrow(z),CompleteRepeats=sum(z$Eligible),
    MeanCindex=safe_mean(z$Cindex),SDCindex=safe_sd(z$Cindex),SECindex=safe_sd(z$Cindex)/sqrt(sum(is.finite(z$Cindex))),
    MeanMeanAUC=safe_mean(z$MeanAUC),SDMeanAUC=safe_sd(z$MeanAUC),SEMeanAUC=safe_sd(z$MeanAUC)/sqrt(sum(is.finite(z$MeanAUC))),
    MeanCindexRank=safe_mean(z$CindexRank),MedianCindexRank=safe_median(z$CindexRank),BestCindexRank=safe_min(z$CindexRank),WorstCindexRank=safe_max(z$CindexRank),
    Top1_Cindex=sum(z$CindexRank<=1,na.rm=TRUE),Top5_Cindex=sum(z$CindexRank<=5,na.rm=TRUE),Top10_Cindex=sum(z$CindexRank<=10,na.rm=TRUE),
    MeanAUCRank=safe_mean(z$MeanAUCRank),MedianAUCRank=safe_median(z$MeanAUCRank),Top5_AUC=sum(z$MeanAUCRank<=5,na.rm=TRUE),Top10_AUC=sum(z$MeanAUCRank<=10,na.rm=TRUE),
    MeanFinalGeneN=safe_mean(z$MeanFinalGeneN),stringsAsFactors=FALSE
  )
}))

mean_jaccard<-function(gene_strings){
  sets<-lapply(gene_strings,function(s){g<-unique(strsplit(as.character(s),";",fixed=TRUE)[[1]]);g[!is.na(g)&g!=""]})
  if(length(sets)<2)return(NA_real_)
  cmb<-combn(seq_along(sets),2)
  jj<-apply(cmb,2,function(k){u<-union(sets[[k[1]]],sets[[k[2]]]);if(!length(u))return(NA_real_);length(intersect(sets[[k[1]]],sets[[k[2]]]))/length(u)})
  safe_mean(jj)
}
feature_stability<-do.call(rbind,lapply(combos$Model,function(m){
  z<-oof_manifest[oof_manifest$Model==m,,drop=FALSE]
  data.frame(
    Model=m,FoldFits=nrow(z),PRDX4_FoldCount=sum(z$FinalContainsPRDX4,na.rm=TRUE),
    PRDX4_FoldFraction=mean(z$FinalContainsPRDX4,na.rm=TRUE),
    MeanGeneJaccard=mean_jaccard(z$FinalGenes),stringsAsFactors=FALSE
  )
}))
model_summary<-merge(model_summary,feature_stability,by="Model",all.x=TRUE,sort=FALSE)
model_summary$MeanCindexSummaryRank<-NA_real_;model_summary$MeanAUCSummaryRank<-NA_real_
ii_c<-which(is.finite(model_summary$MeanCindex));ii_a<-which(is.finite(model_summary$MeanMeanAUC))
model_summary$MeanCindexSummaryRank[ii_c]<-rank(-model_summary$MeanCindex[ii_c],ties.method="min")
model_summary$MeanAUCSummaryRank[ii_a]<-rank(-model_summary$MeanMeanAUC[ii_a],ties.method="min")
model_summary$DescriptiveRankSum<-model_summary$MeanCindexSummaryRank+model_summary$MeanAUCSummaryRank
model_summary<-model_summary[order(model_summary$DescriptiveRankSum,-model_summary$MeanCindex,-model_summary$MeanMeanAUC),,drop=FALSE]
model_summary$DescriptiveRank<-seq_len(nrow(model_summary))

finite_c_idx<-which(is.finite(model_summary$MeanCindex))
if(!length(finite_c_idx))stop("No model has a finite repeated-CV mean C-index")
best_c_idx<-finite_c_idx[which.max(model_summary$MeanCindex[finite_c_idx])]
best_c_model<-model_summary$Model[best_c_idx]
best_c_mean<-model_summary$MeanCindex[best_c_idx]
best_c_se<-model_summary$SECindex[best_c_idx]
one_se_threshold<-best_c_mean-best_c_se
model_summary$Within1SEofBestCindex<-is.finite(model_summary$MeanCindex)&model_summary$MeanCindex>=one_se_threshold
write.csv(model_summary,file.path(out_dir,"Tables","08_model_stability_summary.csv"),row.names=FALSE)
write.csv(head(model_summary,20),file.path(out_dir,"Tables","09_descriptive_top20_stability.csv"),row.names=FALSE)

## Gene-level selection frequencies across all repeat/fold fits.
gene_freq_rows<-list();zz<-1L
for(m in combos$Model){
  z<-oof_manifest[oof_manifest$Model==m,,drop=FALSE]
  gs<-lapply(z$FinalGenes,function(s)unique(strsplit(as.character(s),";",fixed=TRUE)[[1]]))
  tab<-sort(table(unlist(gs)),decreasing=TRUE)
  if(length(tab)){
    gene_freq_rows[[zz]]<-data.frame(Model=m,Gene=names(tab),FoldCount=as.integer(tab),FoldFraction=as.integer(tab)/nrow(z),stringsAsFactors=FALSE)
    zz<-zz+1L
  }
}
gene_frequency<-if(length(gene_freq_rows))do.call(rbind,gene_freq_rows)else data.frame()
write.csv(gene_frequency,file.path(out_dir,"Tables","10_gene_selection_frequency.csv"),row.names=FALSE)

## Paired repeat-level comparison of the best mean-C-index model with M050.
m050_id<-"M050_SuperPC_Ridge"
best_rep<-repeat_metrics[repeat_metrics$Model==best_c_model,c("Repeat","Cindex","MeanAUC"),drop=FALSE]
m050_rep<-repeat_metrics[repeat_metrics$Model==m050_id,c("Repeat","Cindex","MeanAUC"),drop=FALSE]
colnames(best_rep)[2:3]<-c("Best_Cindex","Best_MeanAUC");colnames(m050_rep)[2:3]<-c("M050_Cindex","M050_MeanAUC")
paired<-merge(best_rep,m050_rep,by="Repeat",all=FALSE)
paired$DeltaCindex_BestMinusM050<-paired$Best_Cindex-paired$M050_Cindex
paired$DeltaMeanAUC_BestMinusM050<-paired$Best_MeanAUC-paired$M050_MeanAUC
paired$BestModel<-best_c_model
write.csv(paired,file.path(out_dir,"Tables","11_bestMeanC_vs_M050_by_repeat.csv"),row.names=FALSE)

m050_sum<-model_summary[model_summary$Model==m050_id,,drop=FALSE]
if(!nrow(m050_sum))stop("M050 missing from stability summary")
final<-data.frame(
  Item=c(
    "E-MTAB sample N","E-MTAB event N","Frozen candidate gene N","Repeated CV repeats","Folds per repeat","Total outer validation folds",
    "Historical model combinations","Repeat/fold failure rows","Best model by mean repeated-CV C-index","Best mean C-index","One-SE C-index threshold",
    "M050 mean C-index","M050 SD C-index","M050 mean MeanAUC","M050 mean C-index rank","M050 Top-5 C-index count","M050 PRDX4 fold fraction",
    "M050 mean gene Jaccard","M050 within 1SE of best C-index","SEQC expression opened in this audit","SEQC survival outcomes opened in this audit"
  ),
  Value=c(
    length(emtab_time),sum(emtab_event==1),length(candidate_genes),repeats_n,k_folds,repeats_n*k_folds,
    nrow(combos),nrow(oof_fail),best_c_model,best_c_mean,one_se_threshold,
    m050_sum$MeanCindex,m050_sum$SDCindex,m050_sum$MeanMeanAUC,m050_sum$MeanCindexRank,m050_sum$Top5_Cindex,m050_sum$PRDX4_FoldFraction,
    m050_sum$MeanGeneJaccard,m050_sum$Within1SEofBestCindex,"FALSE","FALSE"
  ),stringsAsFactors=FALSE
)
write.csv(final,file.path(out_dir,"Tables","12_FINAL_STABILITY_AUDIT_SUMMARY.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out_dir,"Logs","sessionInfo.txt"))

cat("\n============================================================\n")
cat("SC28 RECALC 05 completed.\n")
print(final)
cat("\nIMPORTANT: This script does NOT choose or lock a final model.\n")
cat("Upload the ENTIRE SC28_RECALC_05 output folder as ZIP for final audit.\n")
cat("End:",as.character(Sys.time()),"\n")
cat("============================================================\n")
sink()
