rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20260806)

## ============================================================
## PRDX4_RECALC_04_Clean101ML_EMTAB_CV_SEQC_External.R
##
## Clean model-development design
##   BALANCED_40 candidate panel
##   -> E-MTAB-179 only: stratified 5-fold out-of-fold model comparison
##   -> lock the best of the historical 101 combinations
##   -> refit that one pipeline on all E-MTAB-179
##   -> evaluate the locked model once in SEQC
##
## Important safeguards
##   1. Historical ComboID is preserved; M050 = SuperPC -> Ridge.
##   2. PRDX4/LDHA are NOT appended. The audited BALANCED_40 is used as saved.
##   3. SEQC survival outcomes are not read until the best pipeline is locked.
##   4. Fold checkpoints allow safe resume after interruption.
## ============================================================

root_dir <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"
seed_main <- 20260806
k_folds <- 5L
resume_existing_folds <- TRUE

emtab_expr_file <- file.path(root_dir, "SC13A_MANUAL_V6_EMTAB179_ADF_TabFix", "Tables", "05_EMTAB179_expr_gene_symbol_matched.rds")
emtab_clin_file <- file.path(root_dir, "SC13A_MANUAL_V5_EMTAB179_SemicolonFix", "Tables", "04_EMTAB179_clin_matched_raw.csv")
seqc_ready_file <- file.path(root_dir, "01_Preprocess", "SEQC_GSE49711_ready_FIX2", "SEQC_external_validation_ready_input_FIX2.rds")
sc12a_panel_rds <- file.path(root_dir, "SC12A_PRDX4_CandidatePool_Refinement_ModelInput", "Tables", "07_SC12A_SEQC_model_input_PRDX4_panels.rds")

out_dir <- file.path(root_dir, "SC28_RECALC_04C_Clean101ML_EMTAB_CV_SEQC_External_20260808")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Checkpoints"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "Logs", "SC28_RECALC_04_log.txt")
sink(log_file, split = TRUE)
on.exit(try(sink(), silent = TRUE), add = TRUE)

cat("============================================================\n")
cat("SC28 RECALC 04: Clean 101ML with E-MTAB-only model selection\n")
cat("Start:", as.character(Sys.time()), "\n")
cat("Seed:", seed_main, "  folds:", k_folds, "\n")
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
  if(is.null(out))return(NULL); oo<-orient_scores(out$rt,out$rv,time,event); out$rt<-oo$train;out$rv<-oo$valid;names(out$rt)<-rownames(x);names(out$rv)<-rownames(xv);out
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
for(f in c(emtab_expr_file,emtab_clin_file,seqc_ready_file,sc12a_panel_rds)) if(!file.exists(f)) stop("Missing input: ",f)
emtab_expr<-collapse_duplicate_genes(as_expr_matrix(readRDS(emtab_expr_file)))
emtab_clin<-read.csv(emtab_clin_file,check.names=FALSE)
if(!"Sample"%in%colnames(emtab_clin))stop("E-MTAB clinical lacks Sample")
cm<-intersect(colnames(emtab_expr),emtab_clin$Sample); emtab_expr<-emtab_expr[,cm,drop=FALSE];emtab_clin<-emtab_clin[match(cm,emtab_clin$Sample),,drop=FALSE]
emtab_time<-safe_num(emtab_clin$OVERALLSURVIVAL); emtab_event<-clean_event(emtab_clin$DEATHOFDISEASE); keep<-is.finite(emtab_time)&emtab_time>0&!is.na(emtab_event)
emtab_expr<-emtab_expr[,keep,drop=FALSE];emtab_clin<-emtab_clin[keep,,drop=FALSE];emtab_time<-emtab_time[keep];emtab_event<-emtab_event[keep];emtab_samples<-colnames(emtab_expr)

sc12<-readRDS(sc12a_panel_rds); if(is.null(sc12$panels$BALANCED_40))stop("BALANCED_40 missing from SC12A")
candidate_query<-unique(toupper(trimws(as.character(sc12$panels$BALANCED_40))));candidate_query<-candidate_query[!is.na(candidate_query)&candidate_query!=""]
cat("BALANCED_40 saved gene count:",length(candidate_query)," PRDX4:","PRDX4"%in%candidate_query," LDHA:","LDHA"%in%candidate_query,"\n")

## Load SEQC expression for platform harmonization only. DO NOT read outcomes here.
seqc_obj<-readRDS(seqc_ready_file); if(is.null(seqc_obj$expr)||is.null(seqc_obj$clin))stop("SEQC object lacks expr/clin")
seqc_expr<-collapse_duplicate_genes(as_expr_matrix(seqc_obj$expr))
candidate_genes<-candidate_query[candidate_query%in%rownames(emtab_expr)&candidate_query%in%rownames(seqc_expr)]
if(length(candidate_genes)<10)stop("Too few cross-platform BALANCED_40 genes: ",length(candidate_genes))
write.csv(data.frame(Gene=candidate_query,In_EMTAB=candidate_query%in%rownames(emtab_expr),In_SEQC=candidate_query%in%rownames(seqc_expr),In_Both=candidate_query%in%rownames(emtab_expr)&candidate_query%in%rownames(seqc_expr)),file.path(out_dir,"Tables","01_candidate_platform_coverage.csv"),row.names=FALSE)
cat("Cross-platform candidate genes:",length(candidate_genes),"\n")

## Historical 101 numbering.
selectors<-c("All","UnivCox","Lasso","Enet","RidgeTop","StepCox","CoxBoost","RSF","GBM","SuperPC")
learners<-c("Cox","StepCox","Lasso","Enet","Ridge","CoxBoost","RSF","GBM","SuperPC","SurvivalSVM")
combos<-expand.grid(Selector=selectors,Learner=learners,stringsAsFactors=FALSE);combos$ComboID<-seq_len(nrow(combos));combos<-rbind(combos,data.frame(Selector="All",Learner="PLSRcox_extra",ComboID=101));combos$Model<-paste0("M",sprintf("%03d",combos$ComboID),"_",combos$Selector,"_",combos$Learner)
write.csv(combos,file.path(out_dir,"Tables","02_standard_101_plan_preserved.csv"),row.names=FALSE)
if(combos$Model[50]!="M050_SuperPC_Ridge")stop("Historical combination numbering mismatch")

## ---------------------------- 5-fold OOF ----------------------------
fold_id<-make_folds(emtab_event,k_folds,seed_main)
write.csv(data.frame(Sample=emtab_samples,OS_time=emtab_time,OS_event=emtab_event,Fold=fold_id),file.path(out_dir,"Tables","03_EMTAB_5fold_assignment.csv"),row.names=FALSE)

for(fold in seq_len(k_folds)) {
  ck<-file.path(out_dir,"Checkpoints",paste0("fold_",fold,"_results.rds"))
  if(resume_existing_folds && file.exists(ck)){cat("Fold",fold,"checkpoint exists; skipping.\n");next}
  cat("\n===== FOLD",fold,"/",k_folds,"=====\n");vi<-which(fold_id==fold);ti<-which(fold_id!=fold)
  sc<-scale_train_apply(emtab_expr[candidate_genes,ti,drop=FALSE],emtab_expr[candidate_genes,vi,drop=FALSE]);tz<-sc$train_z;vz<-sc$valid_z
  sel_cache<-list();for(si in seq_along(selectors)){sm<-selectors[si];cat(" selector",sm,"\n");sel_cache[[sm]]<-select_genes(sm,candidate_genes,tz,emtab_time[ti],emtab_event[ti],seed_main+fold*100+si)}
  risks<-list();mans<-list();fails<-list()
  for(i in seq_len(nrow(combos))){mo<-combos$Model[i];sm<-combos$Selector[i];lr<-combos$Learner[i];cat(" ",mo,"\n");gs<-sel_cache[[sm]]
    if(!length(gs)){fails[[mo]]<-data.frame(Fold=fold,Model=mo,Reason="selector returned no genes");next}
    fit<-if(lr=="PLSRcox_extra")fit_plsr(gs,tz,vz,emtab_time[ti],emtab_event[ti],mo,seed_main+fold*10000+i) else fit_learner(lr,gs,tz,vz,emtab_time[ti],emtab_event[ti],emtab_time[vi],emtab_event[vi],mo,seed_main+fold*10000+i)
    if(is.null(fit)){fails[[mo]]<-data.frame(Fold=fold,Model=mo,Reason="fit/prediction failed");next}
    rv<-as.numeric(fit$rv); valid_ids<-colnames(vz)
    if(length(rv)!=length(valid_ids)){fails[[mo]]<-data.frame(Fold=fold,Model=mo,Reason=paste0("prediction length mismatch: ",length(rv)," vs ",length(valid_ids)));next}
    risks[[mo]]<-data.frame(Fold=fold,Model=mo,Sample=valid_ids,RiskScore=rv)
    mans[[mo]]<-data.frame(Fold=fold,Model=mo,Selector=sm,Learner=lr,SelectedGeneN=length(gs),FinalGeneN=length(unique(fit$genes)),SelectedContainsPRDX4="PRDX4"%in%gs,FinalContainsPRDX4="PRDX4"%in%fit$genes,SelectedGenes=paste(gs,collapse=";"),FinalGenes=paste(fit$genes,collapse=";"))
  }
  bind<-function(x)if(length(x))do.call(rbind,x)else data.frame();saveRDS(list(risk=bind(risks),manifest=bind(mans),fail=bind(fails)),ck)
  cat("Fold",fold,"saved:",ck,"\n")
}

fold_objs<-lapply(seq_len(k_folds),function(f)readRDS(file.path(out_dir,"Checkpoints",paste0("fold_",f,"_results.rds"))))
bind_nonempty<-function(x){x<-x[vapply(x,function(z)is.data.frame(z)&&nrow(z)>0,logical(1))];if(!length(x))return(data.frame());do.call(rbind,x)}
oof_risk<-bind_nonempty(lapply(fold_objs,`[[`,"risk"));oof_manifest<-bind_nonempty(lapply(fold_objs,`[[`,"manifest"));oof_fail<-bind_nonempty(lapply(fold_objs,`[[`,"fail"))
if(!nrow(oof_risk)||!nrow(oof_manifest))stop("No usable OOF predictions were produced")
write.csv(oof_risk,file.path(out_dir,"Tables","04_all_OOF_risk_scores.csv"),row.names=FALSE);write.csv(oof_manifest,file.path(out_dir,"Tables","05_fold_model_manifest.csv"),row.names=FALSE);write.csv(oof_fail,file.path(out_dir,"Tables","06_fold_failures.csv"),row.names=FALSE)

metric_list<-lapply(combos$Model,function(m){z<-oof_risk[oof_risk$Model==m,,drop=FALSE];z<-z[!duplicated(z$Sample),,drop=FALSE];ix<-match(z$Sample,emtab_samples);me<-calc_metrics(emtab_time[ix],emtab_event[ix],z$RiskScore);mn<-mean(oof_manifest$FinalGeneN[oof_manifest$Model==m],na.rm=TRUE);cbind(Model=m,CoverageN=nrow(z),CoverageFrac=nrow(z)/length(emtab_samples),MeanFinalGeneN=mn,me)})
ranking<-do.call(rbind,metric_list);numcols<-setdiff(colnames(ranking),"Model");ranking[numcols]<-lapply(ranking[numcols],safe_num);ranking$Eligible<-ranking$CoverageFrac>=.95 & is.finite(ranking$Cindex)
ranking$RankScore<-0;ranking$RankScore<-ranking$RankScore+ifelse(!is.na(ranking$CoxP)&ranking$CoxP<.05,2,0)+ifelse(!is.na(ranking$LogRankP)&ranking$LogRankP<.05,2,0)+ifelse(!is.na(ranking$HR)&ranking$HR>1,1,0)+ifelse(!is.na(ranking$KMHR)&ranking$KMHR>1,1,0)+ifelse(!is.na(ranking$Cindex)&ranking$Cindex>=.80,3,0)+ifelse(!is.na(ranking$Cindex)&ranking$Cindex>=.70,2,0)+ifelse(!is.na(ranking$Cindex)&ranking$Cindex>=.60,1,0)+ifelse(!is.na(ranking$MeanAUC)&ranking$MeanAUC>=.80,3,0)+ifelse(!is.na(ranking$MeanAUC)&ranking$MeanAUC>=.70,2,0)+ifelse(!is.na(ranking$MeanAUC)&ranking$MeanAUC>=.60,1,0)
ranking$ComplexityPenalty<-ifelse(ranking$MeanFinalGeneN>20,1,0);ranking$FinalRankScore<-ranking$RankScore-ranking$ComplexityPenalty
ranking<-ranking[order(!ranking$Eligible,-ranking$FinalRankScore,ranking$CoxP,ranking$LogRankP,-ranking$Cindex,-ranking$MeanAUC),,drop=FALSE];ranking$OOF_Rank<-seq_len(nrow(ranking));write.csv(ranking,file.path(out_dir,"Tables","07_EMTAB_OOF_101_model_ranking.csv"),row.names=FALSE);write.csv(head(ranking,20),file.path(out_dir,"Tables","08_EMTAB_OOF_top20.csv"),row.names=FALSE)
eligible<-ranking[ranking$Eligible,,drop=FALSE];if(!nrow(eligible))stop("No model has >=95% OOF coverage")
best_model<-eligible$Model[1];best_row<-combos[match(best_model,combos$Model),,drop=FALSE];best_selector<-best_row$Selector;best_learner<-best_row$Learner
m050_row<-ranking[ranking$Model=="M050_SuperPC_Ridge",,drop=FALSE]
cat("\nLOCKED PIPELINE FROM E-MTAB OOF ONLY:",best_model,"\n");cat("M050 OOF rank:",if(nrow(m050_row))m050_row$OOF_Rank else NA,"\n")
write.csv(data.frame(Item=c("Locked pipeline","Selector","Learner","M050 OOF rank","Locked pipeline is M050"),Value=c(best_model,best_selector,best_learner,if(nrow(m050_row))m050_row$OOF_Rank else NA,best_model=="M050_SuperPC_Ridge")),file.path(out_dir,"Tables","09_LOCKED_MODEL_DECISION_BEFORE_SEQC.csv"),row.names=FALSE)

## ---------------------------- lock/refit full E-MTAB ----------------------------
## Scaling is refit on all E-MTAB. SEQC survival remains untouched.
scfull<-scale_train_apply(emtab_expr[candidate_genes,,drop=FALSE],seqc_expr[candidate_genes,,drop=FALSE]);full_z<-scfull$train_z;seqc_z<-scfull$valid_z
set.seed(seed_main+900001);locked_selected<-select_genes(best_selector,candidate_genes,full_z,emtab_time,emtab_event,seed_main+900001)
locked_fit<-if(best_learner=="PLSRcox_extra")fit_plsr(locked_selected,full_z,seqc_z,emtab_time,emtab_event,best_model,seed_main+900002) else fit_learner(best_learner,locked_selected,full_z,seqc_z,emtab_time,emtab_event,rep(1,ncol(seqc_z)),rep(0,ncol(seqc_z)),best_model,seed_main+900002)
if(is.null(locked_fit))stop("Locked pipeline failed during full E-MTAB refit")
saveRDS(list(model=best_model,selector=best_selector,learner=best_learner,candidate_genes=candidate_genes,selected_genes=locked_selected,final_genes=locked_fit$genes,coefficients=locked_fit$coef,center=scfull$center,scale=scfull$scale),file.path(out_dir,"Tables","10_locked_model_object_light.rds"))
write.csv(locked_fit$coef,file.path(out_dir,"Tables","11_locked_model_coefficients_or_importance.csv"),row.names=FALSE)
write.csv(data.frame(Gene=candidate_genes,Center=scfull$center[candidate_genes],Scale=scfull$scale[candidate_genes]),file.path(out_dir,"Tables","12_full_EMTAB_scaling_parameters.csv"),row.names=FALSE)

## ---------------------------- NOW unblind SEQC outcomes ----------------------------
cat("\nSEQC model predictions are locked. Reading SEQC survival outcomes now.\n")
seqc_clin<-as.data.frame(seqc_obj$clin);if(!"Sample"%in%colnames(seqc_clin))seqc_clin$Sample<-rownames(seqc_clin)
te<-find_time_event_cols(seqc_clin);if(is.na(te$time)||is.na(te$event))stop("Cannot identify SEQC survival columns")
seqc_time<-safe_num(seqc_clin[[te$time]]);seqc_event<-clean_event(seqc_clin[[te$event]])
if(length(locked_fit$rv)!=ncol(seqc_z)||length(locked_fit$rt)!=ncol(full_z))stop("Locked-model prediction length mismatch")
sr<-data.frame(Sample=colnames(seqc_z),RiskScore=as.numeric(locked_fit$rv));ix<-match(sr$Sample,seqc_clin$Sample);sr$OS_time<-seqc_time[ix];sr$OS_event<-seqc_event[ix];sr<-sr[complete.cases(sr[,c("OS_time","OS_event","RiskScore")])&sr$OS_time>0,,drop=FALSE]
tr<-data.frame(Sample=colnames(full_z),OS_time=emtab_time[match(colnames(full_z),emtab_samples)],OS_event=emtab_event[match(colnames(full_z),emtab_samples)],RiskScore=as.numeric(locked_fit$rt))
write.csv(tr,file.path(out_dir,"Tables","13_locked_model_EMTAB_scores.csv"),row.names=FALSE);write.csv(sr,file.path(out_dir,"Tables","14_locked_model_SEQC_external_scores.csv"),row.names=FALSE)
perf<-rbind(cbind(Dataset="E-MTAB-179 full refit (apparent)",calc_metrics(tr$OS_time,tr$OS_event,tr$RiskScore)),cbind(Dataset="SEQC external validation",calc_metrics(sr$OS_time,sr$OS_event,sr$RiskScore)))
write.csv(perf,file.path(out_dir,"Tables","15_LOCKED_MODEL_FINAL_PERFORMANCE.csv"),row.names=FALSE)

final<-data.frame(Item=c("E-MTAB sample N","E-MTAB event N","SEQC external sample N","SEQC external event N","BALANCED_40 saved gene N","Cross-platform candidate gene N","OOF folds","OOF eligible model N","Locked model","Locked selector","Locked learner","M050 OOF rank","Locked model is historical M050","SEQC participated in OOF model selection"),Value=c(length(emtab_time),sum(emtab_event==1),nrow(sr),sum(sr$OS_event==1),length(candidate_query),length(candidate_genes),k_folds,sum(ranking$Eligible),best_model,best_selector,best_learner,if(nrow(m050_row))m050_row$OOF_Rank else NA,best_model=="M050_SuperPC_Ridge","FALSE"))
write.csv(final,file.path(out_dir,"Tables","16_FINAL_DECISION_SUMMARY.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out_dir,"Logs","sessionInfo.txt"))
cat("\n============================================================\n");cat("SC28 RECALC 04 completed.\n");print(final);cat("Upload the ENTIRE output folder as ZIP.\n");cat("End:",as.character(Sys.time()),"\n");cat("============================================================\n")
sink()
