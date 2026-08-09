# PRDX4–M070 neuroblastoma reproducibility package

This archive accompanies the manuscript **“Single-cell lactate state-guided prioritization of PRDX4 and a locked multigene risk model jointly characterize metabolic–redox adaptation and immune hyporesponsiveness in neuroblastoma.”**

## What this package contains

- The recovered source scripts that were used for the audited recalculation sequence from MSigDB lactate-state construction through the final common-gene-universe sensitivity analysis.
- The serialized, locked M070 model used for the final external-validation workflow.
- Original stage result archives or selected final output tables, figures, logs, and R session records.
- File manifests and cryptographic checksums.

The source files are preserved under their executed names. They were recovered from the author’s saved analysis artifacts; they have not been silently rewritten as if they were the executed originals.

## Final analysis sequence

1. `PRDX4_RECALC_01_MSigDB_LactateState.R` — retrieves the five MSigDB lactate-related sets, forms the 284-gene union, calculates LactateScore, and assigns exact within-sample LactateHigh/Mid/Low states.
2. `PRDX4_RECALC_02_ProliferationBalanced_Pseudobulk.R` — proliferation balancing, paired pseudobulk analysis, PRDX4 rank and robustness summaries.
3. `PRDX4_RECALC_03_CandidatePool_M050_ImpactAudit.R` — historical candidate-pool audit only. The filename contains M050 because this stage examined the earlier workflow; it did not select the final model.
4. `PRDX4_RECALC_04_Clean101ML_EMTAB_CV_SEQC_External.R` — transitional clean 101-pipeline rerun and platform-coverage audit. This was not the final model-selection evidence.
5. `PRDX4_RECALC_05_RepeatedCV_SelectionAudit.R` — final 10-repeat × 5-fold E-MTAB-179-only model comparison and stability audit. SEQC survival outcomes were not accessed.
6. `PRDX4_RECALC_06_LockM070_SEQC_External.R` — refits the chosen M070 SuperPC–RSF pipeline on all E-MTAB-179 data, serializes the model, records the lock, and documents the initial cross-platform prediction failure.
7. `PRDX4_RECALC_06B_M070_SEQC_CohortZ_External.R` — outcome-blind expression-scale correction using within-SEQC gene-wise z standardization, followed by external validation with the unchanged locked model.
8. `PRDX4_RECALC_07A_M070_Prognostic_ClinicalIncremental.R` — prognostic and clinical incremental-value refresh.
9. `PRDX4_RECALC_07A_FIX1_ClinicalDCA_AgeUnit.R` — final correction for the E-MTAB-179 age unit and decision-curve implementation. Use FIX1 rather than the uncorrected 07A DCA/calibration outputs.
10. `PRDX4_RECALC_07B_FIX3_M070_Pathway_ImmuneMechanism.R` — final pathway, immune, PRDX4, and selected mechanistic association analysis.
11. `PRDX4_RECALC_07C_CommonGeneUniverse_PathwaySensitivity.R` — final common-gene-universe cross-platform sensitivity analysis. This is the terminal analysis stage.

## Frozen model

- Model: M070, SuperPC selector followed by random survival forest.
- Selected features: 20 genes.
- Locked model file: `Models/M070_locked_model_BEFORE_SEQC.rds`.
- MD5: `e719d0e6b626e75cc134e83fd0929083`.
- SHA-256: `2b8aac34f99f11049252a785c2b783628defde476efd74abd117692acaff6909`.

The model must not be refitted, its features must not be reselected, and the cutoff must not be optimized in SEQC.

## Re-running the workflow

The scripts were executed on Windows with the project root:

`C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026`

and the primary server-input folder:

`C:/Users/Administrator/Desktop/NB_Lactate_Bulk_ServerInput_20260522_112238`

If the files are stored elsewhere, update only the path definitions near the top of each script. Run each entire script from line 1 in the order above; do not execute isolated blocks from stateful scripts. Several scripts install or require CRAN/Bioconductor packages. Exact packages and versions from completed stages are recorded in `Session_Info/` and in the original stage logs.

Large primary expression objects are not duplicated in this repository-ready archive. They are public-data derivatives and must be reconstructed or placed at the paths listed in `Documentation/Input_Data_Manifest.csv`.

## Interpretation safeguards

- The final model-selection evidence is RECALC 05, not the transitional RECALC 04 run and not the historical M050 analysis.
- The locked M070 object was serialized before SEQC survival outcomes were accessed.
- RECALC 06B is a post-audit, outcome-blind correction for expression-scale incompatibility; it must not be described as a prespecified calibration rule.
- The corrected age-unit and DCA outputs are from RECALC 07A FIX1.
- The final pathway/immune outputs are from RECALC 07B FIX3 and the final cross-platform sensitivity outputs are from RECALC 07C.
- CellChat and virtual-knockout findings are exploratory and are not causal evidence.

## Reproducibility status

The archive supports source-code inspection, model-integrity verification, result tracing, and rerunning in the author’s original file layout. It is not a one-click containerized workflow because the large public-data input objects and all R package binaries are not embedded. The repository-release candidate includes an MIT License and citation metadata; the authors must approve the license and replace the GitHub/Zenodo placeholders before public release.
