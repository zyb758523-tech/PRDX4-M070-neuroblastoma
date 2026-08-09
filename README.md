[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21863956.svg)](https://doi.org/10.5281/zenodo.21863956)

# PRDX4–M070 neuroblastoma analysis

This repository contains the audited R code, frozen M070 model, selected outputs, stage archives, session records, and integrity checks supporting the manuscript:

> **Single-cell lactate state-guided prioritization of PRDX4 and a locked multigene risk model jointly characterize metabolic–redox adaptation and immune hyporesponsiveness in neuroblastoma**

## Repository status

- Archived public release: v1.0.0
- Final prognostic model: M070 (`SuperPC` selector followed by random survival forest)
- Model selection: 10 repeats × 5-fold cross-validation in E-MTAB-179 only
- External validation: SEQC, 498 patients and 105 overall-survival events
- Frozen model MD5: `e719d0e6b626e75cc134e83fd0929083`
- Frozen model SHA-256: `2b8aac34f99f11049252a785c2b783628defde476efd74abd117692acaff6909`

SEQC survival outcomes were not used for model selection, feature selection, model fitting, expression harmonization, or cutoff optimization. The final SEQC analysis used outcome-blind, within-cohort gene-wise z standardization to address the audited raw-scale incompatibility between platforms, followed by evaluation with the unchanged frozen M070 model and the E-MTAB-179-derived cutoff.

## Start here

- [`00_README_FIRST.md`](00_README_FIRST.md): complete workflow, rerun instructions, and interpretation safeguards
- [`Documentation/Stage_Manifest.csv`](Documentation/Stage_Manifest.csv): stage-by-stage provenance
- [`Documentation/Input_Data_Manifest.csv`](Documentation/Input_Data_Manifest.csv): required public data and local derivative inputs
- [`Documentation/SHA256SUMS.txt`](Documentation/SHA256SUMS.txt): integrity checks for repository files
- [`CODE_AVAILABILITY.md`](CODE_AVAILABILITY.md): manuscript-ready software availability fields
- [`PUBLIC_RELEASE_GUIDE_CN.md`](PUBLIC_RELEASE_GUIDE_CN.md): Chinese GitHub and Zenodo release instructions

## Final analysis sequence

Run the scripts in `Code/` in the following order:

1. `PRDX4_RECALC_01_MSigDB_LactateState.R`
2. `PRDX4_RECALC_02_ProliferationBalanced_Pseudobulk.R`
3. `PRDX4_RECALC_03_CandidatePool_M050_ImpactAudit.R` — historical audit only
4. `PRDX4_RECALC_04_Clean101ML_EMTAB_CV_SEQC_External.R` — transitional rerun/platform audit
5. `PRDX4_RECALC_05_RepeatedCV_SelectionAudit.R` — final model-selection evidence
6. `PRDX4_RECALC_06_LockM070_SEQC_External.R`
7. `PRDX4_RECALC_06B_M070_SEQC_CohortZ_External.R`
8. `PRDX4_RECALC_07A_M070_Prognostic_ClinicalIncremental.R`
9. `PRDX4_RECALC_07A_FIX1_ClinicalDCA_AgeUnit.R` — use FIX1 for final clinical/DCA outputs
10. `PRDX4_RECALC_07B_FIX3_M070_Pathway_ImmuneMechanism.R`
11. `PRDX4_RECALC_07C_CommonGeneUniverse_PathwaySensitivity.R`

## Data

The study uses public datasets, including GSE137804, E-MTAB-179, and SEQC/GSE49711, together with the MSigDB gene sets documented in the manuscript. Large primary expression objects are not duplicated in this repository. Their accessions, roles, and expected local derivative files are listed in `Documentation/Input_Data_Manifest.csv`.

No directly identifying clinical information is included. Patient/sample-level files in the repository use public cohort identifiers and derived analysis values.

## Environment

The workflow was executed in R on Windows. Completed-stage package versions and platform details are preserved in `Session_Info/` and in each archived result log. The scripts use explicit project-root paths; users running the workflow elsewhere should update only the path definitions near the beginning of each script.

## License

The analysis code is prepared for release under the MIT License. Public release should occur only after the authors have confirmed that this license is consistent with institutional requirements and the terms of the upstream data resources.

## Citation

Citation metadata are provided in `CITATION.cff`. After the first Zenodo archive is published, replace the repository and DOI placeholders in `CITATION.cff`, `CODE_AVAILABILITY.md`, and the manuscript with the final GitHub URL and version-specific Zenodo DOI.

