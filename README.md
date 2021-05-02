# MastersResearch

All scripts executed for this study are contained in this repo in the `scripts/FINAL` folder.

All scripts source the file `functions.R` (in `scripts/FINAL`)

**Pipeline Overview**

**I.** The indices of the train-test-tune splits of the data were created with:

```split_data.R #creating the file data/splits.csv contained in data)```

* In reality, `split_data.R` is run twice:
    * Once to split the SB younger subsample (output file: `data/younger.splits.csv`)
    * Once to split the SB older subsample (output file: `data/older.splits.csv`)
* Older and younger samples for both SB and ELFK are defined as follows:
    1. Compute the `floor()` of every subject's age.
    2. Compute the median (floor of) the age vector in **ELFK**
    3. Younger {SB, ELFK} samples are those whose (floored) age is `<=` median (floored) **ELFK** age.
    4. Older {SB, ELFK} samples are those with a floored age `>` median (floored) **ELFK** age.

**II.** After splitting the data by age, the following scripts conduct analyses:

- `02-sb-GBM-CVPVI_unadjusted_older.R`: this script reads in the older {SB, ELFK} data for {training, testing} gradient boosting models
- `02-sb-GBM-CVPVI_unadjusted_younger.R`: this script reads in the younger {SB, ELFK} data for {training, testing} gradient boosting models
- `02-LogReg-CVPVI_unadjusted_older.R`: this script reads in the older {SB, ELFK} data for {training, testing} logistic regressions
- `02-LogReg-CVPVI_unadjusted_younger.R`: this script reads in the younger {SB, ELFK} data for {training, testing} logistic regressions

**NOTE.** The OLD version used the following script to run the pipeline (having not counterbalanced by age):

- `02_LogReg-CVPVI_AdjustNothing.R` (at 2021-01-12 15:36:21.670587024)
- `02_LogReg-CVPVI_AgeSexResidualize.R` (at 2021-01-12 15:37:40.330639448)
- `02_LogReg-CVPVI_AllResidualize.R` (at 2021-01-12 15:39:02.058693917)
- `02_sb-GBM-CVPVI_AgeSexResidualize.R` (at 2021-01-12 15:40:28.966751838)
- `02_sb-GBM-CVPVI_AllResidualize.R` (at 2021-01-12 15:41:54.542808871)
- `02_sb-GBM-CVPVI_noResidualize.R` (at 2021-01-12 15:43:39.314878697)


**III.** Bootstrap confidence intervals are then computed using:

- `AUC_and_K-S-D_BootCI.R` (at 2021-01-12 20:58:35.339318938)

**NOTE.** The old version (with the 6 old not-age-specific pipeline scripts) ran the following (at the following time)

- `AUC_and_K-S-D_BootCI.R` (at 2021-01-12 20:58:35.339318938)

**IV. Results.** The outputs of the 4 pipeline scripts (one for each age group and model type) are as follows:

- `SB-GBM-CVPVI_unadjusted-older_2021-05-02.rds`
- `SB-GBM-CVPVI_unadjusted-younger_2021-05-02.rds` 
- `SB-LogReg-CVPVI_unadjusted-older_2021-05-02.rds`
- `SB-LogReg-CVPVI_unadjusted-younger_2021-05-02.rds`


**NOTE.** The most recent old results are the following:

- `SB-LogReg-CVPVI_ADJUSTNOTHING_2021-01-12.rds` (created on Jan 12, 2021 at 15:36)
- `SB-LogReg-CVPVI_AgeSexAdjusted_2021-01-12.rds` (created on Jan 12, 2021 at 15:37)
- `SB-LogReg-CVPVI_AllResidualized_2021-01-12.rds` (created on Jan 12, 2021 at 15:39)
- `SB-GBM_CVPVI_AgeSexAdjusted_2021-01-12.rds` (created on Jan 12, 2021 at 15:40)
- `SB-GBM_CVPVI_AllResidualized_2021-01-12.rds` (created on Jan 12, 2021 at 15:41)
- `SB-GBM_CVPVI_ADJUSTNOTHING_2021-01-12.rds` (created on Jan 12 2021, at 15:43)
- `BOOTSTRAP_CVPVI-PIPELINE_2021-01-12.rds` (created on Jan 12, 2021 at 20:59:05.687339296)

**V. Processing Results.** 

The script `process_results.R` reads in (hard-coded) results, such as those in `IV` scores the models, conducts permutation tests of model performance, and calculates CVPVI for the ensemble.
