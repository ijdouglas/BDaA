# MastersResearch

All scripts executed for this study are contained in this repo in the `scripts/FINAL` folder.

All scripts source the file `functions.R` (in `scripts/FINAL`)

**Pipeline Overview**

**I.** The indices of the train-test-tune splits of the data were created with:

- `split_data.R` (creating the file `splits.csv` contained in `data/FINAL`)

**II.** The final results were derived from running the following six scripts:

- `02_LogReg-CVPVI_AdjustNothing.R` (at 2021-01-12 15:36:21.670587024)
- `02_LogReg-CVPVI_AgeSexResidualize.R` (at 2021-01-12 15:37:40.330639448)
- `02_LogReg-CVPVI_AllResidualize.R` (at 2021-01-12 15:39:02.058693917)
- `02_sb-GBM-CVPVI_AgeSexResidualize.R` (at 2021-01-12 15:40:28.966751838)
- `02_sb-GBM-CVPVI_AllResidualize.R` (at 2021-01-12 15:41:54.542808871)
- `02_sb-GBM-CVPVI_noResidualize.R` (at 2021-01-12 15:43:39.314878697)

**III.** Bootstrap confidence intervals are then computed using:

- `AUC_and_K-S-D_BootCI.R` (at 2021-01-12 20:58:35.339318938)

**Results**

The above pipeline resulted in the following files

- `SB-LogReg-CVPVI_ADJUSTNOTHING_2021-01-12.rds` (created on Jan 12, 2021 at 15:36)
- `SB-LogReg-CVPVI_AgeSexAdjusted_2021-01-12.rds` (created on Jan 12, 2021 at 15:37)
- `SB-LogReg-CVPVI_AllResidualized_2021-01-12.rds` (created on Jan 12, 2021 at 15:39)
- `SB-GBM_CVPVI_AgeSexAdjusted_2021-01-12.rds` (created on Jan 12, 2021 at 15:40)
- `SB-GBM_CVPVI_AllResidualized_2021-01-12.rds` (created on Jan 12, 2021 at 15:41)
- `SB-GBM_CVPVI_ADJUSTNOTHING_2021-01-12.rds` (created on Jan 12 2021, at 15:43)
- `BOOTSTRAP_CVPVI-PIPELINE_2021-01-12.rds` (created on Jan 12, 2021 at 20:59:05.687339296)
- `BOOTSTRAP_CVPVI-PIPELINE_2021-01-12.rds` (created on Jan 12, 2021 at 18:50) 
