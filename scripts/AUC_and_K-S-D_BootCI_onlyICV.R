# Bootstrap Confidence Intervals For ROCAUC and K-S Statistic
# The bootstrap will resample, with replacement, models, so that the predicted
# probabilities for ALL subjects, pertaining to a given model, will be resampled
# as a whole dataset, and added to a bootstrapped dataset with the predictions
# from other models (also for all subjects), as well as, due to sampling with
# replacement, the same chunk of predictions (for all subjects) for models
# that are resampled more than once. When resampling is finished, average
# predicted probabilities for each subject are computed, resulting in a single
# vector of predicted probabilities, containing one (bootstrapped-average)
# prediction per subject. This vector is compared to the vector of true group
# labels (also, by definition, containing one value per subject) in order to 
# compute ROCAUC. In addition, the vector is used in combination with that 
# true-group label vector to compute the Kolmogorov-Smirnov test statistic.

# 0. Setup
library(tidyverse)
library(MLmetrics)
library(parallel)
D_stat = function(y_pred, y_true) {
  all_pred_prob = y_pred
  pi_idx = y_true == 1
  comp_idx = y_true == 0
  ks_test <- suppressWarnings(ks.test(x = all_pred_prob[pi_idx],
                                      y = all_pred_prob[comp_idx],
                                      'less'))
  ks_test$stat
}

# 1. Set the working directory to be the same as the location of 02_ProcessResults.Rmd
# setwd(getwd())

# 2. Read in the GBM results and logistic regression results; setup some global variables
res_list = list(
  'GBM_onlyICV.Residualize' = readRDS('../output/SB-GBM+onlyICV-Residualize+noIFC_2021-01-10.rds')
)
L.res_list <- list(
  'GLM_onlyICV.Residualize' = readRDS('../output/SB-LogReg+onlyICV-Residualize+noIFC_2021-01-10.rds')
)
keys <- names(res_list) %>% sub(pattern = "^G.M_", replacement = "", x = .)

# 3. Get the ROCAUC for each
average_performance_boot_final <- map_dfr(keys, .f = ~{
  gbm <- map_dfr(pluck(res_list, paste0('GBM_', .x)), # Returns a list, one for each resample
                 function(x) x$Observed$Model$Predictions) %>%
    group_by(IDENT_SUBID) %>% mutate(preproc = .x, iteration = 1:n(), model = 'GBM')
  glm <- map_dfr(pluck(L.res_list, paste0('GLM_', .x)),
                 function(x) x$Observed$Model$Predictions) %>%
    group_by(IDENT_SUBID) %>% mutate(preproc = .x, iteration = 1:n(), model = 'GLM')
  rbind(gbm, glm) %>%
    remove_rownames()
}) %>% mutate(x_axis = paste(model, preproc, sep = "_")) %>%
  {
    assign('average_performance_boot', ., pos = .GlobalEnv)
    .
  } %>%
  .$x_axis %>% unique %>%
  map_dfr(., function(Z) {
    tmp <- average_performance_boot %>% filter(x_axis == Z)
    #assign("tmp", tmp, pos = .GlobalEnv)
    out <- mclapply(X = 1:1000, mc.cores = 12, FUN = function(seed) {
      set.seed(seed)
      boot_idx <- sample(100, replace = T)
      boot_dat <- map_dfr(boot_idx, function(i) {tmp %>% filter(iteration==i)})
      dim(boot_dat) ->> dim.boot_dat
      boot_perf <- boot_dat %>%
        dplyr::select(IDENT_SUBID, GROUP, preproc, model, x_axis, pred_prob) %>%
        group_by(IDENT_SUBID, GROUP, preproc, model, x_axis) %>%
        # Collapse over iterations, one result per subject now
        summarise(pred_prob = mean(pred_prob)) %>%
        # Now ungroup and summarize into one scalar ROCAUC
        ungroup() %>%
        summarise(rocauc = MLmetrics::AUC(pred_prob, GROUP)) %>%
        mutate(bootstrap = seed, x_axis = Z) %>%
        dplyr::select(bootstrap, rocauc, x_axis)
      boot_perf
    })
    return(reduce(out, rbind))
  })

# 3. Get the K-S test statistic
average_D_boot_final <- map_dfr(keys, .f = ~{
  gbm <- map_dfr(pluck(res_list, paste0('GBM_', .x)), # Returns a list, one for each resample
                 function(x) x$Observed$Model$Predictions) %>%
    group_by(IDENT_SUBID) %>% mutate(preproc = .x, iteration = 1:n(), model = 'GBM')
  glm <- map_dfr(pluck(L.res_list, paste0('GLM_', .x)),
                 function(x) x$Observed$Model$Predictions) %>%
    group_by(IDENT_SUBID) %>% mutate(preproc = .x, iteration = 1:n(), model = 'GLM')
  rbind(gbm, glm) %>%
    remove_rownames()
}) %>% mutate(x_axis = paste(model, preproc, sep = "_")) %>%
  {
    assign('average_D_boot', ., pos = .GlobalEnv)
    .
  } %>%
  .$x_axis %>% unique %>%
  map_dfr(., function(Z) {
    tmp <- average_performance_boot %>% filter(x_axis == Z)
    #assign("tmp", tmp, pos = .GlobalEnv)
    out <- mclapply(X = 1:1000, mc.cores = 12, FUN = function(seed) {
      set.seed(seed)
      boot_idx <- sample(100, replace = T)
      boot_dat <- map_dfr(boot_idx, function(i) {tmp%>%filter(iteration==i)})
      dim(boot_dat) ->> dim.boot_dat
      boot_perf <- boot_dat %>%
        dplyr::select(IDENT_SUBID, GROUP, preproc, model, x_axis, pred_prob) %>%
        group_by(IDENT_SUBID, GROUP, preproc, model, x_axis) %>%
        # Collapse over iterations, one result per subject now
        summarise(pred_prob = mean(pred_prob)) %>%
        # Now ungroup and summarize into one scalar ROCAUC
        ungroup() %>%
        summarise(D_stat = D_stat(pred_prob, GROUP)) %>%
        mutate(bootstrap = seed, x_axis = Z) %>%
        dplyr::select(bootstrap, D_stat, x_axis)
      boot_perf
    })
    return(reduce(out, rbind))
  })

OUTPUT <- left_join(average_performance_boot_final, 
                    average_D_boot_final, 
                    by = c('bootstrap', 'x_axis'))

saveRDS(OUTPUT, paste0("../output/BOOTSTRAP_ROCAUC_AND_K-S-D_", Sys.Date(),".rds"))
