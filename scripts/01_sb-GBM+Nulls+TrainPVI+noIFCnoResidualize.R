# This version does not residualize, drops IFC subjects
#setup
print('ian j. douglas; contact: idouglas@utexas.edu')
print('This script runs the GBM pipeline without residualizing, and drops IFC subjects.')
t0 <- Sys.time()
todays_date <- as.character(Sys.Date())
output_path = paste0("../output/SB-GBM_final_ADJUSTNOTHING_noIFC_",
                     todays_date, 
                     ".rds")
print(paste0("Today's date: ", todays_date, ". Starting script today at ", t0, "."))
print(paste0('When complete, results will be written to ', output_path))
nSplits <- 100
nNulls <- 100
set.seed(100)
# Part A: load requirements and data ----
## Loading packages
suppressMessages(library(MLmetrics))
suppressMessages(library(tidyverse))
suppressMessages(library(mgcv))
suppressMessages(library(gbm))
suppressMessages(library(parallel))
suppressMessages(library(vip))
#suppressMessages(library(fastshap))
source('functions.R')

## Read in the data
SB = readRDS("../data/SB.rds")
ELFK = readRDS('../data/ELFK.rds')
# Read in the list of IFC subjects
ifc_list = read.csv('../data/ifc_list.csv')$id
# Filter them out for ELFK
ELFK = ELFK %>%
  filter(!IDENT_SUBID %in% ifc_list)
## Read in the 100 partitions, generated randomly via 'split_data.R'
splits <- suppressWarnings(read_csv("../data/splits.csv"))

## Preprocessing: NONE FOR THIS PIPELINE

# Part B: Run the pipeline ----
# # 1. Split the data into train and test sets nSplits times (return all in a list)
ttSplits = purrr::map(unique(splits$iter), ~{
  frac <- splits %>% filter(iter == .x) %>% slice(1) %>% .$frac
  .IDX <- splits %>% filter(iter == .x) %>% .$IDX
  out <- list(
    train = SB[.IDX, ], # Using SB data from global env
    tunedata = SB[.IDX[-1:-(length(.IDX)*frac)], ], # Using SB data from global env
    trainfrac = frac,
    train.indices = .IDX,
    test = ELFK # Using ELFK data from global env
  )
  out # returned
})
# ttSplits = tuneSplits(nSplits, SB) %>%
#   # and for each one append ELFK as test data
#   lapply(., function(x) {.x <- x; .x$test <- ELFK; .x}) 
# print('Random splits generated')
# print(paste0('Initiating classification pipeline at ', Sys.time(), '.'))
# print(paste0('When complete, results will write out to ', output_path))

# 2. Fit the GBM to each, get CVPVI, get 100 nulls, and 100 null CVPVI
pipeline_results = mclapply(X=ttSplits, mc.cores = 13, FUN=function(Z) {
  # 1 of 2: Observed model
  # Fit the model with gbm::gbm() directly
  gbm_mod <- gbm(
    GROUP ~ .-subjWeight-IDENT_SUBID, 
    data = get_design(Z$train), # supply only the train data
    weights = subjWeight,
    distribution = "bernoulli", # model the binary outcome with a logit link fn
    n.trees = 750, # a large number to ensure the minimum criterion is reached
    #shrinkage = .1, # defaults
    interaction.depth = 6, # allow for interactions
    #n.minobsinnode = 1, # defaults
    train.fraction = Z$trainfrac, # construct the gbm on the first train.fraction*nrow(data) rows
    bag.fraction = 1 # do not resample from the train set, we manually created a test set.
  )
  # Pre-defined function to extract data frame of test set stats
  results = get_gbm_results(model = gbm_mod, Test = Z$test)
  # Pre-defined function to extract CVPVI for each variable
  cvpvi = get_gbm_cvpvi(gbm_mod, Test = Z$test, null =FALSE)
  # Additionally, compute the permutation importance on the train data
  train_pvi = get_gbm_cvpvi(gbm_mod, Test = Z$train[1:(nrow(Z$train)*Z$trainfrac), ], null = F)
  # # Additionally, get the SHAP values as they pertained to the train data
  # shap <- suppressMessages(
  #   vi_shap(object=gbm_mod, 
  #           feature_names = attr(gbm_mod$Terms, 'term.labels'),
  #           train = get_design(Z$train),
  #           nsim = 100,
  #           pred_wrapper = predict.gbm_bestTree)
  # )
  
  # 2 of 2: Null models.
  # Now get the 100 null results for each of these
  # null.results = purrr::map(1:nNulls, ~{
  #   set.seed(.x)
  #   # Make the null train (and tune) data:
  #   null_train_ind <- sample(1:(length(Z$train.indices)*Z$trainfrac))
  #   null_tune_ind <- sample(((length(Z$train.indices)*Z$trainfrac)+1):length(Z$train.indices))
  #   null_train <- get_design(Z$train)
  #   null_train$GROUP <- null_train$GROUP[c(null_train_ind, null_tune_ind)]
  #   # fit a null model
  #   null_gbm <- gbm(
  #     GROUP ~ .-subjWeight-IDENT_SUBID, 
  #     data = null_train, # supply only the train data
  #     weights = subjWeight,
  #     distribution = "bernoulli", # model the binary outcome with a logit link fn
  #     n.trees = 750, # a large number to ensure the minimum criterion is reached
  #     #shrinkage = .1, # defaults
  #     interaction.depth = 6, # allow for interactions
  #     #n.minobsinnode = 1, # defaults
  #     train.fraction = Z$trainfrac, # construct the gbm on the first train.fraction*nrow(data) rows
  #     bag.fraction = 1 # do not resample from the train set, we manually created a test set.
  #   )
  #   # nullify the test response vector
  #   null_idx_test <- sample(1:nrow(Z$test))
  #   # Shuffle the y variable, along with the confounds in the SAME order
  #   # (just break correlation with the X variables)
  #   nullified = Z$test %>% 
  #     mutate_at(.vars = vars(all_of(names(get_labels(Z$test, 'GROUP')))), 
  #               .funs = ~.[null_idx_test]) # put them all in the same shuffled order
  #   # Score the null model, to the null test data
  #   null.model_results = get_gbm_results(model = null_gbm, Test = nullified)
  #   null.cvpvi = get_gbm_cvpvi(null_gbm, Test = nullified, null = TRUE)
  #   null.train_pvi = get_gbm_cvpvi(null_gbm, Test = null_train[1:(nrow(Z$train)*Z$trainfrac), ], null = T)
  #   # null.shapvals <- suppressMessages(
  #   #   vi_shap(object=null_gbm, 
  #   #           feature_names = attr(null_gbm$Terms, 'term.labels'),
  #   #           train = null_train,
  #   #           nsim = 30,
  #   #           pred_wrapper = predict.gbm_bestTree)
  #   # )
  #   # Return:
  #   list('null.model.res' = null.model_results, 
  #        'null.cvpvi' = null.cvpvi,
  #        'null.train_pvi' = null.train_pvi,
  #        'null.shuffle' = null_idx_test) 
  # })
  
  # Package results and return
  # Include the number of trees used in the best model
  nt <- gbm.perf(gbm_mod, method = 'test', plot.it = F)
  out <- list(
    'Observed' = list('Model' = results, 
		      'CVPVI' = cvpvi, 
		      'train_PVI' = train_pvi, 
		      'num_trees' = nt),
    # 'Null' = list('Model' = lapply(null.results, function(x) x$null.model.res),
    #               'CVPVI' = lapply(null.results, function(x) x$null.cvpvi),
    #               'tPVI' = lapply(null.results, function(x) x$null.train_pvi),
    #               'null.shuffle' = lapply(null.results, function(x) x$null.shuffle)),
    # return some information about the test and tune data
    'data' = list(train.indices = Z$train.indices, trainfrac = Z$trainfrac)
  )
  return(out)
})


# Part C: Save the result ----
saveRDS(pipeline_results, output_path)
print('Done')
Sys.time() - t0
