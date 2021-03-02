# Script drops IFC, residualizes brain features for age, sex, and ICV. Computes
# PVI with train data, and a null distribution for that, plus the parametere
# estimates, plus the CVPVI.

#setup
print('ian j. douglas; contact: idouglas@utexas.edu')
t0 <- Sys.time()
todays_date <- as.character(Sys.Date())
output_path = paste0("../output/SB-LogReg+NullModels_noICVResidualize+noIFC",
                     todays_date, 
                     ".rds")
print(paste0("Today's date: ", todays_date, ". Starting script today at ", t0, "."))
nSplits <- 100
nNulls <- 100

# Part A: load requirements and data ----
## Loading packages
suppressMessages(library(MLmetrics))
suppressMessages(library(tidyverse))
suppressMessages(library(mgcv))
#suppressMessages(library(gbm))
suppressMessages(library(parallel))
suppressMessages(library(vip))
source('functions.R')

## Reading in the data
SB = readRDS("../data/SB.rds") %>% mutate_at('IDENT_SUBID', as.character)
ELFK = readRDS('../data/ELFK.rds') %>% mutate_at('IDENT_SUBID', as.character)
# Drop the IFC subjects
ifc_list = read.csv('../data/ifc_list.csv')$id
# Filter them out for ELFK
ELFK = ELFK %>%
  filter(!IDENT_SUBID %in% ifc_list)
## Read in the 100 partitions, generated randomly via 'split_data.R'
splits <- suppressWarnings(read_csv("../data/splits.csv"))

# Conduct residualization for ALL 3 covariates
SB = cbind(SB %>% select(-starts_with(match = c("L.", 'R.'))),
           imap(SB %>% select(starts_with(match = c("L.", 'R.'))), ~Residualize_noICV(.y, SB)))
ELFK = cbind(ELFK %>% select(-starts_with(match = c("L.", 'R.'))),
             imap(ELFK %>% select(starts_with(match = c("L.", 'R.'))), ~Residualize_noICV(.y, ELFK)))

# Part B: Run the pipeline ----
# 1. Get the list of splits
ttSplits = purrr::map(unique(splits$iter), ~{
  frac <- splits %>% filter(iter == .x) %>% slice(1) %>% .$frac
  .IDX <- splits %>% filter(iter == .x) %>% .$IDX
  out <- list(
    train = SB[.IDX[1:(length(.IDX)*frac)], ], # Using SB data from global env
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
print(paste0('Initiate modeling. When complete, results will save to: ', output_path))
# 2. Fit the logistic reg to each, get CVPVI, get 100 nulls, and 100 null CVPVI
pipeline_results = mclapply(X=ttSplits, mc.cores = 9, FUN=function(Z) {
  # Fit the model with glm() directly
  glm_mod <- glm(
    GROUP ~ .-subjWeight, 
    data = get_design(Z$train) %>% select(-one_of('IDENT_SUBID')), # supply only the train data
    weights = subjWeight,
    family = "quasibinomial"
  )
  # Pre-defined function to extract data frame of test set stats
  results = get_glm_results(model = glm_mod, Test = Z$test %>% mutate_at('IDENT_SUBID', as.character))
  # Pre-defined function to extract CVPVI for each variable
  cvpvi = get_glm_cvpvi(glm_mod, Test = Z$test %>% mutate_at('IDENT_SUBID', as.character), null =FALSE)
  # Additionally, compute the permutation importance on the train data
  train_pvi = get_glm_cvpvi(glm_mod, Test = Z$train %>% mutate_at('IDENT_SUBID', as.character), null = F)
  # Collect the parameter estimates as well
  train_coef = coef(glm_mod)
  # Now get the 100 null results for each of these
  # null.results = purrr::map(1:nNulls, ~{
  #   set.seed(.x)
  #   null_train_ind <- sample(nrow(Z$train))
  #   null_train <- get_design(Z$train)
  #   # make the group variable random
  #   null_train$GROUP <- null_train$GROUP[null_train_ind]
  #   null_train$IDENT_SUBID <- as.character(null_train$IDENT_SUBID)
  #   # fit a null model
  #   null_glm <- glm(
  #     GROUP ~ .-subjWeight, 
  #     data = null_train %>% select(-one_of('IDENT_SUBID')), # supply only the train data
  #     weights = subjWeight,
  #     family = "quasibinomial" # model the binary outcome with a logit link fn
  #   )
  #   # nullify the test response vector
  #   null_idx_test <- sample(nrow(Z$test))
  #   # Shuffle the y variable, along with the confounds in the SAME order
  #   # (just break correlation with the X variables)
  #   nullified = Z$test %>% 
  #     mutate_at(.vars = vars(all_of(names(get_labels(Z$test, 'GROUP')))), 
  #               .funs = ~.[null_idx_test]) %>% # putting them all in the same shuffled order
  #     mutate_at('IDENT_SUBID', as.character)
  #   # Score the null model, to the null test data
  #   null.model_results = get_glm_results(model = null_glm, Test = nullified)
  #   null.cvpvi = get_glm_cvpvi(null_glm, Test = nullified %>% mutate_at('IDENT_SUBID', as.character), null = TRUE)
  #   null.train_pvi = get_glm_cvpvi(null_glm, Test = null_train %>% mutate_at('IDENT_SUBID', as.character), null = T)
  #   list('null.model.res' = null.model_results, 
  #        'null.cvpvi' = null.cvpvi,
  #        'null.train_pvi'= null.train_pvi,
  #        'null.coef' = coef(null_glm), # return the model coefficients
  #        'null.shuffle' = null_idx_test) # returned
  # })
  # Package and return
  out <- list(
    'Observed' = list('Model' = results, 
                      'CVPVI' = cvpvi,
                      'train_PVI' = train_pvi, 
                      'train_coef' = train_coef),
    # 'Null' = list('Model' = lapply(null.results, function(x) x$null.model.res),
    #               'CVPVI' = reduce(lapply(null.results, function(x) x$null.cvpvi), rbind),
    #               'tPVI' = reduce(lapply(null.results, function(x) x$null.train_pvi), rbind),
    #               'tCoefs' = reduce(lapply(null.results, function(x) x$null.coef), rbind)[,-1], # drop the intercept
    #               'null.shuffle' = lapply(null.results, function(x) x$null.shuffle)),
    # return some information about the test and tune data
    'data' = list(train.indices = Z$train.indices)
  )
  return(out)
})


# Part C: Save the result ----
saveRDS(pipeline_results, output_path)
print('Done')
print('Elapsed:')
Sys.time() - t0
