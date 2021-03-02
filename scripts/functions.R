# Script containing subfunctions used by 01_GBM-Classify-Pipeline.R

# 1. A wrapper rbind to bypass incompatible rownames
Rbind = function(x, y) {
  rownames(x) <- rownames(y) <- NULL
  rbind(x, y)
}

# 2. A data cleaning function to be used inside of the function that splits the data
fixup = function(x) { 
  x %>% 
    group_by(IDENT_SUBID) %>% mutate(subjWeight = 1 / n()) %>% ungroup() %>%
    # Reorder columns. Confounds and labels are put first, then SGMV features
    dplyr::select(IDENT_SUBID, GROUP, brain_age_yrs, ICV, GENDER_FEMALE, subjWeight, matches("^L\\.|^R\\."))
}

# 3. Make splits of the data into train, tune, and test
makeSplits = function(n, .data, preproc = FALSE)
{ 
  data <- .data
  dataPartitions = vector(mode = 'list', length = 0)
  data <- data %>% group_by(IDENT_SUBID) %>% mutate(subjWeight = 1 / n()) %>% ungroup() 
  pi.index = which(as.logical(data$GROUP)) # converts 1 to TRUE, passes along to which()
  comp.index = which(!as.logical(data$GROUP)) # converts 1 to F, reverses, sends to which()
  pi.subid = data$IDENT_SUBID[pi.index] # subject ids for the pi group
  comp.subid = data$IDENT_SUBID[comp.index] # same for comps
  pi.weights = (1/n_distinct(pi.subid)) * data$subjWeight[pi.index] # define sampling weights
  comp.weights = (1/n_distinct(comp.subid)) * data$subjWeight[comp.index] # define sampling weights
  # define the train set to be balanced with respect to the groups
  for (i in 1:n) {
    # sample into a training set:
    train.index = c( 
      # create a sample with 70% of the number of people (subject ids)
      sample(pi.index, size = round(n_distinct(pi.subid)*.75), prob = pi.weights), # 70 / 30 split
      sample(comp.index, size = round(n_distinct(pi.subid)*.75), prob = comp.weights) # downsample
    )
    traindata = data[sample(train.index), ] # shuffle the order
    # test data is now those SUBJECTS (not rows) 
    testdata = data[!data$IDENT_SUBID %in% traindata$IDENT_SUBID, ] 
    # Now make a modification to the train data so that the first 70% of
    # the rows contain the true train subset that will be heldout for an internal tuning loop
    # by the gbm() function.
    ids <- unique(sample(traindata$IDENT_SUBID)) # shuffle the order again
    indexer <- 1
    tunedata <- traindata[0, ] # Create an empty data frame with the same column names
    while(nrow(tunedata)/nrow(traindata) < .33) {
      tunedata <- Rbind(tunedata, 
                        traindata %>% dplyr::filter(IDENT_SUBID == ids[indexer]))
      traindata <- traindata %>% dplyr::filter(IDENT_SUBID != ids[indexer])
      indexer <- indexer + 1
    }
    # Data cleaning of the results using a function defined in functions.R
    traindata <- traindata %>% fixup
    tunedata <- tunedata %>% fixup
    testdata <- testdata %>% fixup
    
    # Now, ensure that no one in tune is in train
    stopifnot(isFALSE(any(traindata$IDENT_SUBID %in% tunedata$IDENT_SUBID)))
    # All okay, then Rbind train and tune (Rbind is defined in functions.R)
    traindata <- Rbind(traindata, tunedata)
    # output data frames
    dataPartitions[[paste0('Resample_', i)]] <- list('train' = traindata,
                                                     'tunedata' = tunedata,
                                                     'trainfrac' = (nrow(traindata)-nrow(tunedata))/nrow(traindata),
                                                     'test' = testdata, 
                                                     'train.indices' = train.index)
    print(paste0("Partition ", i, " created."))
  }
  return(dataPartitions)
}

# FUNCTIONS FOR FITTING MODELS, CROSS-VALIDATION, VARIABLE IMPORTANCE COMPUTATION:
# 4. Function to conduct grid search cross-validation to tune the model at iteration n of 250
grid_search_CV <- function(parameter_grid, .ttSplit) {
  apply(parameter_grid, 1, function(params) {
    model <- try(gbm(
      GROUP ~.-IDENT_SUBID-brain_age_yrs-ICV-GENDER_FEMALE-subjWeight, 
      data = .ttSplit$train, # supply only the train data
      weights = subjWeight,
      distribution = "bernoulli", # model the binary outcome with a logit link fn
      n.trees = 750, # a large number to ensure the minimum criterion is reached
      shrinkage = params["lambda"],
      interaction.depth = params["depth"],
      n.minobsinnode = params["nodemin"],
      train.fraction = .ttSplit$trainfrac, # construct the gbm on the first train.fraction*nrow(data) rows
      bag.fraction = 1 # do not resample from the train set, we manually created a test set.
    ))
    if (class(model) == "try-error") {
      return(list(Tree = NA_integer_, TestAUC = -Inf))
    } # this will also end the anonymous function inside the apply call.
    # Otherwise, if the model fit successfully:
    # now for each tree in this one gbm model, get the AUC to our test set.
    tree.num <- gbm.perf(model, method = "test", plot.it = F)
    # If the best model is the last tree, then assume it hadn't yet converged
    # Add 16 trees at a time, and stop the model when the best model is at least 15 trees in the past
    # If the number of trees reaches 2000, then abort, and return the results of 2000th tree
    while((model$n.trees - tree.num) < 15) {
      model <- gbm.more(model, n.new.trees = 16)
      tree.num <- gbm.perf(model, method = "test", plot.it = F)
      if (model$n.trees >= 2000) {break}
    } # end of while loop
  
    # Now best is either the normally converging results, or the extended results from "while"
    # We need to return the best model, the AUC, as well as # of the best tree (not in that order)
    # Finally score on the TUNE set, which will later be used to determine the best grid search iteration
    tuneScore <- model$valid.error[tree.num]
    list(TheModel = model, BestTree = tree.num, TuneAUC = tuneScore)
  }) # ends the anonymous function inside the apply() and writes into "grid_search"
}

# 5. Convenience wrapper for predict() to get the right outputs for vi_permute
predict.gbm_bestTree <- function(object, newdata) 
{
  predict(object, newdata, type = "response", n.trees = gbm.perf(object, method = "test", plot.it = F))
}

# 6. conduct permutation and variable importance calculation for the *true* model
get_cvpvi = function(model_list, .ttSplit, .pred_wrapper) {
  vi_permute(
    object = model_list$TheModel, # supply the model
    feature_names = attr(model_list$TheModel$Terms, 'term.labels'), # just predictors
    train = .ttSplit$test, # supply the TEST data instead
    target = .ttSplit$test$GROUP, # test set response
    metric = "auc", # AUC-based importance
    reference_class = 1, # PI are coded 1 (since using "auc")
    type = "ratio", # PVI = auc/auc*
    nsim = 100, # permute each column 1000 times to compute its importance
    keep = TRUE, # KEEP ALL PERMUTATIONS
    sample_frac = 1, # use the whole test set
    pred_wrapper = .pred_wrapper # instruct each permutation to predict from this fn
  )
}
#7. Calculate the null model, null predictions, and null CVPVI
nullify = function(.ttSplit = ttSplit, 
                   model_list = THEBESTMODEL.list, 
                   parameter_grid=param_grid, 
                   .bestModNum=bestModNum,
                   .pred_wrapper = predict.gbm_bestTree,
                   .num_nulls = num_nulls) 
{
  purrr::map(setNames(as.list(1:.num_nulls), paste0("null_", 1:.num_nulls)), ~{
    nullTrain <- .ttSplit$train[1:round(nrow(.ttSplit$train)*.ttSplit$trainfrac), ] %>% mutate_at('GROUP', sample)
    nullTest <- .ttSplit$test %>% mutate_at('GROUP', sample)
    nullModel <- gbm(GROUP ~.-IDENT_SUBID-brain_age_yrs-ICV-GENDER_FEMALE-subjWeight, 
                     data = nullTrain, 
                     distribution = 'bernoulli',
                     n.trees = model_list$BestTree, 
                     shrinkage = parameter_grid$lambda[.bestModNum],
                     interaction.depth = parameter_grid$depth[.bestModNum],
                     n.minobsinnode = parameter_grid$nodemin[.bestModNum],
                     train.fraction = 1, bag.fraction = 1)
    # Before computing and saving the null var imps, also obtain the predictions
    nullModRes <- predict(nullModel, newdata=nullTest, type='response',n.trees=nullModel$n.trees)
    # Calculate the null variable importances
    nullPVI <- vi_permute(
      object = nullModel, # supply the model
      feature_names = attr(nullModel$Terms, "term.labels"), # just predictors
      train = nullTest, # supply the TEST data instead
      target = nullTest$GROUP, # test set response, nullified
      metric = "auc", # AUC-based importance
      reference_class = 1, # PI are coded 1 (since using "auc")
      type = "ratio", # PVI = auc/auc*
      nsim = 10, # permute each column just 10 times to compute its importance, b/c is null anyway
      keep = FALSE,
      sample_frac = 1, # use the whole test set
      pred_wrapper = .pred_wrapper # instruct each permutation to predict from this fn
    )[,"Importance"] # just return the null values
    return(list('y_hat'=nullModRes, 'nullPVI'=nullPVI))
  })
}
# 8. The master function
# just a function of the data, parameter grid, and size of null dist
gbmCV <- function(ttSplit, param_grid, num_nulls) 
{
  grid_search <- grid_search_CV(parameter_grid = param_grid, .ttSplit = ttSplit)
  # Now we have one model for each hyperparameter combination in "grid_search"
  # Grid Search CV is done. Now we simply delete all the worse models, and output the final model.
  # Identify the best parameter combination:
  bestModNum <- which.max(sapply(grid_search, function(x) {x$TuneAUC}))
  THEBESTMODEL.list <- grid_search[[bestModNum]]
  # Finally, compute the AUC for the test set
  TheTestAUC <- MLmetrics::AUC(predict(THEBESTMODEL.list$TheModel, newdata = ttSplit$test,
                                       type = "response", n.trees = THEBESTMODEL.list$BestTree),
                               y_true = ttSplit$test$GROUP)
  THEBESTMODEL.list$TestAUC <- TheTestAUC
  
  # Now that we know bestModNum is where the best model occurred, output it and related results
  # We also want to compute the variable importances and output them along with this model.
  # First we need to define a helper function
  # Variable importances
  pvi = get_cvpvi(model_list = THEBESTMODEL.list, .ttSplit = ttSplit, .pred_wrapper = predict.gbm_bestTree)
  # Null distribution of the model at each iteration.
  
  # Train a model with a permuted train response, compute variable importances to test data with permuted test response.
  null_res <- nullify(.ttSplit = ttSplit,  model_list = THEBESTMODEL.list,
                      parameter_grid=param_grid, .bestModNum=bestModNum,
                      .pred_wrapper = predict.gbm_bestTree,
                      .num_nulls = num_nulls)
  # Concatenate the 1000 predictions and varImps into N x 1000 (and p x 1000) matrices
  # Just fixing up the outputs from nullify before returning the results of gbmCV
  null_res <- list(
    'nullPreds'=reduce(lapply(null_res, function(x) x$y_hat), cbind) %>% as.matrix,
    'nullImp'=map_dfc(.x = lapply(null_res, function(x) x$nullPVI) %>% as.matrix, 
    .f = ~matrix(.x, ncol = 1))
  )
  
  # Finally return all outputs
  bestModelPkg <- list(
    "GBM"= list('BestModelTestSetPreds' = predict(THEBESTMODEL.list$TheModel, 
                                                  newdata = ttSplit$test, 
                                                  type='response', 
                                                  n.trees = THEBESTMODEL.list$BestTree),
                'BestModel_BestNumTrees' = THEBESTMODEL.list$BestTree,
                'BestModel_TestAUC' = THEBESTMODEL.list$TestAUC,
                # We also want to know what the best parameters were
                'BestParams' = as.data.frame(param_grid) %>% slice(bestModNum)),
    # And the variable importance info
    "PVI" = list('raw_scores'=attr(pvi, "raw_scores")),
    # Add the null information
    'null'=null_res
  )
  return(bestModelPkg) # one will be returned for each ttSplit
} # end of gbmCV function

# 9. an additional function simply to create tune and train data (no test) for ELFK pipeline
tuneSplits = function(n, .data, reproducibility_seed = 1)
{ 
  # Set the seed now once, before splitting data in a sequential for-loop below
  set.seed(reproducibility_seed)
  # A. Read in the supplied arguments to the function call, general setup
  data <- .data
  dataPartitions = vector(mode = 'list', length = 0)
  # Define the 'subjWeight' variable (indicating the .5 or 1.0 weights)
  #data <- data %>% group_by(IDENT_SUBID) %>% mutate(subjWeight = 1 / n()) %>% ungroup() 
  
  # B. The next three lines of code will:
  # 1. Begin by downsampling SB COMPs to create a balanced dataset
  # 2. Compute which indices of the data are the COMPs; then record their subject IDs
  # 3. Then generate samplings weights so that each SUBJECT (for the COMPs) is selected with the same probability for train
  comp.index = which(!as.logical(data$GROUP)) # converts 1 to F, reverses, sends to which()
  comp.subid = data$IDENT_SUBID[comp.index] # same for comps
  comp.weights = (1/n_distinct(comp.subid)) * data$subjWeight[comp.index] # define sampling weights
  # Nothing happens with the PI group, but their code is left here for posterity
  # One exception to the above is the pi.index variable, which is defined to later simply select all PI subjs
  pi.index = which(as.logical(data$GROUP)) # converts 1 to TRUE, passes along to which()
  #pi.subid = data$IDENT_SUBID[pi.index] # subject ids for the pi group
  #pi.weights = (1/n_distinct(pi.subid)) * data$subjWeight[pi.index] # define sampling weights
  
  # C. For each requested splits (the arg `n` to the function call), sample the comps into the train data
  ## Here, we define a balanced train dataset by selecting all PI, and an equal number of COMP scans
  for (i in 1:n) {
    # sample into a training set:
    train.index = c( 
      # create a sample with 70% of the number of people (subject ids)
      pi.index, # take all
      sample(comp.index, size = length(pi.index), prob = comp.weights) # downsample COMPs
    )
    traindata = data[sample(train.index), ] # shuffle the order
    
    # D. DEFINING THE TUNE SET
    # Do so by making a modification to the train data so that the first 70% of
    # the rows contain the true train subset that will be heldout for an internal tuning loop
    # by the gbm() function. To do so, we first select one PI subject, then a COMP, and we
    # continue to alternate so that the train and tune data remain as balanced as possible. Notably,
    # the COMP subjects that were not sampled into train in step C are dropped entirely.
    train_ids <- unique(traindata$IDENT_SUBID) # shuffle the order and take the train data subject IDs to guide selection
    train_pi.ids = traindata$IDENT_SUBID[which(as.logical(traindata$GROUP))] %>% unique
    train_comp.ids = traindata$IDENT_SUBID[which(!as.logical(traindata$GROUP))] %>% unique
    indexer <- 1
    tunedata <- traindata[0, ] # Create an empty data frame with the same column names
    while(nrow(tunedata)/nrow(traindata) < .33) {
      if ((indexer / 2) %>% as.character %>% grepl('\\.5$', x = .)) { # 'if the indexer is even not odd'
        tunedata <- Rbind(tunedata,
                          traindata %>% dplyr::filter(IDENT_SUBID == train_comp.ids[indexer]))
      } else {
        tunedata <- Rbind(tunedata,
                          traindata %>% dplyr::filter(IDENT_SUBID == train_pi.ids[indexer]))
      } # Delete the moved subject from the train data:
      traindata <- traindata %>% dplyr::filter(!IDENT_SUBID %in% tunedata$IDENT_SUBID)
      indexer <- indexer + 1
    }
    
    # E. Data cleaning of the results using a function defined in functions.R
    traindata <- traindata %>% fixup
    tunedata <- tunedata %>% fixup
    
    # Now, ensure that no one in tune is in train
    stopifnot(isFALSE(any(traindata$IDENT_SUBID %in% tunedata$IDENT_SUBID)))
    # All okay, then Rbind train and tune (Rbind is defined in functions.R)
    traindata <- Rbind(traindata, tunedata)
    # output data frames
    dataPartitions[[paste0('Resample_', i)]] <- list('train' = traindata,
                                                     'tunedata' = tunedata,
                                                     'trainfrac' = (nrow(traindata)-nrow(tunedata))/nrow(traindata),
                                                     'train.indices' = train.index)
    print(paste0("Partition ", i, " created."))
  }
  return(dataPartitions)
}

# 10. Get the features. Optionally also request and ID variable by name
get_X = function(data, id = NULL)
{
  data %>% 
    select(matches("^L\\.|^R\\.", ignore.case = F)) %>%
    # I despise tibbles:
    as.matrix %>% as.data.frame %>%
    # If requested, prepend the id (or any other) column too
    cbind(data[id], .)
}

# 11. Extract the labels
get_labels = function(data, append_col = NULL)
{
  data %>% 
    select(-matches("^L\\.|^R\\.|GROUP", ignore.case = F)) %>%
    # I despise tibbles:
    as.matrix %>% as.data.frame %>%
    # If requested, append the id column too
    cbind(., data[append_col])
}

# 12. IMPORTANT: Extract all columns needed for modeling
get_design = function(data)
{
  data %>%
    as.data.frame %>% # just in case
    select(IDENT_SUBID, subjWeight, GROUP, matches("^L\\.|^R\\.|GROUP", ignore.case = F))
}

# Wrapper to fit gbm model and get CVPVI
boost_trees = function(.ttSplit, .Test, 
                       # Set hyperparameters to gbm() defaults
                       .shrinkage = .1, 
                       .depth = 6, 
                       min_node_size=1, 
                       # provide option to fit null model
                       null = FALSE)
{
  mod = gbm(
    GROUP ~ .-subjWeight-IDENT_SUBID, 
    data = get_X(.ttSplit$train), # supply only the train data
    weights = subjWeight,
    distribution = "bernoulli", # model the binary outcome with a logit link fn
    n.trees = 750, # a large number to ensure the minimum criterion is reached
    shrinkage = .shrinkage,
    interaction.depth = .depth,
    n.minobsinnode = min_node_size,
    train.fraction = ttSplit[[1]]$trainfrac, # construct the gbm on the first train.fraction*nrow(data) rows
    bag.fraction = 1 # do not resample from the train set, we manually created a test set.
  )
  out = mclapply(1:num.iterations, function(i) {
    set.seed(i)
    split = tuneSplits(1, Train)
    .train = split[[1]]$train %>% select(-brain_age_yrs, -GENDER_FEMALE, -ICV) %>%
      mutate_at('IDENT_SUBID', as.character)
    mod = gbm(
      GROUP ~ .-subjWeight-IDENT_SUBID, 
      data = .train, # supply only the train data
      weights = subjWeight,
      distribution = "bernoulli", # model the binary outcome with a logit link fn
      n.trees = 750, # a large number to ensure the minimum criterion is reached
      shrinkage = .shrinkage,
      interaction.depth = .depth,
      n.minobsinnode = min_node_size,
      train.fraction = ttSplit[[1]]$trainfrac, # construct the gbm on the first train.fraction*nrow(data) rows
      bag.fraction = 1 # do not resample from the train set, we manually created a test set.
    )
    # Capture the 'best' tree according to tune error
    best_tree = mod %>% gbm.perf(method = 'test', plot = F)
    # Scores
    pred_proba = predict(mod, newdata = Test, n.trees = best_tree, type = 'response')
    pred_class = round(pred_proba)
    res = data.frame(
      iteration = i,
      tune_error = mod$valid.error[best_tree],
      test_err.rate = 1 - MLmetrics::Accuracy(pred_class, Test$GROUP),
      test_ROCAUC = MLmetrics::AUC(pred_proba, Test$GROUP),
      test_confusion = MLmetrics::ConfusionMatrix(pred_class, y_true=Test$GROUP)
    )
    res # returned
  }, mc.cores = 3)
  return(out)
}

get_gbm_results = function(model, Test)
{
  # Scores (requires predict.gbm_bestTree())
  pred_proba = predict.gbm_bestTree(model, newdata = get_design(Test))
  pred_class = round(pred_proba)
  predframe = data.frame(get_labels(Test, append_col = 'GROUP') %>% select(-subjWeight), 
                         pred_prob = pred_proba, pred_group = pred_class,
                         stringsAsFactors = F)
  res = data.frame(
    test_err.rate = 1 - MLmetrics::Accuracy(pred_class, Test$GROUP),
    test_ROCAUC = MLmetrics::AUC(pred_proba, Test$GROUP)
  )
  return(list('Predictions'=predframe, 'Scores'=res))
}

get_glm_results = function(model, Test)
{
   pred_proba = predict(model, newdata = Test %>% select(-one_of('IDENT_SUBID')), type = 'response')
   pred_class = round(pred_proba)
   predframe = data.frame(get_labels(Test, append_col = 'GROUP') %>% select(-subjWeight),
   			  pred_prob = pred_proba, pred_group = pred_class,
   			  stringsAsFactors = F)
   res = data.frame(
     test_err.rate = 1 - MLmetrics::AUC(pred_proba, Test$GROUP),
     test_ROCAUC = MLmetrics::AUC(pred_proba, Test$GROUP)
   )
   return(list('Predictions'=predframe, 'Scores'=res))
}

get_gbm_cvpvi = function(model, 
                         Test, 
                         .pred_wrapper = predict.gbm_bestTree,
                         null = FALSE) {
  # requires:
  # predict.gbm_bestTree <- function(object, newdata) 
  # {
  #   predict(object, newdata, type = "response", n.trees = gbm.perf(object, method = "test", plot.it = F))
  # }
  vi_permute(
    object = model, # supply the model
    feature_names = attr(model$Terms, 'term.labels'), # just predictors
    train = get_design(Test), # supply the TEST data instead
    target = Test$GROUP, # test set response
    metric = "auc", # AUC-based importance
    reference_class = 1, # PI are coded 1 (since using "auc")
    type = "ratio", # PVI = auc/auc*
    nsim = ifelse(!null, 100, 30),
    keep = FALSE, # Not needed
    sample_frac = 1, # use the whole test set
    pred_wrapper = .pred_wrapper # instruct each permutation to predict from this fn
  ) %>%
    as.data.frame %>% select(Variable, Importance) # drop off the standard deviation
}

get_glm_cvpvi = function(model, Test, null = FALSE)
{
   vi_permute(
     object = model,
     feature_names = names(get_X(Test)),
     train = get_design(Test) %>% select(-one_of('IDENT_SUBID')),
     target = Test$GROUP,
     metric = 'auc',
     reference_class = 1,
     type = 'ratio',
     nsim = ifelse(!null, 100, 30),
     keep = F,
     sample_frac = 1,
     pred_wrapper = function(object, newdata) predict(object, newdata %>% select(-one_of('IDENT_SUBID')), type = 'response')
   ) %>%
     as.data.frame %>% select(Variable, Importance)
}

# Repeat the above two for logistic regression
get_logreg_results = function(model, Test)
{
  # Scores (requires predict.gbm_bestTree())
  pred_proba = predict(model, newdata = get_X(Test, id=c('GROUP', 'subjWeight')), type = 'response')
  pred_class = round(pred_proba)
  predframe = data.frame(get_labels(Test, append_col = 'GROUP') %>% select(-subjWeight), 
                         pred_prob = pred_proba, pred_group = pred_class,
                         stringsAsFactors = F)
  res = data.frame(
    test_err.rate = 1 - MLmetrics::Accuracy(pred_class, Test$GROUP),
    test_ROCAUC = MLmetrics::AUC(pred_proba, Test$GROUP)
  )
  return(list('Predictions'=predframe, 'Scores'=res))
}

get_logreg_cvpvi = function(model, 
                         Test, 
                         .pred_wrapper = function(object, newdata) predict(object, 
                                                                           newdata, 
                                                                           type = 'response'),
                         null = FALSE) {
  vi_permute(
    object = model, # supply the model
    feature_names = names(get_X(Test)), # just predictors
    train = get_X(Test, id = c('GROUP', 'subjWeight')), # supply the TEST data instead
    target = Test$GROUP, # test set response
    metric = "auc", # AUC-based importance
    reference_class = 1, # PI are coded 1 (since using "auc")
    type = "ratio", # PVI = auc/auc*
    nsim = ifelse(!null, 100, 30),
    keep = FALSE, # Not needed
    sample_frac = 1, # use the whole test set
    pred_wrapper = .pred_wrapper # instruct each permutation to predict from this fn
  ) %>%
    as.data.frame %>% select(Variable, Importance) # drop off the standard deviation
}

Residualize = function(col_name, train) {
  # USAGE: #imap(DATA[grep("^R\\.|^L\\.", names(DATA))], ~Residualize(.y, DATA))
  # This function is for residualizing brain regions with respect to regressors
  # First, create a subject weight column based on the number of times the subject appears in the data
  # create the weight variable here:
  dat <- train
  # Check data class for the sex factor
  if (!is.factor(dat$GENDER_FEMALE)) {dat$GENDER_FEMALE <- factor(dat$GENDER_FEMALE)}
  # rename the predictor to-be-adjusted to 'x' for uniformity below
  model.data <- dat[c(col_name, "subjWeight", 'brain_age_yrs','GENDER_FEMALE','ICV')]
  names(model.data)[1] <- 'COL' # for uniformity
  # Test for a significant correlation with the brain region in question
  ageSig = cor.test(model.data$COL, model.data$brain_age_yrs)$p.value < .05
  anovaSig = summary(lm(COL ~ GENDER_FEMALE, data = model.data, weights = subjWeight))$coef["GENDER_FEMALE1", "Pr(>|t|)"] < .05
  icvSig = cor.test(model.data$COL, model.data$ICV)$p.value < .05
  # collect these results
  allSigs <- setNames(c(ageSig, anovaSig, icvSig), nm=c('brain_age_yrs','GENDER_FEMALE','ICV')) # define in same order as the columns
  
  # If any significant correlation is detected, proceed with regression to residualize (if linear mod is significant):
  # define Pvalues to mirror allSigs, but set it first (here) to FALSE for both variables
  lm_model <- NULL
  Pvalues <- FALSE
  lmVarNameSig <- character(0)
  if (any(allSigs)) {
    # Fitting the linear model and extracting residuals.
    
    # This will be done in a while loop, and re-run automatically if a term in the model is not sig at p < .05
    # That term will be dropped prior to re-fitting.
    # If all terms are non-signiicant, the model is not fit, the variable is not adjusted by these covariates.
    # When the final model with only significant terms is found the residuals are extracted.
    while(isFALSE(all(Pvalues))) { # if any are false, then we continue with model fitting
      # Full model starting point (and there after subsetted by an updated allSigs):
      input_data <- model.data[, c("COL", "subjWeight", c("brain_age_yrs", "GENDER_FEMALE", "ICV")[allSigs])]
      
      # Fit the model to the surviving predictors:
      lm_model <- lm(formula = COL ~ .-subjWeight, data = input_data, weights = subjWeight)
      mod.sum <- summary(lm_model)
      # Record the names of the significant variables:
      predictors <- Filter(function(x) x %in% c("brain_age_yrs", "GENDER_FEMALE", "ICV"), names(lm_model$model))
      Pvalues <- mod.sum$coef[, ncol(mod.sum$coef)][-1] < .05 # minus intercept
      lmVarNameSig <- unique(c(lmVarNameSig, predictors[Pvalues]))
      # update the list that has boolean entry for each variable in the full covariate set:
      allSigs <- c("brain_age_yrs", "GENDER_FEMALE", "ICV") %in% lmVarNameSig
      # Thus, the loop will end if all variables are significant, so that isFALSE(all(S.PV, P.PV)) == FALSE,
      # or if no variables are entered into the model at all. In this case S.PV and P.PV will both be equal to
      # `logical(0)`, which means that all(logical(0), logical(0)) is actually/ironically TRUE, ending the loop.
      # Otherwise, if anything is FALSE, the while loop will continue (and that variable won't be in the next formula)
    }# end of while loop.
    allSigs<-setNames(allSigs, nm=c("brain_age_yrs", "GENDER_FEMALE", "ICV"))
    
    # Finally, if there were any significant predictors in the model, only they will have been retained.
    # Check again with S.PV and P.PV and then extract the residuals if so.
    if (any(Pvalues)) { # here, any(logical(0)) = FALSE...perfect! (if `any()` is "any are TRUE")
      train.resids <- lm_model$residuals
    } else { # otherwise just return the normalized columns
      train.resids <- as.vector(model.data$COL)
    }
    # end of if (any(allSigs)), so add an else for the case where nothing was sig:
  } else {
    # In this scenario also return the normalized column
    train.resids <- as.vector(model.data$COL)
  }
  
  # Write out the z-scored residuals:
  #return(NULL)
  return(as.vector(unclass(scale(train.resids))))
}

Residualize_noICV = function(col_name, train) {
  # USAGE: #imap(DATA[grep("^R\\.|^L\\.", names(DATA))], ~Residualize(.y, DATA))
  # This function is for residualizing brain regions with respect to regressors
  # First, create a subject weight column based on the number of times the subject appears in the data
  # create the weight variable here:
  dat <- train
  # Check data class for the sex factor
  if (!is.factor(dat$GENDER_FEMALE)) {dat$GENDER_FEMALE <- factor(dat$GENDER_FEMALE)}
  # rename the predictor to-be-adjusted to 'x' for uniformity below
  model.data <- dat[c(col_name, "subjWeight", 'brain_age_yrs','GENDER_FEMALE')]
  names(model.data)[1] <- 'COL' # for uniformity
  # Test for a significant correlation with the brain region in question
  ageSig = cor.test(model.data$COL, model.data$brain_age_yrs)$p.value < .05
  anovaSig = summary(lm(COL ~ GENDER_FEMALE, data = model.data, weights = subjWeight))$coef["GENDER_FEMALE1", "Pr(>|t|)"] < .05
  # collect these results
  allSigs <- setNames(c(ageSig, anovaSig), nm=c('brain_age_yrs','GENDER_FEMALE')) # define in same order as the columns
  
  # If any significant correlation is detected, proceed with regression to residualize (if linear mod is significant):
  # define Pvalues to mirror allSigs, but set it first (here) to FALSE for both variables
  lm_model <- NULL
  Pvalues <- FALSE
  lmVarNameSig <- character(0)
  if (any(allSigs)) {
    # Fitting the linear model and extracting residuals.
    
    # This will be done in a while loop, and re-run automatically if a term in the model is not sig at p < .05
    # That term will be dropped prior to re-fitting.
    # If all terms are non-signiicant, the model is not fit, the variable is not adjusted by these covariates.
    # When the final model with only significant terms is found the residuals are extracted.
    while(isFALSE(all(Pvalues))) { # if any are false, then we continue with model fitting
      # Full model starting point (and there after subsetted by an updated allSigs):
      input_data <- model.data[, c("COL", "subjWeight", c("brain_age_yrs", "GENDER_FEMALE")[allSigs])]
      
      # Fit the model to the surviving predictors:
      lm_model <- lm(formula = COL ~ .-subjWeight, data = input_data, weights = subjWeight)
      mod.sum <- summary(lm_model)
      # Record the names of the significant variables:
      predictors <- Filter(function(x) x %in% c("brain_age_yrs", "GENDER_FEMALE"), names(lm_model$model))
      Pvalues <- mod.sum$coef[, ncol(mod.sum$coef)][-1] < .05 # minus intercept
      lmVarNameSig <- unique(c(lmVarNameSig, predictors[Pvalues]))
      # update the list that has boolean entry for each variable in the full covariate set:
      allSigs <- c("brain_age_yrs", "GENDER_FEMALE") %in% lmVarNameSig
      # Thus, the loop will end if all variables are significant, so that isFALSE(all(S.PV, P.PV)) == FALSE,
      # or if no variables are entered into the model at all. In this case S.PV and P.PV will both be equal to
      # `logical(0)`, which means that all(logical(0), logical(0)) is actually/ironically TRUE, ending the loop.
      # Otherwise, if anything is FALSE, the while loop will continue (and that variable won't be in the next formula)
    }# end of while loop.
    allSigs<-setNames(allSigs, nm=c("brain_age_yrs", "GENDER_FEMALE"))
    
    # Finally, if there were any significant predictors in the model, only they will have been retained.
    # Check again with S.PV and P.PV and then extract the residuals if so.
    if (any(Pvalues)) { # here, any(logical(0)) = FALSE...perfect! (if `any()` is "any are TRUE")
      train.resids <- lm_model$residuals
    } else { # otherwise just return the normalized columns
      train.resids <- as.vector(model.data$COL)
    }
    # end of if (any(allSigs)), so add an else for the case where nothing was sig:
  } else {
    # In this scenario also return the normalized column
    train.resids <- as.vector(model.data$COL)
  }
  
  # Write out the z-scored residuals:
  #return(NULL)
  return(as.vector(unclass(scale(train.resids))))
}

Residualize_onlyICV = function(col_name, train) {
  # USAGE: #imap(DATA[grep("^R\\.|^L\\.", names(DATA))], ~Residualize(.y, DATA))
  # This function is for residualizing brain regions with respect to regressors
  # First, create a subject weight column based on the number of times the subject appears in the data
  # create the weight variable here:
  dat <- train
  # Check data class for the sex factor
  if (!is.factor(dat$GENDER_FEMALE)) {dat$GENDER_FEMALE <- factor(dat$GENDER_FEMALE)}
  # rename the predictor to-be-adjusted to 'x' for uniformity below
  model.data <- dat[c(col_name, "subjWeight", 'ICV')]
  names(model.data)[1] <- 'COL' # for uniformity
  # Test for a significant correlation with the brain region in question
  #ageSig = cor.test(model.data$COL, model.data$brain_age_yrs)$p.value < .05
  #anovaSig = summary(lm(COL ~ GENDER_FEMALE, data = model.data, weights = subjWeight))$coef["GENDER_FEMALE1", "Pr(>|t|)"] < .05
  icvSig = cor.test(model.data$COL, model.data$ICV)$p.value < .05
  # collect these results
  allSigs <- setNames(c(icvSig), nm=c('ICV')) # define in same order as the columns
  
  # If any significant correlation is detected, proceed with regression to residualize (if linear mod is significant):
  # define Pvalues to mirror allSigs, but set it first (here) to FALSE for both variables
  lm_model <- NULL
  Pvalues <- FALSE
  lmVarNameSig <- character(0)
  if (any(allSigs)) {
    # Fitting the linear model and extracting residuals.
    
    # This will be done in a while loop, and re-run automatically if a term in the model is not sig at p < .05
    # That term will be dropped prior to re-fitting.
    # If all terms are non-signiicant, the model is not fit, the variable is not adjusted by these covariates.
    # When the final model with only significant terms is found the residuals are extracted.
    while(isFALSE(all(Pvalues))) { # if any are false, then we continue with model fitting
      # Full model starting point (and there after subsetted by an updated allSigs):
      input_data <- model.data[, c("COL", "subjWeight", c("ICV")[allSigs])]
      
      # Fit the model to the surviving predictors:
      lm_model <- lm(formula = COL ~ .-subjWeight, data = input_data, weights = subjWeight)
      mod.sum <- summary(lm_model)
      # Record the names of the significant variables:
      predictors <- Filter(function(x) x %in% c("ICV"), names(lm_model$model))
      Pvalues <- mod.sum$coef[, ncol(mod.sum$coef)][-1] < .05 # minus intercept
      lmVarNameSig <- unique(c(lmVarNameSig, predictors[Pvalues]))
      # update the list that has boolean entry for each variable in the full covariate set:
      allSigs <- c("ICV") %in% lmVarNameSig
      # Thus, the loop will end if all variables are significant, so that isFALSE(all(S.PV, P.PV)) == FALSE,
      # or if no variables are entered into the model at all. In this case S.PV and P.PV will both be equal to
      # `logical(0)`, which means that all(logical(0), logical(0)) is actually/ironically TRUE, ending the loop.
      # Otherwise, if anything is FALSE, the while loop will continue (and that variable won't be in the next formula)
    }# end of while loop.
    allSigs<-setNames(allSigs, nm=c("ICV"))
    
    # Finally, if there were any significant predictors in the model, only they will have been retained.
    # Check again with S.PV and P.PV and then extract the residuals if so.
    if (any(Pvalues)) { # here, any(logical(0)) = FALSE...perfect! (if `any()` is "any are TRUE")
      train.resids <- lm_model$residuals
    } else { # otherwise just return the normalized columns
      train.resids <- as.vector(model.data$COL)
    }
    # end of if (any(allSigs)), so add an else for the case where nothing was sig:
  } else {
    # In this scenario also return the normalized column
    train.resids <- as.vector(model.data$COL)
  }
  
  # Write out the z-scored residuals:
  #return(NULL)
  return(as.vector(unclass(scale(train.resids))))
}

# This CVPVI function generates the predictions of the
CVPVI = function(Model, Test, Permutations, FeatureNames, pred_func) 
{
  mclapply(X=FeatureNames, mc.cores = 3, FUN=function(x) {
    rowMeans(sapply(as.data.frame(Permutations), function (p) {
      new_dat <- Test %>% mutate_at(x, ~.[p])
      pred_func(Model, new_dat)
    }))
  }) %>% reduce(cbind) %>% as.data.frame %>% setNames(nm = FeatureNames)
}