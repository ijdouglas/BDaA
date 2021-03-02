# This script splits the data a specified number of times into train and tune sets.
# It creates balanced train and tune sets with respect to the grouping variable;
# (i.e.) group size imbalance is NOT preserved (i.e. this is not stratification).

# The result is a spreadsheet in which each column contains the row number of the
# subject that has been sampled into the training data. The row numbers range from
# 1 to 226, pertaining to the 226 subjects in the full training sample (SB). Each
# row contains 180 of the 226 indices, because there were 90 subjects in the smaller
# group. Each column contains a unique combination of row indices. 
# The column title (variable name) of each column encodes the fraction of these 180
# indices that pertain to the TRAIN data, while the remainder pertain to the TUNE
# data, as follows: "trainfrac[###]", where ### is replaced with said fraction.

# This script takes three arguments (in this order)
# 1. The full training sample file path (the SB data; relative to working dir)
# 2. The input file path/name (will be relative to working dir by default)
# 3. The number of partitions to create

# Read in the file that contains the function which splits the data as above
source('functions.R')
## Loading packages
suppressMessages(library(tidyverse))
# Read in the command line arguments (details above)
input_file = as.character(commandArgs(trailingOnly = T)[1])
output_file = as.character(commandArgs(trailingOnly = T)[2])
n_splits = as.numeric(commandArgs(trailingOnly = T)[3])

# Figure out the file type of the training data to read it in accordingly
# This is not exhaustive!
ext <- substr(x = input_file, start = nchar(input_file)-2, stop = nchar(input_file))
if (grepl("rds", ext, ignore.case = T)) {
  read_func <- function(x) readRDS(x)
} else if (grepl("csv", ext, ignore.case = T)) {
  read_func <- function(x) read.csv(x, stringsAsFactors = F)
}
# Read in the training data
TRAIN <- read_func(input_file)

# Split the data
splits <- tuneSplits(n = n_splits, .data = TRAIN, reproducibility_seed = 1)

# Convert the results to a csv interpretable as described above
col_names <- sapply(splits, function(x) paste0('trainfrac[', x$trainfrac, ']')) %>%
  as.vector
out <- as.data.frame(do.call('cbind', args = map_dfc(splits, ~.$train.indices)))
names(out) <- col_names

# Write it out
ext <- substr(x = output_file, start = nchar(output_file)-2, stop = nchar(output_file))
if (grepl("rds", ext, ignore.case = T)) {
  write_func <- function(.object, x) saveRDS(.object, x)
} else if (grepl("csv", ext, ignore.case = T)) {
  write_func <- function(.object, x) write.csv(.object, x, row.names = F)
}

write_func(out, output_file)