#!/usr/bin/env Rscript

# Estimate allele lengths and find outliers at STR loci

suppressPackageStartupMessages({
  library('optparse', quietly = TRUE)
})

option_list = list(
  make_option("--dir",
              help = paste("Working directory containing STR/locus count data and median coverage",
                           "[default: %default]"),
              type = "character",
              default = '.'),
  make_option("--out",
              help = paste("Prefix for all output files",
                           "[default: %default]"),
              type = "character",
              default = ''),
  make_option("--model",
              help = paste("Data to produce linear regression model (provided with STRetch)",
                           "[default: %default]"),
              type = "character",
              default = 'STRcov.model.csv')
)

parser = OptionParser(
  usage = paste("%prog [OPTIONS]",
                "Estimate allele lengths and find outliers at STR loci.\n",
                sep = "\n"),
  option_list = option_list
)

# Parse command line arguments
INVALID_MESSAGE <- "Invalid command line arguments. Use --help for help."
tryCatch({
  suppressWarnings(
    arguments <- parse_args(parser)
  )},
  error = function(e) {
    message(INVALID_MESSAGE)
    quit(status = 2)
  }
)

## Main
base.filename = arguments$out
data.dir = arguments$dir
STRcov.model.csv = arguments$model

locuscov.files = list.files(data.dir, 'locus_counts.txt$', full.names = TRUE)
STRcov.files = list.files(data.dir, 'STRcoverage.txt$', full.names = TRUE)
genomecov.files = list.files(data.dir, 'median_cov$', full.names = TRUE)
#metadata = read.delim(paste(data.dir, "/ataxia_wgs_metadata.txt", sep = ""), stringsAsFactors=FALSE) #XXX make this an optional input?

# Only load packages once files have been read, to save time during program startup when debugging.
suppressPackageStartupMessages({
  library('plyr', quietly = TRUE)
  library('dplyr', quietly = TRUE)
  library('tidyr', quietly = TRUE)
  library('reshape2', quietly = TRUE)
  library('magrittr', quietly = TRUE)
})

get.sample = function(filename) {
  dir.split = tail(strsplit(filename, '/', fixed = TRUE)[[1]], n = 1)
  strsplit(dir.split, '(\\.|_)', fixed = FALSE)[[1]][1]
}

# Parse all STR coverage
parse.STRcov = function(filename) {
  sim.id = get.sample(filename)
  
  cov.data = read.table(filename, stringsAsFactors = FALSE, col.names = c('chrom', 'start', 'end', 'coverage'))
  cov.data$sample = sim.id
  cov.data$repeatunit = sapply(cov.data$chrom, function(x){ strsplit(x, '-', fixed = TRUE)[[1]][2] } )
  cov.data = cov.data[,c('sample', 'repeatunit', 'coverage')]
  
  return(cov.data)
}

# Parse all locus coverage
parse.locuscov = function(filename) {
  sim.id = get.sample(filename)
  locuscov.data = read.table(filename, stringsAsFactors = FALSE, header=TRUE)
  locuscov.data$sample = sim.id
  locuscov.data$locus = paste(locuscov.data$STR_chr, locuscov.data$STR_start, locuscov.data$STR_stop, sep='-')
  locuscov.data$repeatunit = locuscov.data$motif
  locuscov.data$locuscoverage = locuscov.data$count
  return(locuscov.data[,c('sample', 'locus', 'repeatunit', 'locuscoverage')])
}

# Assumes median coverage is the top left value in the text file
parse.genomecov = function(filename) {
  sim.id = get.sample(filename)
  mediancov = read.table(filename)[1,1]
  stopifnot(is.numeric(mediancov))
  
  genomecov.data = data.frame(sample = sim.id, genomecov = mediancov)
  return(genomecov.data[,c('sample','genomecov')])
}

get.genomecov = function(sample, database) {
  database$genomecov[which(database$sample == sample)]
}

# Look for repeat units where the number of reads mapping to the decoy can't be 
# explained by those mapping to all loci with that repeat unit
count.differences = function(locuscov.df, STRcov.df){
  # Sum the counts over all loci for each repeat unit
  locus.totals = group_by(locuscov.df, repeatunit) %>% dplyr::summarise(totalSTRcov = sum(locuscoverage))
  # Merge these with the total counts mapping to the STR decoy
  STRcov.df.totals = merge(STRcov.df, locus.totals, by="repeatunit", all.y = TRUE)
  # Take the difference to get the number of reads mapping to the STR decoy that aren't assigned to a locus
  STRcov.df.totals$difference = STRcov.df.totals$coverage - STRcov.df.totals$totalSTRcov
  return(STRcov.df.totals)
}

## Parse coverage
# parse.STRcov
STRcov.data = ldply(lapply(STRcov.files, parse.STRcov), data.frame)
STRcov.data$coverage = as.numeric(STRcov.data$coverage)
# parse.locuscov
locuscov.data = ldply(lapply(locuscov.files, parse.locuscov), data.frame)
locuscov.data$locuscoverage = as.numeric(locuscov.data$locuscoverage)
# Fill zeros in locuscov
locuscov.data.wide = spread(locuscov.data, sample, locuscoverage)
locuscov.data.wide[is.na(locuscov.data.wide)] = 0
sample.cols = which(names(locuscov.data.wide) %in% unique(locuscov.data$sample))
locuscov.data = gather(locuscov.data.wide, sample, locuscoverage, sample.cols)

# Normalise by median coverage
genomecov.data = ldply(lapply(genomecov.files, parse.genomecov), data.frame)
genomecov.data$genomecov = as.numeric(genomecov.data$genomecov)
# Check coverage is available for all samples
stopifnot(unique(STRcov.data$sample) %in% genomecov.data$sample)

#genomecov.data = group_by(genomecov.data, sample) %>% summarise(genomecov = sum(genomecov)) # future feature? only requred if multiple files per sample
#genomecov.data = genomecov.data[genomecov.data$sample %in% male.affected, ]
factor = 100

STRcov.data$coverage_norm = factor * (STRcov.data$coverage + 1) / sapply(STRcov.data$sample, get.genomecov, genomecov.data)
STRcov.data$coverage_log = log2(STRcov.data$coverage_norm)
locuscov.data$locuscoverage_norm = factor * (locuscov.data$locuscoverage + 1) / sapply(locuscov.data$sample, get.genomecov, genomecov.data)
locuscov.data$locuscoverage_log = log2(locuscov.data$locuscoverage_norm)

locuscov.data$locusSTR = paste(locuscov.data$repeatunit, locuscov.data$locus, sep = ' ')

all.differences = ungroup(group_by(locuscov.data, sample) %>% do({count.differences(., STRcov.data[STRcov.data$sample == .$sample[1],])}))
# Normalise by median coverage
all.differences$difference_norm = factor * (all.differences$difference + 1) / sapply(all.differences$sample, get.genomecov, genomecov.data)
all.differences$difference_log = suppressWarnings(
  log2(all.differences$difference_norm)
  )

# Assign decoy counts to each locus, based on what proportion of the counts for that repeat unit they already have
# Get total counts assigned to all loci of a given repeat unit (i.e. see above)
# Calculate the proportion from each locus.
# Group by sample, then calculate proportion reads to assign to each locus
STRcovtotals = group_by(locuscov.data, sample, repeatunit) %>% dplyr::summarise(totalSTRcov = sum(locuscoverage))
locuscov.totals = merge(locuscov.data, STRcovtotals)
locuscov.totals$locuscoverage_prop = locuscov.totals$locuscoverage/locuscov.totals$totalSTRcov
# Assign that proportion of the total decoy counts for that repeat unit to each locus.
locuscov.totals = merge(locuscov.totals, all.differences)
#XXX Note: this is assigning a proportion of the difference counts, not of the original STR decoy counts.
locuscov.totals$diff_assigned = locuscov.totals$locuscoverage_prop * locuscov.totals$difference
locuscov.totals$total_assigned = locuscov.totals$diff_assigned + locuscov.totals$locuscoverage


# For each locus, calculate if that sample is an outlier relative to the others

# Normalise by median coverage
locuscov.totals$total_assigned_norm = factor * (locuscov.totals$total_assigned + 1) / sapply(locuscov.totals$sample, get.genomecov, genomecov.data)
locuscov.totals$total_assigned_log = log2(locuscov.totals$total_assigned_norm)

total_assigned_wide = acast(locuscov.totals, locus ~ sample, value.var = "total_assigned_log")
locuscoverage_log_wide = acast(locuscov.totals, locus ~ sample, value.var = "locuscoverage_log")
## Like a z score, except using the median and IQR instead of mean and sd.
#XXX change this to total_assigned_wide
z = apply(locuscoverage_log_wide, 1, function(x) {(x - median(x)) / IQR(x)})
#z = apply(total_assigned_wide, 1, function(x) {(x - median(x)) / IQR(x)})
z.long = melt(z, varnames = c('sample', 'locus'), value.name = 'outlier')
locuscov.totals = merge(locuscov.totals, z.long)

# Predict size (in bp) using the ATXN8 linear model (produced from data in decoySTR_cov_sim_ATXN8_AGC.R) 
# Read in the raw data for this model from a file
STRcov.model = read.csv(STRcov.model.csv) 
# Model is build from log2 data (to reduce heteroscedasticity), then converted back
# *2 to convert from repeat units to bp
ATXN8.lm.data = data.frame(allele = log2(STRcov.model$allele2), coverage = log2(STRcov.model$coverage_norm))
ATXN8.lm = lm(data = ATXN8.lm.data, formula = allele ~ coverage)
predict = 2^predict(ATXN8.lm, data.frame(coverage = locuscov.totals$total_assigned_log), interval="confidence")
locuscov.totals = cbind(locuscov.totals, predict)

# Write all samples to a single file
write.data = locuscov.totals[,c('sample', 'locus', 'repeatunit',
                                'totalSTRcov', 'locuscoverage', 'total_assigned', 'outlier', 'fit', 'lwr', 'upr')]
write.data = write.data[order(write.data$outlier, decreasing = T),] #sort by outlier score
write.csv(x = write.data, file = paste(c(base.filename, 'locuscov.totals.csv'), collapse=''), quote = FALSE, row.names=FALSE)

# Write individual files for each sample, remove rows where locuscoverage == 0
samples = unique(write.data$sample)
for (sample in samples) {
  write.csv(write.data[write.data$sample == sample & write.data$locuscoverage != 0,], 
            file = paste(c(base.filename, '_', sample, 'locuscov.totals.csv'), collapse =''), 
            quote = FALSE, row.names=FALSE)
}
