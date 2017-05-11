
# Estimate allele lengths and find outliers at STR loci

# Print traceback on failure (for debugging)
#options(error=traceback)
options(stringsAsFactors=FALSE)
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
              help = paste("Prefix for all output files (suffix will be STRs.tsv)",
                           "[default: %default]"),
              type = "character",
              default = ''),
  make_option("--model",
              help = paste("Data to produce linear regression model (provided with STRetch)",
                           "[default: %default]"),
              type = "character",
              default = 'STRcov.model.csv'),
  make_option("--control",
              help = paste("Input file for median and standard deviation estimates at each ",
                           "locus from a set of control samples. This file can be produced by ",
                           "this script using the emit option. If this option is not set, all ", 
                           "samples in the current batch will be used as controls by default."),
              type = "character",
              default = ''),
  make_option("--emit",
              help = "Output file for median and standard deviation estimates at each locus (tsv).",
              type = "character",
              default = '')
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
emit.file = arguments$emit
control.file = arguments$control

#XXX Change this so input files are provided as arguments, and check the number of files
locuscov.files = list.files(data.dir, 'locus_counts$', full.names = TRUE)
STRcov.files = list.files(data.dir, 'STR_counts$', full.names = TRUE)
genomecov.files = list.files(data.dir, 'median_cov$', full.names = TRUE)
results.suffix = 'STRs.tsv'

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
  return(locuscov.data[,c('sample', 'locus', 'repeatunit', 'reflen', 'locuscoverage')])
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

# Calculate a z score, except using Huber's M-estimator to estimate median and SD
hubers.z = function(x, k=1.5) {
  hubers.est = MASS::hubers(x, k)
  (x - hubers.est$mu)/sqrt(hubers.est$s)
}

# Emit Huber's M-estimator median and SD estimates for use as a control set
hubers.est = function(x, k=1.5) {
  hubers.est = MASS::hubers(x, k)
  return(c(median = hubers.est$mu, SD = sqrt(hubers.est$s)))
}

## Parse coverage
# parse.STRcov
STRcov.data = ldply(lapply(STRcov.files, parse.STRcov), data.frame)
STRcov.data$coverage = as.numeric(STRcov.data$coverage)
# parse.locuscov
locuscov.data = ldply(lapply(locuscov.files, parse.locuscov), data.frame)
locuscov.data$locuscoverage = as.numeric(locuscov.data$locuscoverage)
# Check for multiple rows with the same sample/locus combination and throw an error if found
crosstab = table(locuscov.data$locus, locuscov.data$sample)
multi.loci = row.names(crosstab)[apply(crosstab, 1, function(x) {any(x>1)})]
if(length(multi.loci) > 0) {
    stop('The locus count input data contains multiple rows with the same sample/locus combination. ',
      'This is usually caused by two loci at the same position in the STR annotation bed file. ', 
      'Check these loci:\n', paste(multi.loci, collapse = '\n'))
}

# Fill zeros in locuscov
locuscov.data.wide = spread(locuscov.data, sample, locuscoverage, fill = 0)
sample.cols = which(names(locuscov.data.wide) %in% unique(locuscov.data$sample))
locuscov.data = gather(locuscov.data.wide, sample, locuscoverage, sample.cols)

# Normalise by median coverage
genomecov.data = ldply(lapply(genomecov.files, parse.genomecov), data.frame)
genomecov.data$genomecov = as.numeric(genomecov.data$genomecov)
# Check coverage is available for all samples
stopifnot(unique(STRcov.data$sample) %in% genomecov.data$sample)

factor = 100

STRcov.data$coverage_norm = factor * (STRcov.data$coverage + 1) / sapply(STRcov.data$sample, get.genomecov, genomecov.data)
STRcov.data$coverage_log = log2(STRcov.data$coverage_norm)
locuscov.data$locuscoverage_norm = factor * (locuscov.data$locuscoverage + 1) / sapply(locuscov.data$sample, get.genomecov, genomecov.data)
locuscov.data$locuscoverage_log = log2(locuscov.data$locuscoverage_norm)

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
locuscov.totals$locuscoverage_prop[is.nan(locuscov.totals$locuscoverage_prop)] = 0 # set to 0 when divide by 0 error
# Assign that proportion of the total decoy counts for that repeat unit to each locus.
locuscov.totals = merge(locuscov.totals, all.differences)
# Note: this is assigning a proportion of the difference counts, not of the original STR decoy counts.
locuscov.totals$diff_assigned = locuscov.totals$locuscoverage_prop * locuscov.totals$difference
locuscov.totals$total_assigned = locuscov.totals$diff_assigned + locuscov.totals$locuscoverage

# For each locus, calculate if that sample is an outlier relative to the others

# Normalise by median coverage
locuscov.totals$total_assigned_norm = factor * (locuscov.totals$total_assigned + 1) / sapply(locuscov.totals$sample, get.genomecov, genomecov.data)
locuscov.totals$total_assigned_log = log2(locuscov.totals$total_assigned_norm)
total_assigned_wide = acast(locuscov.totals, locus ~ sample, value.var = "total_assigned_log")

if (control.file != '') {
  # Calculate z scores using median and SD estimates per locus from a provided control set
  control.estimates = read.table(control.file)
  # Reorder the control estimates by the loci in the sample
  control.estimates.sorted = control.estimates[match(rownames(total_assigned_wide), rownames(control.estimates)),]
  # Replace all loci not in the control with the default values
  control.estimates.sorted[rowSums(is.na(control.estimates.sorted)) == ncol(control.estimates.sorted),] = control.estimates["null_locus_counts",]
  # Calculate z scores
  z = t(apply(total_assigned_wide, 2, function(x){(x - control.estimates.sorted$median)/control.estimates.sorted$SD}))
  #z.long = melt(z, varnames = c('locus', 'sample'), value.name = 'outlier')
} else {
  # Calculate a z score using Huber's M-estimator to calculate median and SD across all samples
  z = apply(total_assigned_wide, 1, hubers.z)
  
}
z.long = melt(z, varnames = c('sample', 'locus'), value.name = 'outlier')
locuscov.totals = merge(locuscov.totals, z.long)

# Calculate p values based on z scores (one sided)
pvals = apply(z, 1, function(x) {pnorm(x, lower.tail=FALSE)})
adj_pvals = apply(pvals, 2, function(x) {p.adjust(x, method = "BH")})
pvals.long = reshape2::melt(adj_pvals, varnames = c('locus', 'sample'), value.name = 'p_adj')
locuscov.totals = merge(locuscov.totals, pvals.long)

# Save median and SD of all loci to file if requested (for use as a control set for future data sets)
if (emit.file != '') {
  # Add a null locus that has 0 reads for all individuals (so just uses coverage)
  null_locus_counts = log2(factor * (0 + 1)/genomecov.data$genomecov)
  
  hubers.estimates = data.frame(t(apply(rbind(null_locus_counts, total_assigned_wide), 1, hubers.est)))
  hubers.estimates$n = ncol(total_assigned_wide)
  
  write.table(hubers.estimates, file=emit.file, sep = "\t", quote = FALSE)
}

# Predict size (in bp) using the ATXN8 linear model (produced from data in decoySTR_cov_sim_ATXN8_AGC.R) 
# Read in the raw data for this model from a file
# Note: coverage_norm = (STR coverage/median coverage) * 100
# allele2 is the length of the longer allele in bp inserted relative to reference
STRcov.model = read.csv(STRcov.model.csv) 
# Model is build from log2 data (to reduce heteroscedasticity), then converted back
ATXN8.lm.data = data.frame(allele = log2(STRcov.model$allele2), coverage = log2(STRcov.model$coverage_norm))
ATXN8.lm = lm(data = ATXN8.lm.data, formula = allele ~ coverage)
predict = 2^predict(ATXN8.lm, data.frame(coverage = locuscov.totals$total_assigned_log), interval="confidence")
locuscov.totals = cbind(locuscov.totals, predict)

# Get the estimated size in terms of repeat units (total, not relative to ref)
locuscov.totals$repeatUnits = locuscov.totals$fit/nchar(locuscov.totals$repeatunit) + locuscov.totals$reflen
#XXX These values not currently reported, so no point calculating them
#locuscov.totals$repeatUnitsLwr = locuscov.totals$lwr/nchar(locuscov.totals$repeatunit) + locuscov.totals$reflen
#locuscov.totals$repeatUnitsUpr = locuscov.totals$upr/nchar(locuscov.totals$repeatunit) + locuscov.totals$reflen
# Rename some columns
names(locuscov.totals)[names(locuscov.totals) == 'fit'] <- 'bpInsertion'
names(locuscov.totals)[names(locuscov.totals) == 'coverage'] <- 'decoy_coverage'

# Split locus into 3 columns: chrom start end
locus.cols = matrix(unlist(strsplit(locuscov.totals$locus, '-')), ncol = 3, byrow = T)
colnames(locus.cols) = c('chrom', 'start', 'end')
locuscov.totals = cbind(locuscov.totals, locus.cols)

# Specify output data columns
write.data = locuscov.totals[,c('chrom', 'start', 'end',
                                'sample', 'repeatunit', 'reflen',
                                'locuscoverage', 'decoy_coverage', 'total_assigned',
                                'outlier', 'p_adj', 
                                'bpInsertion', 'repeatUnits'
                                )]
#sort by outlier score then estimated size (bpInsertion), both decending
write.data = write.data[order(write.data$outlier, write.data$bpInsertion, decreasing=T),]
write.data = unique(write.data) #XXX Remove duplicate rows - why do they occur? - no duplicated rows occurring currently

# Write individual files for each sample, remove rows where locuscoverage == 0
samples = unique(write.data$sample)
for (sample in samples) {
  write.table(write.data[write.data$sample == sample & write.data$locuscoverage != 0,], 
            file = paste(c(base.filename, sample, '.', results.suffix), collapse =''), 
            quote = FALSE, row.names=FALSE, sep = '\t')
}

# Write all samples to a single file
write.table(x = write.data, file = paste(c(base.filename, results.suffix), collapse=''), quote = FALSE, row.names=FALSE, sep = '\t')


