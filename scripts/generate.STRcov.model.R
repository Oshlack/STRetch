# Generate the file STRcov.model.csv
# Contains the relationship between the allele size and the normalised coverage
# Specifically the number of bp inserted relative to the reference and 
# the STR-decoy coverage normalised by median coverage calculated with goleft
# Columns:
# sample,repeatunit.x,coverage,coverage_norm,coverage_log,chrom,start,stop,genotype,allele1,allele2,repeatunit.y,locus,allele_mean

#library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)

data.dir = '/Users/hd_vlsci/Documents/git/micro-genotyper/repeat_genotyper_R/data/ataxia_wgs/ATXN8-sim'
bed.files = list.files(data.dir, 'sorted.bed$', full.names = TRUE)
STRcov.files = list.files(data.dir, 'STR_counts$', full.names = TRUE)
locuscov.files = list.files(data.dir, 'locus_counts$', full.names = TRUE)
outfilename = 'scripts/STRcov.model.ATXN8_AGC.csv'

get.sim.no = function(filename) {
  dir.split = tail(strsplit(filename, '/', fixed = TRUE)[[1]], n = 1)
  strsplit(dir.split, '(\\.|_)', fixed = FALSE)[[1]][1]
}

parse.gt = function(gt) {
  split1 = strsplit(gt, '_', fixed = TRUE)[[1]][2]
  repunit = strsplit(gt, '_', fixed = TRUE)[[1]][1]
  gt = strsplit(split1, '/', fixed = TRUE)[[1]]
  gt = sort(as.numeric(gt))
  return(c(gt, repunit))
}

parse.bed.allele = function(filename) {
  sim.id = get.sim.no(filename)
  bed.data = read.table(filename, sep = '\t', stringsAsFactors = FALSE)
  gt = t(sapply(bed.data[ ,4], parse.gt))
  id.gt = cbind(sim.id, gt)
  rownames(id.gt) = NULL
  bed.data = cbind(bed.data, id.gt)
  colnames(bed.data) = c("chrom", "start", "stop", "genotype", "sample", "allele1", "allele2", "repeatunit")
  bed.data$locus = paste(bed.data$chrom, bed.data$start, bed.data$stop, sep='-')
  return(bed.data)
}

# Parse all STR coverage
parse.STRcov = function(filename) {
  sim.id = get.sim.no(filename)
  cov.data = read.table(filename, stringsAsFactors = FALSE, col.names = c('chrom', 'start', 'end', 'coverage'))
  cov.data$sample = sim.id
  cov.data$repeatunit = sapply(cov.data$chrom, function(x){ strsplit(x, '-', fixed = TRUE)[[1]][2] } )
  return(cov.data[,c('sample', 'repeatunit', 'coverage')])
}

# Parse all locus coverage
parse.locuscov = function(filename) {
  sim.id = get.sim.no(filename)
  locuscov.data = read.table(filename, stringsAsFactors = FALSE, header=TRUE)
  locuscov.data$sample = sim.id
  locuscov.data$locus = paste(locuscov.data$STR_chr, locuscov.data$STR_start, locuscov.data$STR_stop, sep='-')
  locuscov.data$repeatunit = locuscov.data$motif
  locuscov.data$locuscoverage = locuscov.data$count
  return(locuscov.data[,c('sample', 'locus', 'repeatunit', 'locuscoverage')])
}

# Main

bed.data = ldply(lapply(bed.files, parse.bed.allele))
# as character required to avoid values becoming factor levels
bed.data$allele1 = as.numeric(as.character(bed.data$allele1)) 
bed.data$allele2 = as.numeric(as.character(bed.data$allele2))
bed.data$allele_mean = (bed.data$allele1 + bed.data$allele2)/2

## Parse coverage
# parse.STRcov
STRcov.data = ldply(lapply(STRcov.files, parse.STRcov), data.frame)
#STRcov.data$coverage = as.numeric(STRcov.data$coverage)

# parse.locuscov
locuscov.data = ldply(lapply(locuscov.files, parse.locuscov), data.frame)
#locuscov.data$locuscoverage = as.numeric(locuscov.data$locuscoverage)
# Fill zeros in locuscov
locuscov.data.wide = spread(locuscov.data, sample, locuscoverage)
locuscov.data.wide[is.na(locuscov.data.wide)] = 0
sample.cols = which(names(locuscov.data.wide) %in% unique(locuscov.data$sample))
locuscov.data = gather(locuscov.data.wide, sample, locuscoverage, sample.cols)

# Normalise counts 
median_cov = 30 # Should be obtained from data
factor = 100
STRcov.data$coverage_norm = factor * (STRcov.data$coverage + 1) / median_cov
STRcov.data$coverage_log = log2(STRcov.data$coverage_norm)
locuscov.data$locuscoverage_norm = factor * (locuscov.data$locuscoverage + 1) / median_cov
locuscov.data$locuscoverage_log = log2(locuscov.data$locuscoverage_norm)

# Extract the relevant info for ATXN8
my.RU = 'AGC'
my.locus = 'chr13-70713515-70713561'
ref.allele = 15.3

STRcov.data.AGC = STRcov.data[STRcov.data$repeatunit == my.RU,]
STRcov.bed.AGC = merge(STRcov.data.AGC, bed.data, by='sample')

locuscov.data.AGC = locuscov.data[locuscov.data$locus == my.locus,]
locuscov.bed.AGC = merge(locuscov.data.AGC, bed.data, by='sample')

# Normalise allele length to bp insertion relative to the reference
STRcov.bed.AGC$allele1 = (locuscov.bed.AGC$allele1 - ref.allele) * nchar(my.RU)
STRcov.bed.AGC$allele2 = (locuscov.bed.AGC$allele2 - ref.allele) * nchar(my.RU)

# Write to file

write.csv(x = STRcov.bed.AGC, file = outfilename, quote = FALSE, row.names=FALSE)