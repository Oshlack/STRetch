# Remove duplicate loci in TRF annotation bed file

suppressPackageStartupMessages({
  library('optparse', quietly = TRUE)
})

option_list = list(
  make_option(c("-i","--input"),
              help = "Input annotated bed file in the format: chrom start end motif reflen",
              type = "character"),
  make_option(c("-o","--output"),
              help = "Output bed file (will be tab delimited)",
              type = "character")
)

parser = OptionParser(
  usage = paste("%prog [OPTIONS]",
                "Remove duplicate loci in TRF annotation bed file.\n",
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
infile = arguments$i
outfile = arguments$o

STRloci = read.table(infile, col.names = c('chrom', 'start', 'end', 'motif', 'reflen'), 
                     stringsAsFactors = F)
# Sort by position, then by motif length (decending)
STRloci = STRloci[order(STRloci$chrom, STRloci$start, STRloci$end, -nchar(STRloci$motif)),]
# Remove all dublicated loci (retiaing the first, which will be the one with the longest motif)
cleaned.STRloci = STRloci[!duplicated(STRloci[,1:3]),]
write.table(cleaned.STRloci, file = outfile, 
            quote = F, row.names = F, col.names = F, sep = '\t')