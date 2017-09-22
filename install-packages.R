## Create the personal library if it doesn't exist. Ignore a warning if the directory already exists.
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = TRUE, recursive = TRUE)
## Install packages
install.packages(c('optparse','plyr','dplyr','tidyr','reshape2'), repos="http://cran.rstudio.com/")
