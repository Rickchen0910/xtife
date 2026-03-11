## Script to prepare the cigar dataset from the raw CSV file.
## Run once from the package root directory:
##   source("data-raw/cigar.R")
##
## Requires the usethis package.

cigar <- read.csv(
  file.path("..", "regife-main", "data", "cigar.csv"),
  stringsAsFactors = FALSE
)
names(cigar) <- tolower(names(cigar))

# Verify expected dimensions
stopifnot(nrow(cigar) == 1380L, ncol(cigar) >= 9L)

usethis::use_data(cigar, overwrite = TRUE)
