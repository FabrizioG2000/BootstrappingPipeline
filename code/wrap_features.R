library(tidyverse)
library(GenomicRanges)
options(scipen = 999999999)

# the resolution datatable
res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6, 5e5, 1e5, 5e4, 1e4, 5e3)
names(res_num) <- res_set

#-----------------------------------------
## Utility function
build_GRange_fn <- function(CAGE_file) {

  # Load the data from the CAGE_file as an object called cage_coord_tbl
  cage_coord_tbl <- get(base::load(CAGE_file))
  tmp_obj <- names(mget(base::load(CAGE_file)))
  rm(list = tmp_obj)
  rm(tmp_obj)

  # Filtering the dataframe removing the rows with missing values in the field "start"
  cage_coord_tbl <- cage_coord_tbl %>% filter(!(is.na(start)))

  # Converting the data to genomic ranges
  cage_chr_Grange <- GRanges(
    seqnames = cage_coord_tbl$chr,
    ranges = IRanges(
      start = as.numeric(cage_coord_tbl$start),
      end = as.numeric(cage_coord_tbl$end)
    )
  )
  return(cage_chr_Grange)
}

#-----------------------------------------

# The input file
CAGE_feature_in_file <- "./data/feature/feature_in.Rda"

# The output file
out_file <- "./data/feature/feature_wrapped.Rda"

# The CAGE file converted to a genomic range object.
CAGE_feature_wrapped <- build_GRange_fn(CAGE_feature_in_file)

# Saving the file to disk.
save(CAGE_feature_wrapped, file = out_file)
