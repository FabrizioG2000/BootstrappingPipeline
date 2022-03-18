# Parsing the provided arguments.
#* This is done at the start of the file as the import of the required libraries below is really slow and 
# we want to give the user the feedback immediately
library("optparse")

# The arguments
option_list <- list(
    make_option(c("-d", "--data"),
        type = "character", default = NULL,
        help = "The input dataset [Required]", metavar = "character"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "The output file name. The path is relative to the directory of the input. [Required]", metavar = "character"
    )
)

# The argument parser
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Chacking that the input file has been passed
if (is.null(opt$data) | is.null(opt$output)) {
    stop("Please specify an input and an output file.")
}

# The input file name
suppressMessages(library(crayon))
message(blue("Started Execution of R Script."))

# If no errors are found in the input parameters, proceed with the script.
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
options(scipen = 999999999)
suppressMessages(library(optparse))
message(green("Successfully imported required libraries"))

# the resolution datatable
res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6, 5e5, 1e5, 5e4, 1e4, 5e3)
names(res_num) <- res_set

#-----------------------------------------
## Utility function
wrap_features_fn <- function(input_file) {

  # Load the data from the CAGE_file as an object called cage_coord_tbl
  cage_coord_tbl <- get(base::load(input_file))
  tmp_obj <- names(mget(base::load(input_file)))
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


# The input file, read from the arguments
input_filepath <- opt$data

# Parsing the output path from the provided argument
output_filepath <- opt$output

# Running the function
wrapped_features <- wrap_features_fn(input_filepath)

# Saving the file to disk.
save(wrapped_features, file = output_filepath)

# Print a success message and the output filepath
message(green("Successfully wrapped features."))
message(blue("Data saved to: ", output_filepath))