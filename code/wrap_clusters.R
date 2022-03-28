# Extracts the chromosome name from the file name using a wildcard system.
# For example, given the filename chr1_spec_res.Rda and the pattern %_spec_res.Rda, the function will return chr1
extract_from_wildcard <- function(text, pattern, wildcard = "%") {
  out <- text

  # the position of the wildcard in the text
  i <- gregexpr(pattern = wildcard, text = pattern)[[1]]

  # The string which in the pattern comes before the wildcard
  prefix <- substr(pattern, 1, i - 1)


  if (prefix != "") {
    out <- sub(prefix, "", text)
  }

  # The string which in the pattern comes after the wildcard
  suffix <- substr(pattern, i + 1, nchar(pattern))

  if (suffix != "") {
    # remove the last n character from out where n is the length of the suffix
    out <- sub(suffix, "", out, fixed = TRUE)
  }
  return(out)
}

# Parsing the provided arguments.
#* This is done at the start of the file as the import of the required libraries below is really slow and
# we want to give the user the feedback immediately
library("optparse")

# The arguments
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "The folder containing the chromosome files.", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = NULL,
    help = "The output folder", metavar = "character"
  ),
  make_option(c("-p", "--pattern"),
    type = "character", default = NULL,
    help = "The naming pattern.", metavar = "character"
  ),
  make_option(c("-s", "--save"),
    type = "character", default = NULL,
    help = "The saving pattern.", metavar = "character"
  )
)

# The argument parser
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Chacking that the input file has been passed
if (is.null(opt$input) | is.null(opt$out) | is.null(opt$pattern) | is.null(opt$save)) {
  message(opt$input)
  message(opt$output)
  message(opt$pattern)
  message(opt$save)

  stop("PLease insert the correct parameters.")
}

# The input file name
suppressMessages(library(crayon))
message(blue("Started Execution of R Script."))

# Produce Coordinate tables for BHiCect clusters
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(furrr))
suppressMessages(library(crayon))
options(scipen = 999999999)

# The resolution datatable
res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6, 5e5, 1e5, 5e4, 1e4, 5e3)
names(res_num) <- res_set

# I-O
raw_cluster_folder <- "./data/cluster/raw/" # The input folder
out_folder <- "./data/cluster/wrapped2/" # The output folder
input_pattern <- opt$pattern # The input pattern

# input_pattern <- "%_spec_res.Rda"
save_pattern <- opt$save # The saving pattern
# save_pattern <- "%_wrapped.Rda"

# If the output folder does not exist, create it
if (!file.exists(out_folder)) {
  dir.create(out_folder, recursive = TRUE)
}

# Getting a list of the available chromosomes, extracting from the wildcard
chr_set <- (lapply(list.files(raw_cluster_folder), function(f) {
  return(c(extract_from_wildcard(f, input_pattern), f))
}))

if (length(chr_set) > 0) {
  # Logging the number of chromosomes found
  message(blue(paste0("Found ", length(chr_set), " chromosomes.")))
} else {

  # If no chromosomes were found tell the user and exit
  message(red("No chromosomes found matching the provided pattern."))
  stop()
}

# Planning the execution on 3 worker sessions
plan(multisession, workers = 3)

# Looping over the chromosomes
for (chromo_and_filepath in chr_set) {
  chromo <- chromo_and_filepath[1] # the chromosome name
  filepath <- chromo_and_filepath[2] # the original filepath

  # The full name of the output file for this chromosome.
  chr_out_file <- paste0(out_folder, gsub("%", chromo, save_pattern))

  # loading the chromosome file (loads an object called chr_spec_res)
  load(paste0(raw_cluster_folder, "/", filepath))

  # Converting the chromosome file to a tibble
  chr_cl_tbl <- tibble(chr = chromo, cl = names(chr_spec_res$cl_member), bins = chr_spec_res$cl_member)

  # Add a table representing the resolution at which the cluster was detected to the table
  chr_cl_tbl <- chr_cl_tbl %>% mutate(res = map_chr(cl, function(x) {
    return(strsplit(x, split = "_")[[1]][1])
  }))

  # Adding a GRange object to the table
  chr_cl_tbl <- chr_cl_tbl %>% mutate(GRange = future_pmap(list(chr, bins, res), function(chr, bins, res) {
   return(GRanges(
     seqnames = chr,
     ranges = IRanges(
       start = as.numeric(bins),
       end = as.numeric(bins) + res_num[res] - 1
     )
   ))
  }))

  # Saving the file to disk
  save(chr_cl_tbl, file = chr_out_file)

  # Log the chromosome name and the saving path to console
  message((paste0("Saved ", green(chromo), " to ", chr_out_file)))
}