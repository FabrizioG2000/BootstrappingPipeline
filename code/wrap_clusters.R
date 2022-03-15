# Produce Coordinate tables for BHiCect clusters
library(tidyverse)
library(GenomicRanges)
library(furrr)
install.packages("furrr")
options(scipen = 999999999)

# The resolution datatable
res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6, 5e5, 1e5, 5e4, 1e4, 5e3)
names(res_num) <- res_set

#--------------------------------------------------------------------
raw_cluster_folder <- "./data/cluster/raw/"
out_folder <- "./data/cluster/wrapped/"

raw_suffix <- "_spec_res.Rda"
wrapped_suffix <- "_wrapped.Rda"

# Getting a list of the available chromosomes
chr_set <- unlist(lapply(strsplit(grep("chr", list.files(raw_cluster_folder), value = T), split = "_"), "[", 1))

plan(multisession, workers = 3)

for (chromo in chr_set) {
  # Logging the chromosome name
  message(chromo) 

  # The full name of the output file for this chromosome.
  chr_out_file <- paste0(out_folder, chromo, wrapped_suffix)

  # loading the chromosome file
  load(paste0(raw_cluster_folder, chromo, raw_suffix))
  chr_cl_tbl <- tibble(chr = chromo, cl = names(chr_spec_res$cl_member), bins = chr_spec_res$cl_member)
  chr_cl_tbl <- chr_cl_tbl %>% mutate(res = map_chr(cl, function(x) {
    return(strsplit(x, split = "_")[[1]][1])
  }))

  chr_cl_tbl <- chr_cl_tbl %>% mutate(GRange = future_pmap(list(chr, bins, res), function(chr, bins, res) {
    return(GRanges(
      seqnames = chr,
      ranges = IRanges(
        start = as.numeric(bins),
        end = as.numeric(bins) + res_num[res] - 1
      )
    ))
  }))

  save(chr_cl_tbl, file = chr_out_file)
}