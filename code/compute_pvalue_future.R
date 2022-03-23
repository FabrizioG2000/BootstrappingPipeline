# Empirical p-value process re-factored as functions
library(tidyverse)
library(GenomicRanges)
library(valr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(crayon)
library(pbapply)
library(future)
library(future.apply)

source("code/common.R")


options(scipen = 999999999)
res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
res_num <- c(1e6, 5e5, 1e5, 5e4, 1e4, 5e3)
names(res_num) <- res_set

WORKERS_NUM <- 3

#-----------------------------------------
## Utility functions
feature_annotation_fn <- function(txdb, chr_feature_Grange, fn_file) {
    peakAnno <- annotatePeak(chr_feature_Grange, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = F)

    rn_annotation <- sample(peakAnno@annoStat$Feature, size = length(chr_feature_Grange), prob = peakAnno@annoStat$Frequency / 100, replace = T)
    # check number of peaks from that category
    n5 <- length(grep("5'", as.character(rn_annotation)))
    n3 <- length(grep("3'", as.character(rn_annotation)))
    nexon <- length(grep("Exon", as.character(rn_annotation)))
    nintron <- length(grep("Intron", as.character(rn_annotation)))
    n1kb <- length(grep("1kb", as.character(rn_annotation)))
    n2kb <- length(grep("2kb", as.character(rn_annotation)))
    n3kb <- length(grep("3kb", as.character(rn_annotation)))
    ndown <- length(grep("Down", as.character(rn_annotation)))
    ninter <- length(grep("Inter", as.character(rn_annotation)))
    n_vec <- c(n3, n5, ndown, nexon, ninter, nintron, n1kb, n2kb, n3kb)
    names(n_vec) <- fn_file
    rm(ChIPseekerEnv, envir = globalenv())
    rm(n5, n3, nexon, nintron, n1kb, n2kb, n3kb, ndown, ninter)
    n_vec <- n_vec[n_vec > 0]
    return(n_vec)
}

rn_feature_GRange_build_fn <- function(n_vec, fn_bed_l, hg19_coord, tmp_cage_tbl, fn_file) {

    # The environment from the function
    fn_env <- environment()

    # The empty vector of functions
    rn_fn_coord_l <- vector("list", length(n_vec))
    names(rn_fn_coord_l) <- names(n_vec)

    # looping over the features names
    f <- names(n_vec)[1]
    for (f in names(n_vec)) {

        # How many feature do we need to generate?
        tmp_n <- n_vec[f]
        if (f == fn_file[5]) {
            plan(multisession, workers = WORKERS_NUM)

            # applies from 1 to 100 to
            rn_fn_coord_l[[f]] <- future_lapply(1:100, future.seed = T, future.packages = c("dplyr", "valr"), future.envir = fn_env, function(x) {
                # Try pattern to not abort at the shuffling step
                ## Sampleing prior to shuffling ensures a greater likelihood of success
                # 1) Sample n features from the features table.
                # 2) tell the sampler how big is the chromosome (so that we don't sample outside it)
                # 3) tell the sampler the region in which to sample those features
                # 4) within as we want to shuffle within the chromosome
                rn_pol <- valr::bed_shuffle(tmp_cage_tbl %>% dplyr::sample_n(tmp_n), genome = hg19_coord, excl = fn_bed_l[[f]], within = T, max_tries = 1e6)
                return(rn_pol)
            })

            plan(sequential)
        }

        # If the feature is not chr22_inter_no
        if (f != fn_file[5]) {
            plan(multisession, workers = WORKERS_NUM)

            # Create a list of 100 random situations
            rn_fn_coord_l[[f]] <- future_lapply(1:100, future.seed = T, future.packages = c("dplyr", "valr"), future.envir = fn_env, function(x) {
                # Try pattern to not abort at the shuffling step
                ## Sampleing prior to shuffling ensures a greater likelihood of success
                # 1) Sample n features from the features table.
                # 2) tell the sampler how big is the chromosome (so that we don't sample outside it)
                # 3) tell the sampler the region in which to sample those features
                # 4) within as we want to shuffle within the chromosome
                # out a tibble of size n
                # (this basically moves the intervals from tmp_cage_tbl in places which are contained in the regions defined by fn_bed_l[[f]])
                rn_pol <- try(valr::bed_shuffle(tmp_cage_tbl %>% sample_n(tmp_n), genome = hg19_coord, incl = fn_bed_l[[f]], within = T, max_tries = 1e6), silent = TRUE)
                return(rn_pol)
            })

            plan(sequential)

            # Collect the successful shuffling by eliminating the shuffles producing try-error objects
            rn_fn_coord_l[[f]] <- rn_fn_coord_l[[f]][!(unlist(lapply(rn_fn_coord_l[[f]], function(x) any(class(x) %in% "try-error"))))]
        }

        message(paste0("Finished ", f, " shuffling"))
    }

    warnings()

    # Assemble these blocks into n rnadom results
    bootstrap_count <- 10

    # Create 10 random results by:
    # looping over 10 and joining the results from the random sampling of one of the features.
    # For example, loop over chr22_inter_no, chr22_inter_no_1, chr22_inter_no_2, ..., chr22_inter_no_10
    # And join the results obtained from sampling one of the 100 available chr22_inter_no features
    plan(multisession, workers = WORKERS_NUM)
    rn_peak_coord_tbl_l <- future_lapply(1:bootstrap_count, future.seed = T, future.envir = fn_env, future.packages = c("dplyr"), function(x) {
        return(do.call(bind_rows, lapply(rn_fn_coord_l, function(f) f[[sample(1:length(f), 1)]])))
    })
    plan(sequential)


    # Remove the list of random functional features from this environment
    rm(rn_fn_coord_l, envir = fn_env)

    # convert the list of random functional features into GRanges
    rn_Grange_l <- pblapply(rn_peak_coord_tbl_l, function(x) {
        rnp_Grange <- GRanges(
            seqnames = x$chrom,
            ranges = IRanges(
                start = x$start,
                end = x$end
            )
        )
        return(rnp_Grange)
    })

    # Remove the original list
    rm(rn_peak_coord_tbl_l, envir = fn_env)

    return(rn_Grange_l)
}


## The main function
# chromo is the chromosome name, so for example chr1
# cl_folder is the folder where the clusters are contained
# cl_file is the suffix of the cluster file
# feature_GRange is the GRange object of the features
# fn_repo is "../data/annotation/", the folder where the annotation files are
empirical_pval_compute_fn <- function(chromo, cl_folder, cl_file, feature_Grange, fn_repo, txdb, hg19_coord) {

    # The main environment of the function
    main_fn_env <- environment()

    # the features of this chromosome
    chr_feature_Grange <- feature_Grange[seqnames(feature_Grange) == chromo]

    # Logging the start of the bootstrappinf
    cat(green(chromo), " Bootstrapping started \n")

    # The folder containing the annotations about this chromosome: "../data/annotation/fn_BED/chr1/"
    fn_folder <- paste0(fn_repo, chromo, "/")

    # All the bed files for this chromosome
    fn_file <- grep("BED$", list.files(fn_folder), value = T)

    # Parsing all the bed files
    # The bed files are parsed and put in a dictionary where the key is the file name and the value is the bed file
    fn_bed_l <- lapply(fn_file, function(f) {
        read_bed(paste0(fn_folder, f), n_fields = 3)
    })
    names(fn_bed_l) <- fn_file

    # Copying the txdb and retaining only the rows relative to this chromosome
    txdb_chr <- txdb
    seqlevels(txdb_chr) <- chromo

    # Logging that the annotation has started
    cat(green(chromo), " feature annotation Started \n")

    # Re-importing ChIPseeker
    suppressMessages(library(ChIPseeker))
    message("Imported ChIPseeker")

    # The number of features for each type
    n_vec <- feature_annotation_fn(txdb_chr, chr_feature_Grange, fn_file)

    # Loading the clusters for this chromosome
    cl_chr_tbl <- load_data((paste0(cl_folder, chromo, cl_file)))

    # Removing ChIPseeker package
    detach("package:ChIPseeker", unload = TRUE)
    message("Removed ChIPseeker")

    # Logging that the GRanges have started being built
    cat(green(chromo), " Counting overlaps \n")

    # table collecting the observed CAGE-peak coordinates
    # counting the number of overlaps between the features and the clusters
    # looping over each list of GRanges for each cluster
    plan(multisession, workers = WORKERS_NUM)
    cl_inter_vec <- unlist(future_lapply(cl_chr_tbl$GRange, future.packages = c("GenomicRanges"), function(x) {

        # Counting the overlaps between a given feature and a cluster
        sum(countOverlaps(x, chr_feature_Grange))
    }))
    plan(sequential)

    # Adds a column containing the number of overlaps to the cluster table and filters for the clusters having more than zero overlap
    cl_chr_tbl <- cl_chr_tbl %>%
        mutate(feature_n = cl_inter_vec) %>%
        filter(feature_n > 0)

    # Converting the Granges to lists.
    cl_list <- GRangesList(cl_chr_tbl$GRange)

    # Create this tmp_cage_tbl by:
    # -1) Converting the chr_feature_Grange to a tibbble
    # -2) Selecting only the seqname, start and end columns
    # -3) Renaming the column seqname to chrom
    tmp_cage_tbl <- chr_feature_Grange %>%
        as_tibble() %>%
        dplyr::select(seqnames, start, end) %>%
        dplyr::rename(chrom = seqnames)

    cat(green(chromo), " building random feature coordinates \n")

    # Building the effective random coordinates on which the p-value will be computed
    rn_Grange_l <- rn_feature_GRange_build_fn(n_vec, fn_bed_l, hg19_coord, tmp_cage_tbl, fn_file)

    # Logging the computation of the pvalues.
    cat(green(chromo), " computing empirical p-value \n")

    # counting the overlaps between the clusters and the newly created random features
    # The resulting vector will be a list of 10 (number of bootstrap) lists of overlaps count for each cluster.
    plan(multisession, workers = WORKERS_NUM)
    rn_pval_l <- future_lapply(rn_Grange_l, future.packages = c("GenomicRanges"), function(x) {
        countOverlaps(cl_list, x)
    })
    plan(sequential)

    # Removing the useless granges.
    rm(rn_Grange_l)

    # Stack the pvalues one after the other one column is one tow of the rnpval_l
    rn_count_mat <- do.call(cbind, rn_pval_l)

    # Remove the useless list
    rm(rn_pval_l)

    # unlist converts a list to a vector.
    # Apply over the list of clusters, the function which calculates the p-value
    #
    cl_emp_pval <- unlist(lapply(1:nrow(cl_chr_tbl), function(x) {

        # x is the index
        # sums the number of times the overlaps are greater in the randomly generated one than in the observed ones and dividing
        # everything by the number of bootstrap samples.
        # So if the random have a value which is greater than the observed one in 5/10 cases, this will return 0.5
        (sum(rn_count_mat[x, ] >= cl_chr_tbl$feature_n[x]) + 1) / (ncol(rn_count_mat) + 1)
    }))

    # Add the pvalue as a column named emp.pval in the cluster table.
    return(cl_chr_tbl %>% mutate(emp.pval = cl_emp_pval))
}

#-----------------------------------------
cl_folder <- "./data/cluster/wrapped/"
cl_file <- "_wrapped.Rda"
feature_file <- "./data/feature/feature_wrapped.Rda"
out_file <- "./data/pval_tbl.Rda"

# The features in the shape of a GRange object
feature_Grange <- load_data(feature_file)

hg19_coord <- read_delim("./data/annotation/hg19.genome",
    "\t",
    escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE
)
names(hg19_coord) <- c("chrom", "size")

fn_repo <- "./data/annotation/fn_BED/"

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

chr_set <- list.files(fn_repo)

cl_chr_emp_pval_l <- lapply(chr_set, function(chromo) {
    cl_chr_tbl <- empirical_pval_compute_fn(chromo, cl_folder, cl_file, feature_Grange, fn_repo, txdb, hg19_coord)
    seqlevels(txdb) <- seqlevels0(txdb)

    return(cl_chr_tbl)
})

cl_chr_emp_pval_tbl <- do.call(bind_rows, cl_chr_emp_pval_l)
save(cl_chr_emp_pval_tbl, file = out_file)