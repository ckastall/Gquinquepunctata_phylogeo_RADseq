#!/usr/bin/env Rscript

# Author Chedly Kastally <ckastall@gmail.com>
# Version 1.0
# Copyright (C) 2024 Chedly Kastally <ckastall@gmail.com>
# Modified On 2024-03-13 09:10
# Created  2024-03-13 09:10

# setup ----

dflt_dp_retained_pos_f  <- "catalog_calls_DP_retained.tsv.gz"
dflt_retained_samples_f <- "retained_samples.txt"
dflt_pop_map_f          <- "population_map.txt"
dflt_depth_thresh       <- 6
dflt_ploidy             <- 2
dflt_output_prefix      <- "OutputAllSites"

# Optparse ----

library("optparse")

option_list <-
    list(
         make_option(c("-i", "--input_file"),
         type    = "character",
         default = dflt_dp_retained_pos_f,
         help    = "Input file [Default: %default]"),
         make_option(c("-p", "--population_map"),
         type    = "character",
         default = dflt_pop_map_f,
         help    = "File with sample grouping [Default: %default]"),
         make_option(c("-r", "--retained_samples"),
         type    = "character",
         default = dflt_retained_samples_f,
         help    = "File with list of retained samples [Default: %default]"),
         make_option(c("-t", "--depth_threshold"),
         type    = "numeric",
         default = dflt_depth_thresh,
         help    = "Minimum depth threshold (inclusive) [Default: %default]"),
         make_option(c("-d", "--ploidy"),
         type    = "numeric",
         default = dflt_ploidy,
         help    = "Ploidy [Default: %default]"),
         make_option(c("-o", "--output_prefix"),
         type    = "character",
         default = dflt_output_prefix,
         help    = "Output prefix [Default: %default]")
         )

parser <-
    OptionParser(usage = "%prog -i input_file -o output_prefix",
                 option_list = option_list)

opt                <- parse_args(parser)
input_file         <- opt$input_file
output_prefix      <- opt$output_prefix
retained_samples_f <- opt$retained_samples
population_map     <- opt$population_map
depth_threshold    <- opt$depth_threshold
ploidy             <- opt$ploidy

# read DP table ----

# library("data.table")
dp_retained_pos   <- read.table(input_file, sep = "\t", header = T, comment.char = "")

# fix the header of DP table ----

header <- colnames(dp_retained_pos)
header <- gsub("^X\\.+[0-9]+\\.", "", header)
header <- gsub("\\.DP$", "", header)

colnames(dp_retained_pos) <- header

# optional: subset the samples ----

if (!is.null(retained_samples_f)) {

    ## Check the retained samples
    retained_samples <- readLines(retained_samples_f)
    stopifnot(all(retained_samples %in% header))

    dp_retained_pos <-
        dp_retained_pos[, c(1:2, which(header %in% retained_samples))]

    header <- colnames(dp_retained_pos)
}

# read the population map file ----

pop_dat           <- read.table(population_map, header = F, sep = "\t")
colnames(pop_dat) <- c("sample_id", "population_id")

stopifnot(all(header[-2:-1] %in% pop_dat$sample_id))

sample2pop <- setNames(pop_dat$population_id, pop_dat$sample_id)

unique_pop <- unique(pop_dat$population_id)

# count the number of samples with Depth >= thresholds ----

## The overall case: we use all samples together

result_overall <-
    apply(dp_retained_pos[, -2:-1], 1, function(x) {
        sum(x >= depth_threshold)
    })

## the per population case

sample_id_coln     <- colnames(dp_retained_pos)[-1:-2]
population_id_coln <- sample2pop[sample_id_coln]

list_result_per_pop <-
    lapply(unique_pop, function(c_pop){

        c_df <- dp_retained_pos[, which(sample2pop[colnames(dp_retained_pos)] == c_pop)]

        res <- apply(c_df, 1, function(x) {
            sum(x >= depth_threshold)
        })

        return(res)
    })

all_res <- do.call(cbind, list_result_per_pop)
colnames(all_res) <- unique_pop

# compute the number of possible comparisons ----

## Adjust the count so that we have the number of alleles (ie, adjust for
## ploidy)

all_res_2 <-
    cbind(overall = result_overall * ploidy, all_res * ploidy)

all_combn <- apply(all_res_2, 2, function(x) {
        choose(x, 2)
    })

n_comparisons <-
    as.data.frame(colSums(all_combn))

# format and write the result table ----

final_table <-
    data.frame(population = rownames(n_comparisons),
               n_comparisons = c(n_comparisons[, 1]))

output_table_name <-
    paste0(output_prefix,
           "_AllSites_nComparisons_perGrouping.tsv")

write.table(final_table, file = output_table_name, sep = "\t", row.names = F, quote = F)
