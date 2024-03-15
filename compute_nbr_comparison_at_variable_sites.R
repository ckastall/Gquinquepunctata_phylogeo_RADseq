#!/usr/bin/env Rscript

# Author Chedly Kastally <ckastall@gmail.com>
# Version 1.0
# Copyright (C) 2024 Chedly Kastally <ckastall@gmail.com>
# Modified On 2024-03-13 09:10
# Created  2024-03-13 09:10

# setup ----

dflt_vcf_file      <- "retained_vcf.vcf.gz"
dflt_pop_map_f     <- "population_map.tsv"
dflt_ploidy        <- 2
dflt_output_prefix <- "OutputPolySites"

# Optparse ----

library("optparse")

option_list <-
    list(
         make_option(c("-i", "--vcf_file"),
         type    = "character",
         default = dflt_vcf_file,
         help    = "Input file [Default: %default]"),
         make_option(c("-p", "--population_map"),
         type    = "character",
         default = dflt_pop_map_f,
         help    = "File with sample grouping [Default: %default]"),
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

opt            <- parse_args(parser)
vcf_file       <- opt$vcf_file
output_prefix  <- opt$output_prefix
population_map <- opt$population_map
ploidy         <- opt$ploidy

# read DP table ----

# library("data.table")
vcf_dat_line   <- readLines(vcf_file)
vcf_dat_line <- grep("##", vcf_dat_line, invert = T, value = T)

header <- strsplit(vcf_dat_line[1], "\t")[[1]]
main   <- strsplit(vcf_dat_line[-1], "\t")

# read the population map file ----

pop_dat           <- read.table(population_map, header = F, sep = "\t")
colnames(pop_dat) <- c("sample_id", "population_id")
stopifnot(all(header[-1:-9] %in% pop_dat$sample_id))

sample2pop <- setNames(pop_dat$population_id, pop_dat$sample_id)
unique_pop <- unique(pop_dat$population_id)

# count the number of samples with Depth >= thresholds ----

get_ac0_an <- function(gts) {

    gt_noMiss <- gts[gts != "./."]
    alleles   <- (c(unlist(strsplit(gt_noMiss, "/"))))

    ac0       <- sum(alleles == "0")
    an        <- length(alleles)

    return(c(ac0 = ac0, an = an))
}

list_res <- lapply(main, function(x) {

    res <- x[-1:-9]
    res <- strsplit(res, ":")

    res_gt <- sapply(res, function(y) y[[1]])
    names(res_gt) <- header[-1:-9]

    get_ndiff_ncomp <- function(gts) {
        gt_noMiss <- gts[gts != "./."]
        alleles <- (c(unlist(strsplit(gt_noMiss, "/"))))
        stopifnot(any(alleles != "."))
        ac0 <- sum(alleles == "0")
        an <- length(alleles)
        ndiff <- (ac0 * (an - ac0))
        ncomp <- choose(an, 2)
        return(c(ndiff = ndiff, ncomp = ncomp))
    }

    ## Overall case
    overall_res <- get_ndiff_ncomp(res_gt)
    names(overall_res) <- paste0("overall_", names(overall_res))

    ## per pop case

    res_pop <- lapply(unique_pop, function(c_pop) {
            gts_pop <- res_gt[sample2pop[names(res_gt)] == c_pop]
            res_pop <- get_ndiff_ncomp(gts_pop)
            names(res_pop) <- paste0(c_pop, "_", names(res_pop))
            return(res_pop)
             })

    c(overall_res, unlist(res_pop))

    })

all_res <- do.call(rbind, list_res)
all_sums <- colSums(all_res)

# format and write the result table ----

all_ncomp <- all_sums[(seq_along(all_sums) %% 2) == 0]

## we actually won't use the nDiff: we also need to account for invariant site;
## this is expected to be done in a separate analysis;
## but in case the vcf file had also the invariant site, I leave this in.
all_ndiff <- all_sums[(seq_along(all_sums) %% 2) == 1]

final_table_result <- data.frame(
    population  = gsub("_ncomp", "", names(all_ncomp)),
    n_differences = all_ndiff,
    n_comparisons = all_ncomp
)
rownames(final_table_result) <- NULL

output_table_name <-
    paste0(output_prefix,
           "_OnlyPolymorphicSites_nDiffnComp_perGrouping.tsv")

write.table(final_table_result, file = output_table_name, row.names = F, quote = F, sep = "\t")

if (F) {

    ## Here some code to compute $\pi$ if you also have at hand the data from invariant sites as well 

    ## add the information also from invariant sites
    allSites_counts_file <- "quinqueExclusif_filter_r50p4_AllSites_nComparisons_perGrouping.tsv"
    allSites_counts <- read.table(allSites_counts_file, sep = "\t", header = T)

    allSites_counts$type <- "AllSites"
    final_table_result$type <- "OnlyPolymorphic"

    # compute pi
    final_table_result$pi <- final_table_result$n_differences / allSites_counts$n_comparisons

}
