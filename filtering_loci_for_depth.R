#!/usr/bin/env Rscript

# Author Chedly Kastally <ckastall@gmail.com>
# Version 1.0
# Copyright (C) 2023 Chedly Kastally <ckastall@gmail.com>
# Modified On 2023-02-17 17:13
# Created  2023-02-17 17:13

# Distributed under terms of the MIT license.

## Used to filter the catalog.calls file after conversion using bcftools to
## retain the DP of each sample.

# Optparse ----

# default values

dflt_popfile <- "population_file.tsv"
dflt_infile  <- "example_test_2.tsv.gz"
dflt_pop_threshold <- 0.5
dflt_glo_threshold <- 0.5
dflt_dp_threshold  <- 6

library("optparse")

option_list <-
    list(
         make_option(c("-i", "--input_file"),
         type    = "character",
         default = dflt_infile,
         help    = "Input file [Default: %default]"),
         make_option(c("-p", "--population_file"),
         type    = "character",
         default = dflt_popfile,
         help    = "Population file, tab separated columns with a header including two columns 'sample_id' and 'region'  [Default: %default]"),
         make_option(c("-t", "--max_population_miss_threshold"),
         type    = "numeric",
         default = dflt_pop_threshold,
         help    = "Maximum missing data allowed in a single population, inclusive [Default: %default]"),
         make_option(c("-g", "--max_global_miss_threshold"),
         type    = "numeric",
         default = dflt_glo_threshold,
         help    = "Maximum missing data allowed overall, inclusive --- NOT USED RIGHT NOW [Default: %default]"),
         make_option(c("-d", "--min_depth"),
         type    = "numeric",
         default = dflt_dp_threshold,
         help    = "Minimum depth threshold, inclusive [Default: %default]")
         )

parser <- OptionParser(usage = "%prog -i input_file -p population_file", option_list = option_list)
opt           <- parse_args(parser)

input_file                    <- opt$input_file
population_file               <- opt$population_file
max_global_miss_threshold     <- opt$max_global_miss_threshold
max_population_miss_threshold <- opt$max_population_miss_threshold
min_depth                     <- opt$min_depth

# setup ----

stopifnot(file.exists(input_file))
stopifnot(file.exists(population_file))
stopifnot(is.numeric(min_depth))
stopifnot(is.numeric(max_global_miss_threshold))
stopifnot(is.numeric(max_population_miss_threshold))

# reading the population file ----

pop_dat <- read.table(population_file, sep = "\t", header = T)

pop_expected_colnames <- c("sample_id", "region")

stopifnot(all(pop_expected_colnames %in% colnames(pop_dat)))

sample2group <-
    setNames(
        pop_dat$region,
        pop_dat$sample_id
    )

# reading the input file ----

con <- file(input_file, open = "r")

## Handling the header ----

header <- readLines(con, n = 1)
header <- unlist(strsplit(header, "\t"))
header <- gsub("#?\\s?\\[[0-9]+\\]", "", header)
header <- gsub(":.*", "", header)

expected_header <- c("CHROM", "POS")
stopifnot(all(expected_header %in% header[1:2]))

if (!all(header[-1:-2] %in% pop_dat$sample_id)) {
    cat(sprintf("the following samples from the input file were not found in the population file:\n"))
    cat(sprintf("%s\n\n", paste0(header[-1:-2][which(!(header[-1:-2] %in% pop_dat$sample_id))], collapse = "\n")))
    stop("Aborting, not all samples are in a group.") 
}

unique_groups <- unique(pop_dat$region)

converted_header <- sample2group[header[-2:-1]]

list_index_groups <- sapply(unique_groups, function(x) which(converted_header == x))
list_size_groups  <- sapply(unique_groups, function(x) sum(converted_header == x))

## Process the file line by line ----

res <- vector(mode = "list")
line_count <- 0
i <- 1

# print header
## Only good if I produce a bed file, not the case right now.
# cat("CHROM", "POS0", "POS", sep = "\t")
# cat("\n")

while (TRUE) {

    line <- readLines(con, n = 1)
    if (length(line) == 0) {
        break
    }

    line_count <- line_count + 1
    # cat(sprintf("processing line: %s\n", line_count))

    line_split <- unlist(strsplit(line, "\t"))

    vec <- as.numeric(line_split[-2:-1])

    missing_rate_per_population <-
        sapply(seq_along(list_index_groups), function(x) {
               sum(vec[list_index_groups[[x]]] < min_depth) / list_size_groups[[x]]
        })

    missing_rate_test_population <- missing_rate_per_population <= max_population_miss_threshold

    if (all(missing_rate_test_population)) {
        # cat(line)
	## This would produce a bed file
        # cat(line_split[1], as.numeric(line_split[2]) - 1, line_split[2], sep = "\t")

	## this would produce a list of locusId
        cat(line_split[1], ":", line_split[2], sep = "")
        cat("\n")

    } else {

        # cat(sprintf("SNP excluded at CHROM: %s, POS: %s, population missing thresholds: %s; ", line_split[1], line_split[2], paste0(missing_rate_per_population, collapse = ", ")))
        # cat(sprintf("number of population PASSED the missing rate threshold: %s out of %s\n", sum(missing_rate_test_population), length(missing_rate_test_population)))

    }

}

# cat(sprintf("\n\nDONE\n"))
