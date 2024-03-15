#! /bin/bash
#
# Author Chedly Kastally <ckastall@gmail.com>
# Version 1.0
# Copyright (C) 2023 Chedly Kastally <ckastall@gmail.com>
#
# Distributed under terms of the MIT license.
#

depth_table="${1-catalog.call_DPtable.tsv.gz}"
retained_pos_file="${2-retained_positions_noOverlap_Masked_noHighPoly.bed.gz}"
geno2radtag_pos_file="${3-catalog_converted.bed.gz}"

[[ -f "${depth_table}" ]] || { echo "ERROR: the depth table file was not found at path: ${depth_table} ; aborting." ; exit 1 ; } 
[[ -f "${retained_pos_file}" ]] || { echo "ERROR: the retained position file was not found at path: ${retained_pos_file} ; aborting." ; exit 1 ; } 
[[ -f "${geno2radtag_pos_file}" ]] || { echo "ERROR: the index of genomic positions file was not found at path: ${geno2radtag_pos_file} ; aborting." ; exit 1 ; } 

temp_retained_pos_file="${retained_pos_file%.bed.gz}.tmp.bed"
temp_geno2radtag_pos_file="${geno2radtag_pos_file%.bed.gz}.tmp.bed"
retained_position_converted_file="catalog_subset_converted_retained_positions.bed"

## there is no header in those file
zcat "${retained_pos_file}" | awk '{print $1 "_" $2 "_" $3, $0}' > "${temp_retained_pos_file}"
zcat "${geno2radtag_pos_file}" | awk '{print $1 "_" $2 "_" $3, $0}' > "${temp_geno2radtag_pos_file}"

zgrep -w -F -f "${temp_retained_pos_file}" "${temp_geno2radtag_pos_file}" | \
    awk 'BEGIN{FS="\t";OFS="\t"}; \
        {\
            chrom_pos_radtag=$4;
            split(chrom_pos_radtag, a, ":");
            print a[1] ":" a[2] + 1;
        }' >  "${retained_position_converted_file}"

# some clean-up
rm -f "${temp_retained_pos_file}"
rm -f "${temp_geno2radtag_pos_file}"

tmp_depth_table="${depth_table%.tsv.gz}_tmp.tsv"

zcat "${depth_table}" | \
    awk 'BEGIN{FS=OFS="\t"}(/^#/){print};(!/^#/){print $1 ":" $2, $0}' > "${tmp_depth_table}"

final_table_retained_pos_dp_table="catalog_calls_DP_retained.tsv"

zcat "${depth_table}" | head -n 1 > "${final_table_retained_pos_dp_table}"

zgrep -w -F \
    -f "${retained_position_converted_file}" \
    "${tmp_depth_table}" | \
    cut -f2- >> "${final_table_retained_pos_dp_table}"

gzip ${final_table_retained_pos_dp_table}

# clean-up
rm -f "${tmp_depth_table}"
