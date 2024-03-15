# Motivation

This repository contains code to count the number of sequenced sites when analyzing RAD-seq data and processed with STACKS (v2.62[1]).

Estimating the amount of genomic positions sequenced, both variable and non-variable, is required for several analyses of populations genetics, such as the computation of nucleotide diversity ($\pi$) or for estimating the site frequency spectrum.

While STACKS (v2.62) can output a VCF file retaining all positions of RAD-tags assembled during the analysis, it discards RAD-tags with no variable position found (as of v2.62). However, these RAD-tag with no variant must be accounted for to estimate accurately the amount of the genome sequenced.

[1]: the code present here was not tested against other versions of STACKS.

# Citation

This code was used in:

_Evolution of the geographic range of a European montane leaf beetle in response to climate changes at the end of the Quaternary_

Chedly Kastally <sup>1,2*</sup>, Maeva Sorel <sup>3*</sup>, Flavien Collart <sup>4</sup>, Patrick Mardulyn <sup>3</sup>

$^1$ Department of Forest Sciences, University of Helsinki, Finland;  
$^2$ Viikki Plant Sciences Center, University of Helsinki, Finland;  
$^3$ Evolutionary Biology and Ecology, Université Libre de Bruxelles (ULB), Belgium;  
$^4$ Department of Ecology and Evolution, University of Lausanne, Switzerland;  
$^*$ these authors have contributed equally to this study.

# Requirements

The following files are required during the analysis:

- *catalog.calls*, one of the output files of STACKS's `ref_map.pl` (or `denovo_map.pl`)
- *catalog.fa*, which was produced with *catalog.calls*
- *population_file.tsv*, a population map file with column 1 = the sample ID and column 2 = the grouping ID
- *regular.vcf.gz*, the VCF output of STACKS's script `population` ran following `ref_map.pl` (or `denovo_map.pl`)

The following tools are also required:

- `R`
- `bash` (and various UNIX tools, e.g. `awk` and `zcat`)
- `tabix` and `bgzip`, both included in `htslib` (http://www.htslib.org/)
- `bcftools` (http://www.htslib.org/)
- `bedtools2` (https://github.com/arq5x/bedtools2)

# Steps of the analyses

## 1. Index RAD-tags and their genomic positions

Use the script `stacks_catalog_fa_parser.sh` to convert the *catalog.fa* into an index of genomic positions for all RAD-tags.

```{.bash}

# the script presumes the catalog.fa was compressed --- if it was not, compress it with `gzip`
bash stacks_catalog_fa_parser.sh catalog.fa.gz | \
    gzip -c - > catalog_converted.bed.gz

```

## 2. Convert catalog.calls into a VCF file

Using `bcftools`, convert the catalog.calls file into a VCF file with proper header.

```{.bash}

# create a header from the vcf produced eg by the STACKS's script `population`
zgrep "^#" regular.vcf.gz > header.txt

# decompress the catalog.calls
zcat catalog.calls > catalog.calls.vcf

# adjust the header
bcftools reheader -h header.txt catalog.calls.vcf | \
    bgzip -c > catalog.calls.vcf.bgzip 

# index the vcf file
tabix -p catalog.calls.vcf.bgzip

```

## 3. Generate the table of sequencing depth

Use `bcftools` to convert the VCF file into a table with the depth of sequencing (at the genotype level) for each sample and locus.

```{.bash}

bcftools query -Hf "%CHROM\t%POS[\t%DP]\n" catalog.call.vcf.bgzip | \
            sed -E -e 's/\./0/g' | \
            gzip -c - > catalog.call_DPtable.tsv.gz

```

## 4. Filter loci based on sequencing depth and missing data

Use R script `filtering_loci_for_depth.R` to filter the depth table based on missing data and depth of sequencing using the same filtering as used for the variable sites

```{.bash}

# In this case, we use a minimum depth of 6 for each position and each
# individual and maximum missing data proportion of 50% in each population.

filtering_loci_for_depth.R \
        -i catalog.call_DPtable.tsv.gz \
        -p population_file.tsv \
        --min_depth 6 \
        --max_population_miss_threshold 0.5 | \
    awk 'BEGIN{FS=OFS=":"}{print $1, $2 - 1}' | \
    gzip -c - > retained_positions.txt


```

## 5. Generate a bed file of the retained position

Convert the filtered position from previous step to a bed file with the genomic coordinates obtained in step 1.

```{.bash}

zgrep -E -w \
    -f retained_positions.txt \
    catalog_converted.bed.gz | \
    sort -k 1,1 -k 2,2n -k3,3n | \
    gzip -c - > retained_positions.bed.gz

```


## 6. Removing overlapping positions

RAD-tags can overlap, notably for those RADtag where no variant were found. These need to be removed from the bed file

```{.bash}

zcat retained_positions.bed.gz | \
    cut -f1-3 | \
    uniq -c | \
    sed -E -e 's/^\s+//' | \
    awk 'BEGIN{OFS="\t"}(/^1/){print $2, $3, $4}' | \
    gzip -c - > retained_positions_noOverlap.bed.gz

```

## (optional) Removing positions of the genome using `bedtools2` and a file with the coordinates of the mask

If a file mask.bed contains known positions of the repetitive fraction of the genome that would be better excluded, they can be excluded with `bedtools`:

```{.bash}

bedtools intersect \
    -v \
    -nonamecheck \
    -a mask.bed \
    -b retained_positions_noOverlap.bed.gz | \
    gzip -c - > retained_positions_noOverlap_Masked.bed.gz

```

## (optional) Removing RAD-tags with a high count of polymorphic positions

Identifying the RAD-tags with too many SNPs:

```{.bash}

## NOTE: this assumes the ID column was filled by STACKS with format: "RADtagID:Position"

zgrep -v "^#" regular.vcf.gz | \
    cut -f 3 | \
    awk 'BEGIN{FS=":"}{print $1}' | \
    sort -n | \
    uniq -c | \
    sed -E -e 's/^\s+//' | \
    awk 'BEGIN{OFS="\t"; print "radtag", "nSNPs"}{print $2, $1}' > count_snp_per_radtag.tsv

```

Print the RAD-tags with 8 variants or more:

```{.bash}

awk '(NR>1){if ($2 > 8) {print $1 ":"}}' count_snp_per_radtag.tsv > RAD-tags_to_remove.txt

```

Get the coordinates of those RAD-tags from the converted catalog.fa:

```{.bash}
grep -F \
    -f RAD-tags_to_remove.txt \
    catalog_converted.bed.gz | \
    awk 'BEGIN{OFS="\t"}{print $1, $$2, $3}' | \
    gzip -c - > genomic_coordinates_RAD-tags_to_remove.bed.gz
```

Removing the genomic regions associated with highly polymorphic RAD-tags from the retained genome bed file:

```{.bash}

bedtools intersect \
    -v \
    -nonamecheck \
    -a retained_positions_noOverlap_Masked.bed.gz \
    -b genomic_coordinates_RAD-tags_to_remove.bed.gz | \
    gzip -c - > retained_positions_noOverlap_Masked_noHighPoly.bed.gz

```

## 7. Counting the number of positions sequenced from the final bed file

```{.bash}

awk ' \
    BEGIN{total=0}; \
    {len=$3 - $2;total=total + len}; \
    END{print total}' retained_positions_noOverlap_Masked_noHighPoly.bed.gz

```

## 8. Estimation of nucleotide diversity ($\pi$)

To estimate the nucleotide diversity, we need to look at both variant and invariant positions. And we also need to know if there is missing data or not (in both variant and invariant positions).

$\pi$ (per bp) is equal to the total number of difference between all pairs of sequences over the total number of pairs compared.

To compute the number of comparisons made, we need first to get the number of all sequenced positions (across all RAD tags) and after filtering.

We can use the depth table produced at step 3, and re-filter for all retained positions in the final set (i.e., after also filtering the variant sites) from the last file of filtered positions. 

```{.bash}

extract_retained_positions_from_dp_table.sh \
    catalog.call_DPtable.tsv.gz \
    retained_positions_noOverlap_Masked_noHighPoly.bed.gz \
    catalog_converted.bed.gz

```

A file *catalog_calls_DP_retained.tsv.gz* should have been produced; from this file, we can count the number of allele at each position, and count the number of possible pairs.

The `--depth_treshold` parameter can be used to set as missing genotypes called below that value, here we set it to 6 (i.e. gt with depth of 6 or above are retained).

```{.bash}

./compute_nbr_comparison_at_AllSites.R \
    --input_file "catalog_calls_DP_retained.tsv.gz" \
    --population_map "population_file.tsv" \
    --depth_threshold 6 \
    --ploidy 2 \
    --output_prefix "AllSitesComparisons"

```

Finally, we use the vcf file, containing only polymorphic sites, to count the number of differences across all variable positions.

Notice that if the vcf file already contained all sequenced positions (variant and invariant) then this should provide both the number of comparison and number of difference observed; but often, the invariant positions are not retained in the vcf file.

```{.bash}

./compute_nbr_comparison_at_variable_sites.R \
    --vcf_file "regular.vcf.gz" \
    --population_map "population_file.tsv" \
    --ploidy 2 \
    --output_prefix "PolymorphicSitesComparisons"

```

Finally, one can divide the number of difference obtained from this last step with the number of comparisons from the previous step to obtain $\pi$.
