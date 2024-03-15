#! /bin/sh
#
# Author Chedly Kastally <ckastall@gmail.com>
# Version 1.0
# Copyright (C) 2023 Chedly Kastally <ckastall@gmail.com>
#
# Distributed under terms of the MIT license.
#

# Uses the catalog.fa.gz produced by STACKS to create a table with genomics
# positions of each RADtags, based on the RADtag names.

in_file="${1-catalog.fa.gz}"

[ -f "${in_file}" ] || { echo "ERROR: catalog file not found, aborting." ; exit 1 ; }

zcat "${in_file}" | \
    awk 'BEGIN{OFS="\t";FS=" "};
        (/^>/){
            # print ;
            radtag=$1;
            sub("^>","",radtag);
            contig=$2;
            sub("pos=","",contig);
            pos=contig;
            sub(":[0-9]+:.*","",contig);
            sub("[A-Za-z0-9.]+:","",pos);
            sign=pos
            sub("[0-9]+:","", sign)
            sub(":[+-]","", pos);
            getline;
            length_seq=length($0);
            # print radtag, contig, pos, length_seq, sign
            if (sign == "+") {
                for (i=0; i <= length_seq; i++) {
                    print contig, pos + i - 1, pos + i, radtag ":" i ":" sign   
                }
            } else if (sign == "-") {
                for (i=length_seq; i >= 0; i=i-1) {
                    print contig, pos - i - 1, pos - i, radtag ":" i ":" sign ;
                }
            }
        }'
