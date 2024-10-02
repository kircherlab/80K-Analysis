#!/bin/bash

# input
input_annotation_gtf=/home/kisa/coding/80K_MPRA/gencode_v42/gencode.v42.annotation.gtf.gz
input_annotation_gtf=/home/kisa/coding/80K_MPRA/gencode_v42/gencode.v42.gtf.gz

# output
# only transcripts
only_transcripts=/home/kisa/coding/80K_MPRA/gencode_v42/gencode.v42.transcripts.bed.gz
only_trans_log=/home/kisa/coding/80K_MPRA/gencode_v42/gencode.v42.transcripts.log
# collapsed
collapsed_transcripts=/home/kisa/coding/80K_MPRA/gencode_v42/gencode.v42.transcripts.collapsed.bed.gz
collapsed_trans_log=/home/kisa/coding/80K_MPRA/gencode_v42/gencode.v42.transcripts.collapsed.log


zcat $input_annotation_gtf | \
awk '$3 == "transcript"' | \
awk -v "OFS=\\t" '{{start = $4; if ($7 == "-"){{start = $5}}; print $1,start-1,start,$16"|"$10,".",$7}}' | \
sed 's/[\\";]//g' | \
sort -k1,1 -k2,2n | \
bgzip -c > $only_transcripts 2> $only_trans_log


zcat $only_transcripts | \
sort -k4,4 -k1,1 -k 2,2n | \
awk -v "OFS=\\t" 'BEGIN {{
    id="";chr="";start="";end="";strand=""
}}{{
    if ($4 == id) {{
        end=$3
    }} else {{
        if (NR > 1) {{
            print chr,start,end,id,".",strand;
        }};
        id=$4;chr=$1;start=$2;end=$3;strand=$6;
    }}
}} END {{
    print chr,start,end,id,".",strand
}}' | \
sort -k1,1 -k2,2n | \
bgzip -c > $collapsed_transcripts 2> $collapsed_trans_log