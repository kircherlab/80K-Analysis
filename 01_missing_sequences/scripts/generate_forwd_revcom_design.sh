#!/bin/bash
design_fasta=/fast/groups/ag_kircher/MPRA/IGVF_Y1_design/resources/association_data/design_no_duplicates_sequence_and_header.fa
forw_revc_fasta=/data/gpfs-1/users/kisa11_c/work/coding/80K_analysis/01_missing_sequences/results/fixed_design/design_no_duplicates_sequence_and_header_forw_revc.tsv

paste <(
    cat $design_fasta | awk '{{if ($1 ~ /^>/) {{ gsub(/[\]\[]/,"_"); print ">"substr($1,2)"_forw"}}}}';
    cat $design_fasta | awk '{{if ($1 ~ /^>/) {{ gsub(/[\]\[]/,"_"); print ">"substr($1,2)"_revc"}}}}';
) <(
    cat $design_fasta | awk '{{if ($1 ~ /^[^>]/) {{ seq=substr($1,16,270)}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}';
    cat $design_fasta | awk '{{if ($1 ~ /^[^>]/) {{ seq=substr($1,16,270)}}; if ($1 ~ /^>/ && NR!=1) {{print seq; seq=""}}}} END {{print seq}}' | tr ACGTacgt TGCAtgca | rev;
) > $forw_revc_fasta
