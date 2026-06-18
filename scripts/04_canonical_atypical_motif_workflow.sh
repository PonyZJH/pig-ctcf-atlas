#!/usr/bin/env bash

###############################################################################
# 1. Paths and parameters
###############################################################################
# Classification of canonical and atypical CTCF binding sites using JASPAR 2024
# CTCF motifs and FIMO.
#
# Inputs:
#   results/cbs/1.bed
#   results/motif_classification/CTCF.fasta
#   resources/jaspar2024/ctcf.meme
#
# The original notes started from an existing CTCF.fasta file, so FASTA
# extraction is not added here.

project_dir="$(pwd)"

consensus_bed="${project_dir}/results/cbs/1.bed"
motif_dir="${project_dir}/results/motif_classification"
consensus_fasta="${motif_dir}/CTCF.fasta"
jaspar_ctcf_meme="${project_dir}/resources/jaspar2024/ctcf.meme"

mkdir -p "${motif_dir}"
cd "${motif_dir}" || exit 1

###############################################################################
# 2. FIMO scan with JASPAR 2024 CTCF motifs
###############################################################################
fimo --max-stored-scores 800000 \
    --oc fimo_out2024 \
    "${jaspar_ctcf_meme}" \
    "${consensus_fasta}"

cut -f 3- fimo_out2024/fimo.tsv > peaks_with_motif.txt

###############################################################################
# 3. Intersect motif hits with CBS coordinates
###############################################################################
bedtools intersect \
    -a "${consensus_bed}" \
    -b <(sed '1d' peaks_with_motif.txt) \
    -wa -wb \
    > bed_motif.bed

awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' bed_motif.bed > intersected.bed

sort intersected.bed | uniq -c > counts.txt
awk '{count[$1]++} END {for (i in count) print i, count[i]}' counts.txt \
    | sort -n \
    > countsnum.txt

###############################################################################
# 4. Canonical and atypical CBS sets
###############################################################################
sed 's/\t/_/g' intersected.bed \
    | sort | uniq \
    > Canonical.txt

grep -v -F -f Canonical.txt <(sed 's/\t/_/g' "${consensus_bed}") \
    > Atypical.txt

sed 's/_/\t/g' Canonical.txt \
    | sort -k1,1 -k2,2n \
    > Canonical.bed

sed 's/_/\t/g' Atypical.txt \
    | sort -k1,1 -k2,2n \
    > Atypical.bed
