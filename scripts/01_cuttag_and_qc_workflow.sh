#!/usr/bin/env bash

###############################################################################
# 1. Paths and parameters
###############################################################################
# This file is a simple analysis record for the pig CTCF CUT&Tag processing and
# QC steps. Bash sections should be run in bash. Blocks marked "R script" should
# be run in R.
#
# sample_table columns:
#   raw_prefix  clean_prefix  sample_id
#
# raw_prefix:
#   Prefix before _R1.fq.gz and _R2.fq.gz in raw_dir.
# clean_prefix:
#   Prefix before 1_val_1.fq.gz and 2_val_2.fq.gz after Trim Galore.
# sample_id:
#   Readable sample label used in summary files. Use the same value as
#   clean_prefix if no extra sample label is needed.

project_dir="$(pwd)"

sample_table="" # Sample table: raw_prefix, clean_prefix, sample_id; no header.
raw_dir="" # Directory containing paired-end raw FASTQ files.
clean_dir="${project_dir}/results/clean_fastq"
bam_dir="${project_dir}/results/bam"
qc_dir="${project_dir}/results/qc"

reference_index="" # Bowtie2 index prefix for the study reference genome (not a FASTA file).
gtf="" # Gene annotation GTF matching the study reference genome.
spp_script="" # phantompeakqualtools run_spp_nodups.R file.

threads=10
genome_size="2.68e9"
mapq=30

mkdir -p "${clean_dir}" "${bam_dir}" "${qc_dir}"/{fragment_length,pbc,frip,phantompeak,bigwig,tss}

###############################################################################
# 2. Read trimming
###############################################################################
while read -r raw_prefix clean_prefix sample_id; do
    [ -z "${raw_prefix}" ] && continue
    [ "${raw_prefix#\#}" != "${raw_prefix}" ] && continue

    trim_galore -j "${threads}" --paired \
        "${raw_dir}/${raw_prefix}_R1.fq.gz" \
        "${raw_dir}/${raw_prefix}_R2.fq.gz" \
        --gzip -o "${clean_dir}"
done < "${sample_table}"

###############################################################################
# 3. Alignment, filtering and peak calling
###############################################################################
while read -r raw_prefix clean_prefix sample_id; do
    [ -z "${raw_prefix}" ] && continue
    [ "${raw_prefix#\#}" != "${raw_prefix}" ] && continue

    sample="${clean_prefix}"
    work_dir="${bam_dir}/${sample}.work"
    mkdir -p "${work_dir}"

    cp "${clean_dir}/${sample}1_val_1.fq.gz" "${clean_dir}/${sample}2_val_2.fq.gz" "${work_dir}/"
    cd "${work_dir}" || exit 1

    bowtie2 -p "${threads}" -x "${reference_index}" \
        --local --very-sensitive --no-mixed --no-discordant --phred33 \
        -I 10 -X 700 \
        -1 "${sample}1_val_1.fq.gz" \
        -2 "${sample}2_val_2.fq.gz" \
        -S "${sample}_bowtie2.sam" \
        > "${sample}_bowtie2.txt" 2>&1

    samtools view -F 0x04 -@ "${threads}" "${sample}_bowtie2.sam" \
        | awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' \
        | sort | uniq -c \
        | awk -v OFS="\t" '{print $2, $1/2}' \
        > "${sample}_fragmentLen.txt"

    samtools view -bS -F 0x04 -@ "${threads}" \
        "${sample}_bowtie2.sam" > "${sample}_bowtie2.mapped.predup.bam"

    bedtools bamtobed -bedpe -i "${sample}_bowtie2.mapped.predup.bam" \
        | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' \
        | grep -v 'chrM' \
        | sort | uniq -c \
        | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
        > "${sample}.pbc.qc"

    rm "${sample}_bowtie2.mapped.predup.bam"

    picard SortSam \
        I="${sample}_bowtie2.sam" \
        O="${sample}_bowtie2.sorted.sam" \
        SORT_ORDER=coordinate

    rm "${sample}_bowtie2.sam"

    picard MarkDuplicates \
        INPUT="${sample}_bowtie2.sorted.sam" \
        OUTPUT="${sample}_bowtie2.sorted.rmDup.sam" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        METRICS_FILE="${sample}_sorted.metrics" \
        VALIDATION_STRINGENCY=LENIENT

    rm "${sample}_bowtie2.sorted.sam"

    samtools view -h -q "${mapq}" \
        "${sample}_bowtie2.sorted.rmDup.sam" \
        > "${sample}_bowtie2.sorted.rmDup.qualityScore.sam"

    rm "${sample}_bowtie2.sorted.rmDup.sam"

    samtools view -bS -F 0x04 -@ "${threads}" \
        "${sample}_bowtie2.sorted.rmDup.qualityScore.sam" \
        > "${sample}_bowtie2.mapped.bam"

    rm "${sample}_bowtie2.sorted.rmDup.qualityScore.sam"
    samtools index -@ "${threads}" "${sample}_bowtie2.mapped.bam"

    macs2 callpeak -f BAMPE \
        -t "${sample}_bowtie2.mapped.bam" \
        -g "${genome_size}" \
        --qvalue 0.01 \
        --nomodel --shift 0 \
        -n "${sample}" \
        -B --SPMR --keep-dup all

    cp "${sample}_bowtie2.txt" "${sample}.pbc.qc" "${sample}_fragmentLen.txt" \
        "${sample}_sorted.metrics" "${sample}_bowtie2.mapped.bam"* \
        "${sample}"*xls "${sample}"*narrowPeak "${sample}"*summits.bed \
        "${bam_dir}/"

    cp "${sample}_fragmentLen.txt" "${qc_dir}/fragment_length/"
    cp "${sample}.pbc.qc" "${qc_dir}/pbc/"

    cd - >/dev/null || exit 1
    rm -rf "${work_dir}"
done < "${sample_table}"

###############################################################################
# 4. FRiP score
###############################################################################
while read -r raw_prefix clean_prefix sample_id; do
    [ -z "${raw_prefix}" ] && continue
    [ "${raw_prefix#\#}" != "${raw_prefix}" ] && continue

    sample="${clean_prefix}"
    bam="${bam_dir}/${sample}_bowtie2.mapped.bam"
    peak="${bam_dir}/${sample}_peaks.narrowPeak"

    total_reads=$(samtools view -c "${bam}")
    reads_in_peak=$(
        bedtools sort -i "${peak}" \
            | bedtools merge -i stdin \
            | bedtools intersect -u -a "${bam}" -b stdin -ubam \
            | samtools view -c
    )

    {
        echo "Reads_in_peak: ${reads_in_peak}"
        echo "totalReads: ${total_reads}"
        echo "FRiP value:" $(bc <<< "scale=2;100*${reads_in_peak}/${total_reads}")"%"
    } > "${qc_dir}/frip/${sample}_FRIP_score.txt"
done < "${sample_table}"

###############################################################################
# 5. Phantompeakqualtools enrichment metrics
###############################################################################
while read -r raw_prefix clean_prefix sample_id; do
    [ -z "${raw_prefix}" ] && continue
    [ "${raw_prefix#\#}" != "${raw_prefix}" ] && continue

    sample="${clean_prefix}"
    Rscript "${spp_script}" -tmpdir="./" \
        -c="${bam_dir}/${sample}_bowtie2.mapped.bam" \
        -savp -out="${qc_dir}/phantompeak/${sample}.qual" \
        > "${qc_dir}/phantompeak/${sample}.Rout"
done < "${sample_table}"

###############################################################################
# 6. BigWig and TSS signal
###############################################################################
while read -r raw_prefix clean_prefix sample_id; do
    [ -z "${raw_prefix}" ] && continue
    [ "${raw_prefix#\#}" != "${raw_prefix}" ] && continue

    sample="${clean_prefix}"

    bamCoverage -p "${threads}" \
        --normalizeUsing RPKM \
        -b "${bam_dir}/${sample}_bowtie2.mapped.bam" \
        -o "${qc_dir}/bigwig/${sample}.bw"

    computeMatrix scale-regions \
        -p "${threads}" \
        -b 2500 -a 2500 \
        -R "${gtf}" \
        -S "${qc_dir}/bigwig/${sample}.bw" \
        -o "${qc_dir}/tss/${sample}_TSS.gz" \
        --outFileSortedRegions "${qc_dir}/tss/${sample}_genes.bed" \
        --outFileNameMatrix "${qc_dir}/tss/${sample}_matrix.tab"

    plotHeatmap -m "${qc_dir}/tss/${sample}_TSS.gz" \
        -out "${qc_dir}/tss/${sample}_TSS.png"
done < "${sample_table}"

###############################################################################
# 7. QC summary table
###############################################################################
echo -e "sample\tFRiP_percent\tPBC_mt\tPBC_m0\tPBC_m1\tPBC_m2\tNRF\tPBC1\tPBC2" \
    > "${qc_dir}/qc_summary.tsv"

while read -r raw_prefix clean_prefix sample_id; do
    [ -z "${raw_prefix}" ] && continue
    [ "${raw_prefix#\#}" != "${raw_prefix}" ] && continue

    sample="${clean_prefix}"
    frip=$(awk '/FRiP value:/ {gsub("%","",$3); print $3}' "${qc_dir}/frip/${sample}_FRIP_score.txt")
    awk -v sample="${sample}" -v frip="${frip}" \
        'BEGIN{OFS="\t"} {print sample, frip, $1, $2, $3, $4, $5, $6, $7}' \
        "${qc_dir}/pbc/${sample}.pbc.qc" >> "${qc_dir}/qc_summary.tsv"
done < "${sample_table}"

###############################################################################
# R script: average fragment length plot and sample clustering
###############################################################################
# Run this block in R after the Bash sections above.

library(ggplot2)

project_dir <- getwd()
qc_dir <- file.path(project_dir, "results/qc")
fragment_dir <- file.path(qc_dir, "fragment_length")
fragment_plot <- file.path(qc_dir, "total_fragment_len_mean.pdf")
density_matrix <- file.path(project_dir, "results/cbs/TPM_normalized_density.txt")
cluster_tree <- file.path(qc_dir, "sample_signal_wardD_top10pct.newick")

fragment_files <- list.files(fragment_dir, pattern = "fragmentLen.txt$", full.names = TRUE)
result <- NULL
for (path in fragment_files) {
  tmp <- read.table(path, header = FALSE)
  sample <- sub("_fragmentLen.txt$", "", basename(path))
  colnames(tmp)[1:2] <- c("fragment_length", sample)
  tmp <- tmp[, 1:2]
  if (is.null(result)) {
    result <- tmp
  } else {
    result <- merge(result, tmp, by = "fragment_length", all = TRUE)
  }
}

result$mean_count <- rowMeans(result[, -1, drop = FALSE], na.rm = TRUE)

pdf(fragment_plot, width = 5, height = 4)
ggplot(result, aes(x = fragment_length, y = mean_count)) +
  geom_line(size = 1, color = "blue") +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 800)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    legend.position = "none",
    axis.title = element_text(size = 16, color = "black"),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, color = "black", size = 16),
    axis.text.y = element_text(size = 16, color = "black")
  )
dev.off()

data <- read.table(density_matrix, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
data$sd <- apply(data, 1, sd, na.rm = TRUE)
data <- data[order(data$sd, decreasing = TRUE), ]
datt <- data[1:(round(nrow(data) * 0.1)), 1:(ncol(data) - 1)]
distM <- 1 - cor(datt, method = "spearman", use = "pairwise.complete.obs")
hc <- hclust(as.dist(distM), method = "ward.D")
ape::write.tree(ape::as.phylo(as.dendrogram(hc)), file = cluster_tree)

###############################################################################
# End R script
###############################################################################
