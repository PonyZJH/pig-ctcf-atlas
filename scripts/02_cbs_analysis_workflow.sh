#!/usr/bin/env bash

###############################################################################
# 1. Paths and parameters
###############################################################################
# This file is a simple analysis record for the CTCF binding site (CBS) analyses.
# Bash sections should be run in bash. Blocks marked "R script" should be run in R.
#
# sample_table columns:
#   raw_prefix  clean_prefix  sample_id
#
# tissue_table columns:
#   SampleID  Original  Filename  Organ

project_dir="$(pwd)"

sample_table="" # Sample table: raw_prefix, clean_prefix, sample_id; no header.
tissue_table="" # Metadata: SampleID, Original, Filename, Organ; no header.

bam_dir="${project_dir}/results/bam"
peak_dir="${project_dir}/results/bam"
bw_dir="${project_dir}/results/qc/bigwig"
cbs_dir="${project_dir}/results/cbs"

genome_fasta="" # Study reference genome FASTA.
gtf="" # Gene annotation GTF matching the reference genome.
chrom_sizes="" # Chromosome sizes: chromosome name and length, tab-separated.
cpg_fasta_dir="" # Directory of chromosome FASTA files (*.fa) for cpgplot.
repeat_dir="" # Directory of repeat annotation BED files matching *rmChrM*bed.

threads=10
genome_size=2682231099

bed_dir="${cbs_dir}/bed"
overlap_dir="${cbs_dir}/bed/overlap"
tissue_bed_dir="${cbs_dir}/bed/tissue"
signal_dir="${cbs_dir}/signal"
annotation_dir="${cbs_dir}/annotation"
motif_dir="${cbs_dir}/motif"
cpg_dir="${cbs_dir}/cpg"

mkdir -p "${bed_dir}" "${overlap_dir}" "${tissue_bed_dir}" \
    "${signal_dir}" "${annotation_dir}" "${motif_dir}" "${cpg_dir}"

consensus_bed="${cbs_dir}/1.bed"

###############################################################################
# 2. Average CTCF signal and TSS heatmap
###############################################################################
input_bw=$(ls "${bw_dir}"/*.bw)

wiggletools mean ${input_bw} \
    | wiggletools write_bg "${signal_dir}/total_repMean.bg" -

sort -k1,1 -k2,2n "${signal_dir}/total_repMean.bg" \
    > "${signal_dir}/total_repMean.sorted.bg"

bedGraphToBigWig \
    "${signal_dir}/total_repMean.sorted.bg" \
    "${chrom_sizes}" \
    "${signal_dir}/total_CTCF.bw"

computeMatrix scale-regions \
    -p "${threads}" \
    -b 2500 -a 2500 \
    -R "${gtf}" \
    -S "${signal_dir}/total_CTCF.bw" \
    -o "${signal_dir}/total_CTCF_TSS.gz"

plotHeatmap \
    -m "${signal_dir}/total_CTCF_TSS.gz" \
    -out "${signal_dir}/total_CTCF_TSS.pdf" \
    --colorMap "RdBu_r"

###############################################################################
# 3. High-confidence single-sample peak filtering
###############################################################################
while read -r raw_prefix clean_prefix sample_id; do
    [ -z "${raw_prefix}" ] && continue
    [ "${raw_prefix#\#}" != "${raw_prefix}" ] && continue

    awk '$8 >= 5 {print $1, $2, $3}' "${peak_dir}/${clean_prefix}"*narrowPeak \
        | grep -Ev 'chrM|NW|AEMK|FPKY' \
        > "${bed_dir}/${sample_id}.bed"
done < "${sample_table}"

###############################################################################
# 4. Consensus CBS merging and statistics
###############################################################################
cat "${bed_dir}"/*.bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${consensus_bed}"

awk 'BEGIN{sum=0; min=""; max=0; n=0}
     {len=$3-$2+1; sum+=len; n++; if(min=="" || len<min) min=len; if(len>max) max=len}
     END{printf "CBS_count\t%d\nMean_length_bp\t%.3f\nMin_length_bp\t%d\nMax_length_bp\t%d\n", n, sum/n, min, max}' \
    "${consensus_bed}" > "${cbs_dir}/consensus_cbs_stats.txt"

awk -v genome_size="${genome_size}" \
    'BEGIN{sum=0} {sum += $3-$2+1} END{printf "Covered_bp\t%d\nGenome_size_bp\t%d\nGenome_percent\t%.6f\n", sum, genome_size, 100*sum/genome_size}' \
    "${consensus_bed}" > "${cbs_dir}/consensus_cbs_genome_coverage.txt"

###############################################################################
# 5. Read-count matrix with multiBamSummary
###############################################################################
bam_list=""
while read -r raw_prefix clean_prefix sample_id; do
    [ -z "${raw_prefix}" ] && continue
    [ "${raw_prefix#\#}" != "${raw_prefix}" ] && continue

    ln -sf "${bam_dir}/${clean_prefix}_bowtie2.mapped.bam" \
        "${cbs_dir}/${sample_id}.bam"
    ln -sf "${bam_dir}/${clean_prefix}_bowtie2.mapped.bam.bai" \
        "${cbs_dir}/${sample_id}.bam.bai"
    bam_list="${bam_list} ${cbs_dir}/${sample_id}.bam"
done < "${sample_table}"

multiBamSummary BED-file \
    --BED "${consensus_bed}" \
    -b ${bam_list} \
    -out "${signal_dir}/peak_bam_density.npz" \
    --outRawCounts "${signal_dir}/peak_bam_density.txt"

sed -i.bak "s/'//g; s/.bam//g" "${signal_dir}/peak_bam_density.txt"

###############################################################################
# 6. Presence/absence overlap files
###############################################################################
for bed in "${bed_dir}"/*.bed; do
    sample=$(basename "${bed}" .bed)
    bedtools intersect -a "${consensus_bed}" -b "${bed}" -wa \
        | sort | uniq \
        > "${overlap_dir}/${sample}_overlap.bed"
done

###############################################################################
# R script: signal normalization, occupancy matrix and CBS plots
###############################################################################
# Run this block in R after sections 1-6 above.

library(ggplot2)

project_dir <- getwd()
cbs_dir <- file.path(project_dir, "results/cbs")
signal_dir <- file.path(cbs_dir, "signal")
bed_dir <- file.path(cbs_dir, "bed")
overlap_dir <- file.path(bed_dir, "overlap")
tissue_bed_dir <- file.path(bed_dir, "tissue")
tissue_table <- "" # Metadata: SampleID, Original, Filename, Organ; no header.
genome_size <- 2682231099

raw_counts <- file.path(signal_dir, "peak_bam_density.txt")
consensus_bed <- file.path(cbs_dir, "1.bed")
tpm_matrix <- file.path(cbs_dir, "TPM_normalized_density.txt")
occupancy_matrix <- file.path(cbs_dir, "CTCF_149sample_01.txt")
tissue_matrix <- file.path(cbs_dir, "merged_tissues.txt")

count_table <- read.table(raw_counts, header = TRUE, check.names = FALSE, sep = "\t")
colnames(count_table) <- sub("^#", "", colnames(count_table))
peak_id <- paste(count_table[[1]], count_table[[2]], count_table[[3]], sep = "_")
peak_length_kb <- (count_table[[3]] - count_table[[2]] + 1) / 1000

N <- as.data.frame(count_table[, 4:ncol(count_table), drop = FALSE], check.names = FALSE)
for (i in 1:ncol(N)) {
  N[, i] <- as.numeric(N[, i]) / peak_length_kb
  N[, i] <- (N[, i] / sum(N[, i], na.rm = TRUE)) * 1e6
}
rownames(N) <- peak_id
write.table(N, file = tpm_matrix, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

consensus <- read.table(consensus_bed, header = FALSE)
consensus_id <- paste(consensus$V1, consensus$V2, consensus$V3, sep = "_")
occupancy <- data.frame(row.names = consensus_id)

overlap_files <- list.files(overlap_dir, pattern = "_overlap.bed$", full.names = TRUE)
for (path in overlap_files) {
  sample <- sub("_overlap.bed$", "", basename(path))
  if (file.info(path)$size > 0) {
    dat <- read.table(path, header = FALSE)
    present <- paste(dat$V1, dat$V2, dat$V3, sep = "_")
  } else {
    present <- character(0)
  }
  occupancy[[sample]] <- as.integer(consensus_id %in% present)
}
write.table(occupancy, file = occupancy_matrix, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

metadata <- read.table(tissue_table, header = FALSE, stringsAsFactors = FALSE)
colnames(metadata)[1:4] <- c("SampleID", "Original", "Filename", "Organ")

dir.create(tissue_bed_dir, recursive = TRUE, showWarnings = FALSE)
tissue_data <- data.frame(row.names = rownames(occupancy))
for (organ in unique(metadata$Organ)) {
  tissue_samples <- intersect(metadata$Filename[metadata$Organ == organ], colnames(occupancy))
  if (length(tissue_samples) == 0) next
  tissue_data[[organ]] <- rowSums(occupancy[, tissue_samples, drop = FALSE], na.rm = TRUE)
  present <- rownames(tissue_data)[tissue_data[[organ]] > 0]
  write.table(present, file = file.path(tissue_bed_dir, paste0(organ, "_peak.txt")),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}
write.table(tissue_data, file = tissue_matrix, sep = "\t", quote = FALSE, row.names = TRUE)

tissue_summary <- data.frame()
for (path in list.files(tissue_bed_dir, pattern = "_peak.txt$", full.names = TRUE)) {
  tissue <- sub("_peak.txt$", "", basename(path))
  ids <- readLines(path)
  coords <- do.call(rbind, strsplit(ids, "_"))
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  coords$V2 <- as.numeric(coords$V2)
  coords$V3 <- as.numeric(coords$V3)
  total_len <- sum(coords$V3 - coords$V2 + 1)
  tissue_summary <- rbind(tissue_summary, data.frame(
    Tissue = tissue,
    Number = length(ids),
    Length = total_len,
    Percentage = total_len / genome_size * 100
  ))
}

pdf(file.path(cbs_dir, "tissue_peak_number_ctcf.pdf"), width = 14, height = 6)
ggplot(tissue_summary, aes(Tissue, Number, fill = "Peak_number")) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("#FDB462")) +
  labs(title = NULL, x = NULL, y = "Number of peaks") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        legend.position = c(0.75, 0.9), legend.title = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black", size = 12),
        axis.text.y = element_text(vjust = 0.5, colour = "black", size = 12))
dev.off()

pdf(file.path(cbs_dir, "tissue_peak_len_ctcf.pdf"), width = 14, height = 6)
ggplot(tissue_summary, aes(Tissue, Percentage, fill = "Peak_length")) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("#B3DE69")) +
  labs(title = NULL, x = NULL, y = "Genome (%)") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        legend.position = c(0.75, 0.9), legend.title = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black", size = 12),
        axis.text.y = element_text(vjust = 0.5, colour = "black", size = 12))
dev.off()

consensus$len <- (consensus$V3 - consensus$V2 + 1) / 1000
mean_len <- mean(consensus$len)
consensus$Group <- "Total"

pdf(file.path(cbs_dir, "CTCF_peak_len_distribution.pdf"), width = 4, height = 4)
ggplot(consensus, aes(x = len, fill = Group)) +
  geom_density(size = 1, alpha = 0.4, color = "#FDB462") +
  geom_vline(xintercept = mean_len, linetype = "dashed", color = "#FDB462") +
  scale_fill_manual(values = "#FDB462") +
  annotate("text", x = 5, y = 0.75,
           label = paste("Mean length =", round(mean_len * 1000), "bp"),
           color = "#FDB462") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 14), breaks = seq(0, 14, 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Length of CTCF peaks (kb)") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        legend.position = "none", axis.line = element_line(color = "black"))
dev.off()

data <- read.table(tpm_matrix, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
data$sd <- apply(data, 1, sd, na.rm = TRUE)
data <- data[order(data$sd, decreasing = TRUE), ]
datt <- data[1:(round(nrow(data) * 0.1)), 1:(ncol(data) - 1)]
distM <- 1 - cor(datt, method = "spearman", use = "pairwise.complete.obs")
hc <- hclust(as.dist(distM), method = "ward.D")
ape::write.tree(ape::as.phylo(as.dendrogram(hc)),
                file = file.path(cbs_dir, "Total_CTCF_sample_TPM_ward.D_0.1.newick"))

###############################################################################
# End R script
###############################################################################

###############################################################################
# 7. HOMER annotation and motif analysis
###############################################################################
annotatePeaks.pl "${consensus_bed}" "${genome_fasta}" -gtf "${gtf}" \
    > "${annotation_dir}/1Ref_all_peaks_merge_nonchr_homer_annotation.txt"

sed 1d "${annotation_dir}/1Ref_all_peaks_merge_nonchr_homer_annotation.txt" \
    | awk '{print $8}' \
    | sort | uniq -c \
    > "${annotation_dir}/1Ref_all_peaks_merge_nonchr_homer_summary.txt"

while read -r raw_prefix clean_prefix sample_id; do
    [ -z "${raw_prefix}" ] && continue
    [ "${raw_prefix#\#}" != "${raw_prefix}" ] && continue

    findMotifsGenome.pl \
        "${bam_dir}/${clean_prefix}_bowtie2.mapped.bam" \
        "${genome_fasta}" \
        "${motif_dir}/${sample_id}_motif" \
        -preparsedDir "${motif_dir}" \
        -mset vertebrates
done < "${sample_table}"

###############################################################################
# 8. CpG island, shore and shelf processing
###############################################################################
for fasta in "${cpg_fasta_dir}"/*.fa; do
    chr_name=$(basename "${fasta}" .fa)
    cpgplot \
        -sequence "${fasta}" \
        -window 100 \
        -minlen 200 \
        -minoe 0.6 \
        -minpc 50. \
        -outfile "${cpg_dir}/${chr_name}.cpgplot" \
        -graph none \
        -outfeat "${cpg_dir}/${chr_name}.gff"
done

cat "${cpg_dir}"/*.gff \
    | grep "cpgplot" \
    | grep -v "chrMT" \
    | awk 'BEGIN{OFS="\t"}{print $1,$4,$5}' \
    | sort -k1,1 -k2,2n \
    > "${cpg_dir}/total_cpg.bed"

bedtools flank -i "${cpg_dir}/total_cpg.bed" -g "${chrom_sizes}" -b 2000 \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${cpg_dir}/total_cpgShore.tmp.bed"

awk 'BEGIN{OFS="\t"}{start=$2-2001; if(start<0) start=0; print $1,start,$3+2001}' \
    "${cpg_dir}/total_cpg.bed" \
    | bedtools flank -i - -g "${chrom_sizes}" -b 2000 \
    | sort -k1,1 -k2,2n \
    | awk '$3>$2' \
    | bedtools merge -i - \
    > "${cpg_dir}/total_cpgShelf.tmp.bed"

awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3+1}' \
    "${cpg_dir}/total_cpg.bed" > "${cpg_dir}/cpg_expanded.bed"

bedtools subtract \
    -a "${cpg_dir}/total_cpgShore.tmp.bed" \
    -b "${cpg_dir}/cpg_expanded.bed" \
    > "${cpg_dir}/total_cpgShore.bed"

awk 'BEGIN{OFS="\t"}{start=$2-1; if(start<0) start=0; print $1,start,$3+1}' \
    "${cpg_dir}/total_cpgShore.bed" > "${cpg_dir}/shore_expanded.bed"

bedtools subtract \
    -a "${cpg_dir}/total_cpgShelf.tmp.bed" \
    -b "${cpg_dir}/cpg_expanded.bed" \
    | bedtools subtract -a - -b "${cpg_dir}/shore_expanded.bed" \
    > "${cpg_dir}/total_cpgShelf.bed"

rm "${cpg_dir}/total_cpgShore.tmp.bed" "${cpg_dir}/total_cpgShelf.tmp.bed" \
    "${cpg_dir}/cpg_expanded.bed" "${cpg_dir}/shore_expanded.bed"

###############################################################################
# 9. CpG and repeat overlap
###############################################################################
for bed in "${cpg_dir}/total_cpg.bed" "${cpg_dir}/total_cpgShore.bed" "${cpg_dir}/total_cpgShelf.bed"; do
    label=$(basename "${bed}")
    bedtools intersect -a "${consensus_bed}" -b "${bed}" -wa \
        | sort | uniq \
        > "${annotation_dir}/CTCF_enhancer_overlap_${label}"
done

for bed in "${repeat_dir}"/*rmChrM*bed; do
    label=${bed##*_}
    bedtools intersect -a "${consensus_bed}" -b "${bed}" -wa \
        | sort | uniq \
        > "${annotation_dir}/CTCF_overlap_${label}"
done

wc -l "${annotation_dir}"/*overlap*bed \
    | grep "bed" \
    > "${annotation_dir}/CTCF_peak_other_annotate.txt"

###############################################################################
# R script: annotation summary plot
###############################################################################
# Run this block in R after HOMER, CpG and repeat overlap sections.

library(ggplot2)

project_dir <- getwd()
cbs_dir <- file.path(project_dir, "results/cbs")
annotation_dir <- file.path(cbs_dir, "annotation")
consensus_bed <- file.path(cbs_dir, "1.bed")
denominator <- nrow(read.table(consensus_bed, header = FALSE))

enh <- read.table(file.path(annotation_dir, "1Ref_all_peaks_merge_nonchr_homer_summary.txt"),
                  header = FALSE, stringsAsFactors = FALSE)
enh$V1 <- enh$V1 / denominator * 100
colnames(enh) <- c("Percent", "Label")
enh$Group <- "genome_annotation"

dat2 <- read.table(file.path(annotation_dir, "CTCF_peak_other_annotate.txt"),
                   header = FALSE, stringsAsFactors = FALSE)
dat2$V2 <- basename(dat2$V2)
dat2$V2 <- gsub("CTCF_overlap_total_", "", dat2$V2)
dat2$V2 <- gsub("CTCF_overlap_", "", dat2$V2)
dat2$V2 <- gsub("CTCF_enhancer_overlap_", "", dat2$V2)
dat2$V2 <- gsub(".bed", "", dat2$V2)
dat2$V2 <- gsub("total_cpg$", "CpG_Island", dat2$V2)
dat2$V2 <- gsub("total_cpgShore", "CpG_Shore", dat2$V2)
dat2$V2 <- gsub("total_cpgShelf", "CpG_Shelf", dat2$V2)
dat2$V1 <- dat2$V1 / denominator * 100
colnames(dat2) <- c("Percent", "Label")
dat2$Group <- ifelse(grepl("CpG", dat2$Label), "CPG_island", "TE_annotation")

data <- rbind(enh, dat2)
data$Label <- factor(data$Label, levels = c(
  "intron", "Intergenic", "exon", "TTS",
  "LINE", "SINE", "LTR", "DNA",
  "CpG_Shelf", "CpG_Shore", "CpG_Island"
))

pdf(file.path(annotation_dir, "CTCF_peak_annotation.pdf"), width = 3.4, height = 3)
ggplot(data, aes(x = Label, y = Percent, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("Percentage (%)") +
  theme(
    legend.text = element_text(colour = "black", size = 9),
    panel.border = element_blank(),
    legend.position = c(0.8, 0.85),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black", size = 12),
    axis.text.y = element_text(vjust = 0.5, colour = "black", size = 12)
  )
dev.off()

###############################################################################
# End R script
###############################################################################
