#!/usr/bin/env bash

###############################################################################
# 1. Paths and parameters
###############################################################################
# Haplotype-resolved CTCF binding, HS-CBS and ASB analysis.
# Bash sections should be run in bash. Blocks marked "R script" should be run in
# R from the project root.
#
# The original notes used "le" as a local shorthand for less -S. No alias is
# defined here; keep using it only if it is available in the user's shell.
#
# haplotype_sample_table columns:
#   sample_id  bam_prefix
#
# phased_variant_bed columns:
#   chr  start  end
#
# phased_variant_allele_table columns:
#   chr  position  position  ref  alt

project_dir="$(pwd)"

sample_table="" # Haplotype sample table: sample_id, bam_prefix; no header.
bam_dir="${project_dir}/results/bam"
consensus_cbs="${project_dir}/results/cbs/1.bed"

phase_dir="${project_dir}/results/haplotype_asb"
phase_bam_dir="${phase_dir}/bam"
phase_log_dir="${phase_dir}/logs"
phase_qc_dir="${phase_dir}/qc"
hs_depth_dir="${phase_dir}/hs_cbs/depth"
hs_result_dir="${phase_dir}/hs_cbs/result"
asb_input_dir="${phase_dir}/asb/input"
asb_depth_dir="${phase_dir}/asb/depth"
asb_result_dir="${phase_dir}/asb/result"

reference_fasta="" # Study reference genome FASTA matching BAM and variant coordinates.
phased_vcf="" # Phased VCF for WhatsHap haplotag.
phased_variant_bed="" # Variant BED: chromosome, 0-based start, exclusive end.
phased_variant_allele_table="" # Allele table: chr, position, position, ref, alt; see notes for allele orientation.
jaspar_all_meme="" # JASPAR 2024 CORE vertebrate motifs in MEME format (879 motifs in the original analysis).

threads=5
tmp_root="${TMPDIR:-${phase_dir}/tmp}"

mkdir -p "${phase_bam_dir}" "${phase_log_dir}" "${phase_qc_dir}" \
    "${hs_depth_dir}" "${hs_result_dir}" \
    "${asb_input_dir}" "${asb_depth_dir}" "${asb_result_dir}" \
    "${tmp_root}"

###############################################################################
# 2. WhatsHap haplotag and split H1/H2 BAM files
###############################################################################
while read -r sample_id bam_prefix; do
    [ -z "${sample_id}" ] && continue
    [ "${sample_id#\#}" != "${sample_id}" ] && continue

    work_dir="${tmp_root}/${sample_id}_haplotag"
    mkdir -p "${work_dir}"

    whatshap haplotag \
        -o "${work_dir}/${sample_id}_haplotagged.bam" \
        --ignore-read-groups \
        --reference "${reference_fasta}" \
        "${phased_vcf}" \
        "${bam_dir}/${bam_prefix}_bowtie2.mapped.bam" \
        > "${phase_log_dir}/${sample_id}.whatshap.log" 2>&1

    samtools view -H "${bam_dir}/${bam_prefix}_bowtie2.mapped.bam" > "${work_dir}/header.sam"
    samtools view -@ "${threads}" -h "${work_dir}/${sample_id}_haplotagged.bam" \
        | awk -F '\t' '$0 ~ "HP:i:1"' \
        > "${work_dir}/H1"
    samtools view -@ "${threads}" -h "${work_dir}/${sample_id}_haplotagged.bam" \
        | awk -F '\t' '$0 ~ "HP:i:2"' \
        > "${work_dir}/H2"

    cat "${work_dir}/header.sam" "${work_dir}/H1" > "${work_dir}/H1_head"
    cat "${work_dir}/header.sam" "${work_dir}/H2" > "${work_dir}/H2_head"

    samtools view -bS "${work_dir}/H1_head" > "${phase_bam_dir}/${sample_id}_haplotagged_H1.bam"
    samtools view -bS "${work_dir}/H2_head" > "${phase_bam_dir}/${sample_id}_haplotagged_H2.bam"
    cp "${work_dir}/${sample_id}_haplotagged.bam" "${phase_bam_dir}/"

    rm -rf "${work_dir}"
done < "${sample_table}"

###############################################################################
# 3. Phasing efficiency and haplotype read-count QC
###############################################################################
grep -h "Total alignments processed" "${phase_log_dir}"/*.whatshap.log \
    | awk '{print $4}' \
    > "${phase_qc_dir}/total_alignments_processed.txt"

grep -h "Alignments that could be tagged" "${phase_log_dir}"/*.whatshap.log \
    | awk '{print $6}' \
    > "${phase_qc_dir}/alignments_that_could_be_tagged.txt"

ls "${phase_bam_dir}"/*_H1.bam | while read -r bam; do
    samtools flagstat "${bam}" >> "${phase_qc_dir}/total_H1_reads.txt"
done

ls "${phase_bam_dir}"/*_H2.bam | while read -r bam; do
    samtools flagstat "${bam}" >> "${phase_qc_dir}/total_H2_reads.txt"
done

le "${phase_qc_dir}/total_H1_reads.txt" \
    | grep total \
    | awk '{print $1}' \
    > "${phase_qc_dir}/total_H1_reads.txt2"

le "${phase_qc_dir}/total_H2_reads.txt" \
    | grep total \
    | awk '{print $1}' \
    > "${phase_qc_dir}/total_H2_reads.txt2"

###############################################################################
# R script: phasing efficiency and H1/H2 read-count plots
###############################################################################
library(ggplot2)

project_dir <- getwd()
phase_qc_dir <- file.path(project_dir, "results/haplotype_asb/qc")

used <- read.table(file.path(phase_qc_dir, "total_alignments_processed.txt"))
used$mean <- mean(used$V1)
used$se <- sd(used$V1) / sqrt(nrow(used))
used$group <- "Used_reads"

phased <- read.table(file.path(phase_qc_dir, "alignments_that_could_be_tagged.txt"))
phased$mean <- mean(phased$V1)
phased$se <- sd(phased$V1) / sqrt(nrow(phased))
phased$group <- "Phased_reads"

dat <- rbind(phased, used)

pdf(file.path(phase_qc_dir, "phased_efficiency.pdf"), width = 6, height = 3)
ggplot(dat, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.7), width = 0.15) +
  scale_fill_manual(values = c("#f49128", "#194a55")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        panel.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.title = element_text(size = 16, color = "black"),
        axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(size = 16, color = "black")) +
  labs(x = "", y = "Read Counts")
dev.off()

h1 <- read.table(file.path(phase_qc_dir, "total_H1_reads.txt2"))
h2 <- read.table(file.path(phase_qc_dir, "total_H2_reads.txt2"))
dat <- data.frame(H1_reads = h1$V1, H2_reads = h2$V1)

pdf(file.path(phase_qc_dir, "phased_bias.pdf"), width = 4, height = 3)
ggplot(data = dat, aes(x = H1_reads, y = H2_reads)) +
  geom_point(size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        panel.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.title = element_text(size = 16, color = "black", vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, color = "black", size = 16),
        axis.text.y = element_text(size = 16, color = "black")) +
  labs(x = "H1_reads", y = "H2_reads")
dev.off()
# End R script

###############################################################################
# 4. Haplotype read counts at consensus CBSs
###############################################################################
while read -r sample_id bam_prefix; do
    [ -z "${sample_id}" ] && continue
    [ "${sample_id#\#}" != "${sample_id}" ] && continue

    samtools index "${phase_bam_dir}/${sample_id}_haplotagged_H1.bam"
    samtools index "${phase_bam_dir}/${sample_id}_haplotagged_H2.bam"

    multiBamSummary BED-file \
        --BED "${consensus_cbs}" \
        -b "${phase_bam_dir}/${sample_id}_haplotagged_H1.bam" \
           "${phase_bam_dir}/${sample_id}_haplotagged_H2.bam" \
        -out "${hs_depth_dir}/${sample_id}_deeptools.npz" \
        --outRawCounts "${hs_depth_dir}/${sample_id}_deeptools.txt"
done < "${sample_table}"

###############################################################################
# R script: HS-CBS binomial test, FDR filtering and direction consistency
###############################################################################
project_dir <- getwd()
hs_depth_dir <- file.path(project_dir, "results/haplotype_asb/hs_cbs/depth")
hs_result_dir <- file.path(project_dir, "results/haplotype_asb/hs_cbs/result")
dir.create(hs_result_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(hs_depth_dir, pattern = "_deeptools.txt$", full.names = TRUE)
all_threshold <- data.frame()

for (path in files) {
  raw <- read.table(path, header = TRUE, check.names = FALSE)
  colnames(raw) <- sub("^#", "", colnames(raw))
  dat <- data.frame(
    chr = raw[[1]],
    start = raw[[2]],
    end = raw[[3]],
    H1_count = as.numeric(raw[[4]]),
    H2_count = as.numeric(raw[[5]])
  )
  dat <- dat[dat$H1_count >= 10 | dat$H2_count >= 10, ]
  dat$total <- dat$H1_count + dat$H2_count
  dat$ratio <- dat$H1_count / dat$total
  dat$pvalue <- vapply(seq_len(nrow(dat)), function(i) {
    binom.test(dat$H1_count[i], dat$total[i], 0.5)$p.value
  }, numeric(1))
  dat$name <- sub("_deeptools.txt$", "", basename(path))
  dat$qvalue <- p.adjust(dat$pvalue, method = "BH")

  qvalue_file <- file.path(hs_depth_dir, paste0(dat$name[1], "_qvalue.txt"))
  write.table(dat, file = qvalue_file, col.names = TRUE, row.names = FALSE, quote = FALSE)

  threshold <- dat[(dat$ratio >= 0.65 | dat$ratio <= 0.35) & dat$qvalue < 0.05, ]
  threshold_file <- file.path(hs_depth_dir, paste0(dat$name[1], "_threshold.txt"))
  write.table(threshold, file = threshold_file, col.names = TRUE, row.names = FALSE, quote = FALSE)
  all_threshold <- rbind(all_threshold, threshold)
}

write.table(all_threshold, file = file.path(hs_result_dir, "CTCF_AS_peak_threshold_total.txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE)

if (nrow(all_threshold) > 0) {
  all_threshold$loci <- paste(all_threshold$chr, all_threshold$start, all_threshold$end, sep = "_")
  all_threshold$direction <- ifelse(all_threshold$ratio > 0.6, "H1", "H2")

  N <- NULL
  for (locus in unique(all_threshold$loci)) {
    tmp <- all_threshold[all_threshold$loci == locus, ]
    h1_fraction <- nrow(tmp[tmp$ratio >= 0.65, ]) / nrow(tmp)
    if (h1_fraction > 0.5) {
      tmp$percentage <- h1_fraction
    } else {
      tmp$percentage <- 1 - h1_fraction
    }
    N <- rbind(N, tmp)
  }

  write.table(N, file = file.path(hs_result_dir, "CTCF_AS_peak_direction.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)

  consistent <- N[N$percentage == 1, ]
  write.table(consistent, file = file.path(hs_result_dir, "CTCF_AS_peak_direction_consistent.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)

  consistent_bed <- unique(consistent[, c("chr", "start", "end", "direction")])
  consistent_bed <- consistent_bed[order(consistent_bed$chr, consistent_bed$start), ]
  write.table(consistent_bed,
              file = file.path(hs_result_dir, "CTCF_AS_peak_direction_consistent.bed"),
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

  consistent_uniq <- unique(consistent[, c("chr", "start", "end")])
  consistent_uniq <- consistent_uniq[order(consistent_uniq$chr, consistent_uniq$start), ]
  write.table(consistent_uniq,
              file = file.path(hs_result_dir, "CTCF_AS_peak_direction_consistent_uniq.bed"),
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

threshold_files <- list.files(hs_depth_dir, pattern = "_threshold.txt$", full.names = TRUE)
result <- NULL
for (path in threshold_files) {
  pp <- read.table(path, header = TRUE)
  if (nrow(pp) == 0) next
  rownames(pp) <- paste(pp$chr, pp$start, pp$end, sep = "_")
  pp <- data.frame(V1 = rownames(pp), present = 1)
  colnames(pp)[2] <- sub("_threshold.txt$", "", basename(path))
  if (is.null(result)) {
    result <- pp
  } else {
    result <- merge(result, pp, by = "V1", all = TRUE)
  }
}

if (!is.null(result)) {
  result[is.na(result)] <- 0
  rownames(result) <- result$V1
  result <- result[, 2:ncol(result), drop = FALSE]
  write.table(result, file = file.path(hs_result_dir, "CTCF_AS_peak_01file.txt"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
}
# End R script

###############################################################################
# 5. ASB target variants within HS-CBSs
###############################################################################
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' \
    "${hs_result_dir}/CTCF_AS_peak_direction_consistent_uniq.bed" \
    | sort | uniq \
    > "${hs_result_dir}/CTCF_AS_peak_total.bed"

bedtools intersect \
    -a "${phased_variant_bed}" \
    -b "${hs_result_dir}/CTCF_AS_peak_total.bed" \
    -wa \
    | sort | uniq \
    > "${asb_input_dir}/target_variants.bed"

awk 'BEGIN{OFS=""}{print $1,":",$2+1,"-",$3}' \
    "${asb_input_dir}/target_variants.bed" \
    > "${asb_input_dir}/target_variants.regions.txt"

###############################################################################
# 6. SAMtools mpileup at candidate ASB variants
###############################################################################
while read -r sample_id bam_prefix; do
    [ -z "${sample_id}" ] && continue
    [ "${sample_id#\#}" != "${sample_id}" ] && continue

    samtools index "${phase_bam_dir}/${sample_id}_haplotagged_H1.bam"
    samtools index "${phase_bam_dir}/${sample_id}_haplotagged_H2.bam"
    samtools index "${phase_bam_dir}/${sample_id}_haplotagged.bam"

    : > "${asb_depth_dir}/${sample_id}_depth.txt"
    while read -r region; do
        samtools mpileup \
            -f "${reference_fasta}" \
            -r "${region}" \
            "${phase_bam_dir}/${sample_id}_haplotagged_H1.bam" \
            "${phase_bam_dir}/${sample_id}_haplotagged_H2.bam" \
            >> "${asb_depth_dir}/${sample_id}_depth.txt"
    done < "${asb_input_dir}/target_variants.regions.txt"
done < "${sample_table}"

cd "${asb_depth_dir}" || exit 1
ls *depth.txt | while read -r id; do
    le "${id}" \
        | awk '{print $1" "$2" "$3" "$4" "$5}' \
        | awk '{if($4>=10 || $5>=10)print $0}' \
        > "${id%%.txt*}_filter.txt"
done
cd "${project_dir}" || exit 1

###############################################################################
# R script: ASB binomial test, FDR filtering and direction consistency
###############################################################################
project_dir <- getwd()
asb_depth_dir <- file.path(project_dir, "results/haplotype_asb/asb/depth")
asb_result_dir <- file.path(project_dir, "results/haplotype_asb/asb/result")
dir.create(asb_result_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(asb_depth_dir, pattern = "_filter.txt$", full.names = TRUE)
all_threshold <- data.frame()

for (path in files) {
  data2 <- read.table(path, header = FALSE)
  if (nrow(data2) == 0) next

  colnames(data2)[1:5] <- c("chr", "start", "end", "H1_count", "H2_count")
  data2$total <- data2$H1_count + data2$H2_count
  data2$ratio <- data2$H1_count / data2$total
  data2$pvalue <- vapply(seq_len(nrow(data2)), function(i) {
    binom.test(data2$H1_count[i], data2$total[i], 0.5)$p.value
  }, numeric(1))
  data2$name <- sub("_filter.txt$", "", basename(path))
  data2$qvalue <- p.adjust(data2$pvalue, method = "BH")

  qvalue_file <- file.path(asb_depth_dir, paste0(data2$name[1], "_qvalue.txt"))
  write.table(data2, file = qvalue_file, col.names = TRUE, row.names = FALSE, quote = FALSE)

  threshold <- data2[(data2$ratio >= 0.65 | data2$ratio <= 0.35) & data2$qvalue < 0.05, ]
  threshold_file <- file.path(asb_depth_dir, paste0(data2$name[1], "_threshold.txt"))
  write.table(threshold, file = threshold_file, col.names = TRUE, row.names = FALSE, quote = FALSE)
  all_threshold <- rbind(all_threshold, threshold)
}

write.table(all_threshold, file = file.path(asb_result_dir, "CTCF_ASB_threshold_total.txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE)

if (nrow(all_threshold) > 0) {
  all_threshold$loci <- paste0(all_threshold$chr, "_", all_threshold$start)
  all_threshold$direction <- ifelse(all_threshold$ratio > 0.6, "H1", "H2")

  N <- NULL
  for (locus in unique(all_threshold$loci)) {
    tmp <- all_threshold[all_threshold$loci == locus, ]
    h1_fraction <- nrow(tmp[tmp$ratio >= 0.65, ]) / nrow(tmp)
    if (h1_fraction > 0.5) {
      tmp$percentage <- h1_fraction
    } else {
      tmp$percentage <- 1 - h1_fraction
    }
    N <- rbind(N, tmp)
  }

  write.table(N, file = file.path(asb_result_dir, "CTCF_ASB_direction.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)

  consistent <- N[N$percentage == 1, ]
  write.table(consistent, file = file.path(asb_result_dir, "CTCF_ASB_direction_consistent.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)

  consistent_direction <- unique(data.frame(
    chr = consistent$chr,
    start = consistent$start,
    end = consistent$start,
    direction = consistent$direction
  ))
  consistent_direction <- consistent_direction[order(consistent_direction$chr, consistent_direction$start), ]
  write.table(consistent_direction,
              file = file.path(asb_result_dir, "CTCF_ASB_direction_consistent_direction.bed"),
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}
# End R script

###############################################################################
# 7. Keep ASB variants consistent with haplotype-specific CBS direction
###############################################################################
bedtools intersect \
    -a "${asb_result_dir}/CTCF_ASB_direction_consistent_direction.bed" \
    -b "${hs_result_dir}/CTCF_AS_peak_direction_consistent.bed" \
    -wa -wb \
    | awk '$4==$8' \
    > "${asb_result_dir}/AS_peak_ASB.txt"

awk 'BEGIN{OFS="\t"}{print $1,$2,$2}' "${asb_result_dir}/AS_peak_ASB.txt" \
    | sort | uniq \
    > "${asb_result_dir}/ASB_24bp_loci.bed"

###############################################################################
# 8. Ref/alt 49 bp windows centered on ASB variants
###############################################################################
cd "${asb_result_dir}" || exit 1

bedtools intersect \
    -a ASB_24bp_loci.bed \
    -b "${phased_variant_allele_table}" \
    -wa -wb \
    | awk '$2==$5 {print $1"\t"$2"\t"$3"\t"$7"\t"$8}' \
    > ASB_24bp_loci.fa

awk 'length($4)==1 && length($5)==1' ASB_24bp_loci.fa > ASB_24bp_loci2.fa

awk 'BEGIN{OFS="\t"}{print $1,$2-25,$2-1}' ASB_24bp_loci2.fa > ASB_24bp_tmp1.bed
awk 'BEGIN{OFS="\t"}{print $1,$2,$2+24}' ASB_24bp_loci2.fa > ASB_24bp_tmp2.bed

bedtools getfasta -fi "${reference_fasta}" -bed ASB_24bp_tmp1.bed -fo ASB_24bp_tmp1.fa
bedtools getfasta -fi "${reference_fasta}" -bed ASB_24bp_tmp2.bed -fo ASB_24bp_tmp2.fa

grep -v ">" ASB_24bp_tmp1.fa > ASB_24bp_tmp1.fa2
grep -v ">" ASB_24bp_tmp2.fa > ASB_24bp_tmp2.fa2

paste ASB_24bp_tmp1.fa2 ASB_24bp_loci2.fa ASB_24bp_tmp2.fa2 \
    | awk '{print $1""$5""$7}' \
    > ASB_49bp_ref.fa

paste ASB_24bp_tmp1.fa2 ASB_24bp_loci2.fa ASB_24bp_tmp2.fa2 \
    | awk '{print $1""$6""$7}' \
    > ASB_49bp_alt.fa

paste ASB_24bp_tmp1.bed ASB_24bp_tmp2.bed \
    | awk '{print ">"$1":"$2"-"$6}' \
    > head_name

paste head_name ASB_49bp_ref.fa | awk '{print $1;print $2}' > ASB_49bp_ref_head.fa
paste head_name ASB_49bp_alt.fa | awk '{print $1;print $2}' > ASB_49bp_alt_head.fa

ls *head.fa | while read -r fasta; do
    fimo \
        --oc "./${fasta%%_head.fa*}_fimo_output" \
        --parse-genomic-coord \
        "${jaspar_all_meme}" \
        "${fasta}"
done

le ASB_49bp_alt_fimo_output/fimo.tsv \
    | awk '{print $3"\t"$4"\t"$5"\t"$0}' \
    | grep chr \
    > ASB_49bp_alt_result.bed

le ASB_49bp_ref_fimo_output/fimo.tsv \
    | awk '{print $3"\t"$4"\t"$5"\t"$0}' \
    | grep chr \
    > ASB_49bp_ref_result.bed

bedtools intersect \
    -a ASB_24bp_loci2.fa \
    -b ASB_49bp_ref_result.bed \
    -wa -wb \
    | awk '{print $1"_"$2"_"$4"_"$5"\t"$6"_"$7"_"$8"\t"$9"\t"$10"\t"$14"\t"$15}' \
    | sort | uniq \
    > ASB_49bp_ref_result_overlap_ASB.bed

bedtools intersect \
    -a ASB_24bp_loci2.fa \
    -b ASB_49bp_alt_result.bed \
    -wa -wb \
    | awk '{print $1"_"$2"_"$4"_"$5"\t"$6"_"$7"_"$8"\t"$9"\t"$10"\t"$14"\t"$15}' \
    | sort | uniq \
    > ASB_49bp_alt_result_overlap_ASB.bed

cd "${project_dir}" || exit 1

###############################################################################
# R script: motif disruption summary and ASB plots
###############################################################################
library(dplyr)
library(tidyr)
library(ggplot2)

project_dir <- getwd()
asb_result_dir <- file.path(project_dir, "results/haplotype_asb/asb/result")

ref <- read.table(file.path(asb_result_dir, "ASB_49bp_ref_result_overlap_ASB.bed"))
alt <- read.table(file.path(asb_result_dir, "ASB_49bp_alt_result_overlap_ASB.bed"))

ref$name <- paste(ref$V1, ref$V2, ref$V5, ref$V3, ref$V4, sep = "_")
alt$name <- paste(alt$V1, alt$V2, alt$V5, alt$V3, alt$V4, sep = "_")

ref_filtered <- ref %>%
  group_by(name) %>%
  filter(V6 == max(V6)) %>%
  ungroup() %>%
  as.data.frame()

alt_filtered <- alt %>%
  group_by(name) %>%
  filter(V6 == max(V6)) %>%
  ungroup() %>%
  as.data.frame()

data <- merge(ref_filtered, alt_filtered, by = "name")
data$FC <- data$V6.x - data$V6.y
dat <- data[abs(data$FC) > 3, c(1, 7, 13, 14)]
colnames(dat)[2:3] <- c("ref_TF_score", "alt_TF_score")

dat2 <- separate(dat, col = name,
                 into = c("chr", "loci", "ref", "alt", "Chr", "start", "end", "strand", "motif", "TF"),
                 sep = "_")

ASB <- read.table(file.path(asb_result_dir, "AS_peak_ASB.txt"))
dat2 <- dat2 %>% mutate(V1 = as.character(chr), V2 = as.character(loci))
ASB <- ASB %>% mutate(V1 = as.character(V1), V2 = as.character(V2))
N <- inner_join(dat2, ASB, by = c("V1", "V2"))

dat4 <- N[, c(1:4, 6:13, 17, 19:20)]
colnames(dat4)[13] <- "direction"
colnames(dat4)[14:15] <- c("AS_peak_start", "AS_peak_end")
dat4$name <- paste0(dat4$chr, "_", dat4$loci)

ASB_ratio <- read.table(file.path(asb_result_dir, "CTCF_ASB_direction_consistent.txt"), header = TRUE)
ASB_ratio$loci <- as.character(ASB_ratio$loci)
ASB_ratio_summary <- ASB_ratio %>%
  group_by(loci) %>%
  summarise(ASB_ref_read_mean = mean(H1_count, na.rm = TRUE),
            ASB_alt_read_mean = mean(H2_count, na.rm = TRUE),
            ASB_total_read_mean = mean(total, na.rm = TRUE),
            .groups = "drop") %>%
  as.data.frame()
colnames(ASB_ratio_summary)[1] <- "name"

dat8 <- merge(dat4, ASB_ratio_summary, by = "name")
dat8$ratio <- dat8$ASB_ref_read_mean / dat8$ASB_total_read_mean

dat8 <- dat8 %>%
  mutate(threshold = case_when(
    FC > 0 & ratio > 0.5 ~ "Hap1",
    FC < 0 & ratio < 0.5 ~ "Hap2",
    TRUE ~ "Stable"
  ))

dat8$threshold <- factor(dat8$threshold, levels = c("Hap1", "Hap2", "Stable"))
dat9 <- dat8
dat9$threshold <- gsub("Hap2", "Paternal", dat9$threshold)
dat9$threshold <- gsub("Hap1", "Maternal", dat9$threshold)
dat_colored <- dat9 %>% filter(threshold != "Stable")

png(file.path(asb_result_dir, "TF_score_ASB_score.png"), width = 2000, height = 2000, res = 300)
p <- ggplot(dat9, aes(x = ratio, y = FC, color = threshold)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("#DC143C", "#00008B", "#808080")) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), breaks = seq(-22, 20, 5),
                     labels = as.character(abs(seq(-22, 20, 5))), limits = c(-22, 20)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.1),
                     labels = as.character(seq(0, 1, 0.1)), limits = c(0, 1)) +
  theme(legend.text = element_text(colour = "black", size = 12),
        panel.border = element_blank(),
        legend.position = c(0.9, 0.95),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 13)) +
  ylab("Difference between maternal TF score and paternal TF score") +
  xlab("Average (Maternal read counts / total read counts)") +
  geom_vline(xintercept = 0.5, lty = 3, col = "black", lwd = 0.5) +
  geom_hline(yintercept = 0, lty = 3, col = "black", lwd = 0.5)

p <- p + geom_smooth(data = dat_colored, aes(x = ratio, y = FC),
                     method = "lm", se = FALSE, color = "black", size = 0.8)

cor_test <- cor.test(dat_colored$ratio, dat_colored$FC)
label_text <- paste0("r = ", round(cor_test$estimate, 2),
                     "\n", "p = ", format.pval(cor_test$p.value, digits = 2, scientific = TRUE))
p <- p + annotate("text", x = 0.65, y = 18, label = label_text,
                  size = 5, color = "black", hjust = 0)
print(p)
dev.off()

a <- dat8[(dat8$direction == "H2" & dat8$FC < 0), ]
b <- dat8[(dat8$direction == "H1" & dat8$FC > 0), ]
dat5 <- rbind(a, b)
dat6 <- dat5[(dat5$ref_TF_score > 0 | dat5$alt_TF_score > 0), ]
write.table(dat6, file.path(asb_result_dir, "TF_score_ASB_score.txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE)

dat <- unique(data.frame(name = dat6$name, TF = dat6$TF))
TF_per_loci <- as.data.frame(table(dat$name))
loci_per_TF <- as.data.frame(table(dat$TF))

pdf(file.path(asb_result_dir, "TF_per_loci.pdf"), height = 3, width = 5)
ggplot(TF_per_loci, aes(x = Freq)) +
  geom_histogram(binwidth = 0.5) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Number of allelic TFs per ASB", y = "Count")
dev.off()

pdf(file.path(asb_result_dir, "loci_per_TF.pdf"), height = 3, width = 3.5)
ggplot(loci_per_TF, aes(x = Freq)) +
  geom_histogram(binwidth = 0.5) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Number of allelic binding sites per TF", y = "Count")
dev.off()

ctcf <- dat6[grepl("CTCF.*", dat6$TF), ]

for (motif_id in unique(ctcf$motif)) {
  tmp <- ctcf[ctcf$motif == as.character(motif_id), ]
  tmp$counts <- ifelse(tmp$strand == "+", tmp$loci - tmp$start, tmp$end - tmp$loci)
  names_df <- data.frame(Var1 = 0:max(tmp$end - tmp$start))
  counts <- as.data.frame(table(tmp$counts))
  counts2 <- merge(names_df, counts, by = "Var1", all = TRUE)
  counts2[is.na(counts2)] <- 0
  counts2$prop <- (counts2$Freq / sum(counts2$Freq)) * 100

  p <- ggplot(counts2, aes(x = Var1, y = prop)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "red", size = 2) +
    scale_x_continuous(breaks = 0:max(tmp$end - tmp$start)) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_classic() +
    xlab("") +
    ylab("Percentage")
  ggsave(file.path(asb_result_dir, paste0("motif_", motif_id, ".pdf")),
         plot = p, width = 10, height = 4)

  tmp2 <- tmp[, c("counts", "direction", "FC")]
  summ <- tmp2 %>%
    group_by(counts) %>%
    summarise(n = n(),
              mean = mean(abs(FC), na.rm = TRUE),
              sd = sd(abs(FC), na.rm = TRUE),
              se = sd / sqrt(n),
              .groups = "drop")

  p <- ggplot(summ, aes(x = counts, y = mean)) +
    geom_line(size = 1, color = "#386CB0") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.35, size = 0.35, color = "#386CB0") +
    geom_point(size = 2, color = "red") +
    labs(x = "Position", y = "Differential CTCF binding score") +
    scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10),
                       expand = expansion(mult = c(0, 0.05))) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(size = 0.35),
          axis.line = element_line(size = 0.35),
          axis.ticks = element_line(size = 0.35))
  ggsave(file.path(asb_result_dir, paste0("motif_", motif_id, "_effect.pdf")),
         plot = p, width = 10, height = 4)

  tmp3 <- tmp[, c("counts", "ref", "alt")]
  df <- tmp3 %>%
    mutate(mutation_type = case_when(
      (ref == "A" & alt == "G") | (ref == "G" & alt == "A") |
        (ref == "C" & alt == "T") | (ref == "T" & alt == "C") ~ "Transition",
      (ref == "A" & alt == "C") | (ref == "C" & alt == "A") |
        (ref == "A" & alt == "T") | (ref == "T" & alt == "A") |
        (ref == "G" & alt == "C") | (ref == "C" & alt == "G") |
        (ref == "G" & alt == "T") | (ref == "T" & alt == "G") ~ "Transversion",
      TRUE ~ NA_character_
    ))

  df_summary <- df %>%
    filter(!is.na(mutation_type)) %>%
    group_by(counts, mutation_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(counts) %>%
    mutate(perc = n / sum(n) * 100) %>%
    ungroup()

  p <- ggplot(df_summary, aes(x = factor(counts), y = perc, fill = mutation_type)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(round(perc, 1), "%")),
              position = position_stack(vjust = 0.5), size = 4) +
    scale_y_continuous(expand = expansion(c(0, 0))) +
    scale_fill_manual(values = c("Transition" = "#E41A1C", "Transversion" = "#377EB8")) +
    labs(x = "Position", y = "Percentage (%)", fill = "Mutation type") +
    theme_classic(base_size = 14)
  ggsave(file.path(asb_result_dir, paste0("motif_", motif_id, "_mutation.pdf")),
         plot = p, width = 10, height = 4)
}
# End R script
