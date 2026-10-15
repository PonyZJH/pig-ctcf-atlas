#!/usr/bin/env bash

###############################################################################
# 1. Paths and parameters
###############################################################################
# CTCF cluster and singleton analysis.
# Bash sections should be run in bash. Blocks marked "R script" should be run in
# R from the project root.
#
# The original notes used "le" as a local shorthand for less -S. No alias is
# defined here; keep using it only if it is available in the user's shell.

project_dir="$(pwd)"

individual_bed_dir="${project_dir}/results/cbs/bed"
consensus_cbs="${project_dir}/results/cbs/1.bed"
tissue_table="" # Metadata: SampleID, Original, Filename, Organ; no header.
cluster_dir="${project_dir}/results/ctcf_cluster"

mkdir -p "${cluster_dir}"

###############################################################################
# 2. Cluster individual sample CBSs within 10 kb
###############################################################################
for bed_file in "${individual_bed_dir}"/*.bed; do
    sample=$(basename "${bed_file}" .bed)
    sample_dir="${cluster_dir}/${sample}"
    mkdir -p "${sample_dir}"

    echo "Processing ${bed_file} ..."

    bedtools cluster -i "${bed_file}" -d 10000 > "${sample_dir}/1.cluster.bed"

    bedtools groupby \
        -i "${sample_dir}/1.cluster.bed" \
        -g 4 \
        -c 1,2,3 \
        -o distinct,min,max \
        > "${sample_dir}/cluster_summary.bed"

    awk 'BEGIN{OFS="\t"}{print $2,$3,$4,"cluster_"$1}' \
        "${sample_dir}/cluster_summary.bed" \
        > "${sample_dir}/cluster.bed"

    cut -f4 "${sample_dir}/1.cluster.bed" \
        | sort | uniq -c \
        | awk 'BEGIN{OFS="\t"}{print "cluster_"$2,$1}' \
        > "${sample_dir}/cluster_count.txt"

    awk '$2==1 {print $1}' "${sample_dir}/cluster_count.txt" > "${sample_dir}/single_site.txt"
    awk '$2!=1 {print $1}' "${sample_dir}/cluster_count.txt" > "${sample_dir}/cluster_site.txt"

    awk 'NR==FNR {a[$1]; next} $4 in a' \
        "${sample_dir}/single_site.txt" "${sample_dir}/cluster.bed" \
        | cut -f -3 \
        | sort -k1,1 -k2,2n \
        > "${sample_dir}/non_cluster.bed"

    awk 'NR==FNR {a[$1]; next} $4 in a' \
        "${sample_dir}/cluster_site.txt" "${sample_dir}/cluster.bed" \
        | cut -f -3 \
        | sort -k1,1 -k2,2n \
        > "${sample_dir}/truecluster.bed"

    echo "Finished ${sample}"
done

###############################################################################
# 3. Consensus cluster set and overlap files
###############################################################################
cat "${cluster_dir}"/*/truecluster.bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${cluster_dir}/1.bed"

for sample_dir in "${cluster_dir}"/*; do
    [ -d "${sample_dir}" ] || continue
    sample=$(basename "${sample_dir}")

    bedtools intersect \
        -a "${cluster_dir}/1.bed" \
        -b "${sample_dir}/truecluster.bed" \
        -wa \
        | sort | uniq \
        > "${cluster_dir}/${sample}_overlap.bed"
done

###############################################################################
# R script: cluster 0/1 matrix and tissue aggregation
###############################################################################
project_dir <- getwd()
cluster_dir <- file.path(project_dir, "results/ctcf_cluster")
tissue_table <- "" # Metadata: SampleID, Original, Filename, Organ; no header.

cluster_bed <- read.table(file.path(cluster_dir, "1.bed"), header = FALSE)
cluster_id <- paste(cluster_bed$V1, cluster_bed$V2, cluster_bed$V3, sep = "_")

data <- data.frame(row.names = cluster_id)
overlap_files <- list.files(cluster_dir, pattern = "_overlap.bed$", full.names = TRUE)

for (path in overlap_files) {
  sample <- sub("_overlap.bed$", "", basename(path))
  if (file.info(path)$size > 0) {
    pp <- read.table(path, header = FALSE)
    present <- paste(pp$V1, pp$V2, pp$V3, sep = "_")
  } else {
    present <- character(0)
  }
  data[[sample]] <- as.integer(cluster_id %in% present)
  cat(sample, "\n")
}

write.table(data,
            file = file.path(cluster_dir, "CTCFcluster_149sample_01.txt"),
            sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

metadata <- read.table(tissue_table, header = FALSE, stringsAsFactors = FALSE)
colnames(metadata)[1:4] <- c("SampleID", "Original", "Filename", "Organ")

merged_data <- data.frame(row.names = rownames(data))
for (organ in unique(metadata$Organ)) {
  organ_files <- metadata$Filename[metadata$Organ == organ]
  organ_files <- intersect(organ_files, colnames(data))
  if (length(organ_files) == 0) {
    cat("Warning: No columns found for organ:", organ, "\n")
    next
  }
  organ_samples <- data[, organ_files, drop = FALSE]
  merged_data[[organ]] <- rowSums(organ_samples, na.rm = TRUE)
}

write.table(merged_data,
            file = file.path(cluster_dir, "merged_tissues.txt"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

tissue_output_dir <- file.path(cluster_dir, "tissue")
dir.create(tissue_output_dir, showWarnings = FALSE)

for (i in seq_len(ncol(merged_data))) {
  tissue_name <- colnames(merged_data)[i]
  tissue_data <- merged_data[, i, drop = FALSE]
  tissue_data[tissue_data == 0] <- NA
  tissue_data_non_na <- na.omit(tissue_data)
  non_empty_sites <- rownames(tissue_data_non_na)

  output_file <- file.path(tissue_output_dir, paste0(tissue_name, "_cluster.txt"))
  write.table(non_empty_sites, file = output_file,
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat("Processed:", tissue_name, "\n")
}

site_output_dir <- file.path(cluster_dir, "site")
dir.create(site_output_dir, showWarnings = FALSE)

for (sample in colnames(data)) {
  present_peaks <- rownames(data)[data[[sample]] == 1]
  output_file <- file.path(site_output_dir, paste0(sample, "_cluster.txt"))
  writeLines(present_peaks, output_file)
}
# End R script

###############################################################################
# 4. Cluster-associated CBSs and singleton CBSs
###############################################################################
bedtools intersect \
    -a "${consensus_cbs}" \
    -b "${cluster_dir}/1.bed" \
    -v \
    | sort -k1,1 -k2,2n \
    > "${cluster_dir}/singleton.bed"

bedtools intersect \
    -a "${consensus_cbs}" \
    -b "${cluster_dir}/1.bed" \
    -wa \
    | sort | uniq \
    > "${cluster_dir}/cluster_site.bed"
