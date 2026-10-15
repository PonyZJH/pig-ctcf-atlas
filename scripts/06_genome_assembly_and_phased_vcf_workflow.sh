#!/usr/bin/env bash

###############################################################################
# 1. Inputs and paths
###############################################################################
# Source: raw hifiasm notes. Fill input paths before running in Bash.
# Outputs: prefix.dip.vcf.gz (dipcall), hifi.vcf.gz (DeepVariant), and
# hifi.g.vcf.gz (the accompanying gVCF). Hap1/Hap2 parental identity must be
# established from the assembly configuration, not inferred from file names.
project_dir="$(pwd)"
father_r1="" # Paternal paired-end WGS FASTQ, read 1.
father_r2="" # Paternal paired-end WGS FASTQ, read 2.
mother_r1="" # Maternal paired-end WGS FASTQ, read 1.
mother_r2="" # Maternal paired-end WGS FASTQ, read 2.
hifi_reads="" # Offspring HiFi FASTQ.
ont_reads="" # Offspring ultra-long ONT FASTQ.
reference_fasta="" # Reference FASTA matching downstream CUT&Tag alignments.
hifi_bam="" # Existing reference-aligned, sorted HiFi BAM for DeepVariant.
paftools_js="" # minimap2 paftools.js file.
dipcall_aux_js="" # dipcall dipcall-aux.js file.
output_dir="${project_dir}/results/genome_phasing"
assembly_prefix="${output_dir}/assembly"
hap1_gfa="" # Hap1 primary-contig GFA produced by hifiasm.
hap2_gfa="" # Hap2 primary-contig GFA produced by hifiasm.
threads=10
mkdir -p "${output_dir}"

###############################################################################
# 2. Parental k-mer databases
###############################################################################
yak count -b37 -t10 -o "${output_dir}/pat.yak" \
    <(cat "${father_r1}" "${father_r2}") \
    <(cat "${father_r1}" "${father_r2}")
yak count -b37 -t10 -o "${output_dir}/mat.yak" \
    <(cat "${mother_r1}" "${mother_r2}") \
    <(cat "${mother_r1}" "${mother_r2}")

###############################################################################
# 3. Trio-assisted haplotype assembly and FASTA conversion
###############################################################################
hifiasm -o "${assembly_prefix}" -t 10 \
    -1 "${output_dir}/pat.yak" -2 "${output_dir}/mat.yak" \
    --ul "${ont_reads}" "${hifi_reads}"

# Apply the original GFA-to-FASTA awk example to each haplotype GFA.
awk '/^S/{print ">"$2;print $3}' "${hap1_gfa}" | gzip > "${output_dir}/hap1.fa.gz"
awk '/^S/{print ">"$2;print $3}' "${hap2_gfa}" | gzip > "${output_dir}/hap2.fa.gz"

###############################################################################
# 4. Assembly-based dipcall variants and callable intervals
###############################################################################
minimap2 -c --paf-no-hit -xasm5 --cs -t"${threads}" \
    "${reference_fasta}" "${output_dir}/hap1.fa.gz" > "${output_dir}/prefix.hap1.paf"
gzip "${output_dir}/prefix.hap1.paf"
minimap2 -c --paf-no-hit -xasm5 --cs -t"${threads}" \
    "${reference_fasta}" "${output_dir}/hap2.fa.gz" > "${output_dir}/prefix.hap2.paf"
gzip "${output_dir}/prefix.hap2.paf"

gzip -dc "${output_dir}/prefix.hap1.paf.gz" | sort -k6,6 -k8,8n \
    | k8 "${paftools_js}" call - 2> "${output_dir}/prefix.hap1.var.gz.vst" \
    | gzip > "${output_dir}/prefix.hap1.var.gz"
gzip -dc "${output_dir}/prefix.hap2.paf.gz" | sort -k6,6 -k8,8n \
    | k8 "${paftools_js}" call - 2> "${output_dir}/prefix.hap2.var.gz.vst" \
    | gzip > "${output_dir}/prefix.hap2.var.gz"

gzip -dc "${output_dir}/prefix.hap1.var.gz" | grep '^R' | cut -f2- > "${output_dir}/prefix.hap1.bed"
gzip -dc "${output_dir}/prefix.hap2.var.gz" | grep '^R' | cut -f2- > "${output_dir}/prefix.hap2.bed"
bedtk isec -m "${output_dir}/prefix.hap1.bed" "${output_dir}/prefix.hap2.bed" > "${output_dir}/prefix.dip.bed"

minimap2 -a -xasm5 --cs -t"${threads}" "${reference_fasta}" "${output_dir}/hap1.fa.gz" \
    | gzip > "${output_dir}/prefix.hap1.sam.gz"
minimap2 -a -xasm5 --cs -t"${threads}" "${reference_fasta}" "${output_dir}/hap2.fa.gz" \
    | gzip > "${output_dir}/prefix.hap2.sam.gz"

k8 "${dipcall_aux_js}" samflt "${output_dir}/prefix.hap1.sam.gz" \
    | samtools sort -m4G --threads 4 -o "${output_dir}/prefix.hap1.bam" -
k8 "${dipcall_aux_js}" samflt "${output_dir}/prefix.hap2.sam.gz" \
    | samtools sort -m4G --threads 4 -o "${output_dir}/prefix.hap2.bam" -
htsbox pileup -q5 -evcf "${reference_fasta}" "${output_dir}/prefix.hap1.bam" "${output_dir}/prefix.hap2.bam" \
    | htsbox bgzip > "${output_dir}/prefix.pair.vcf.gz"
k8 "${dipcall_aux_js}" vcfpair "${output_dir}/prefix.pair.vcf.gz" \
    | htsbox bgzip > "${output_dir}/prefix.dip.vcf.gz"
ls -lh "${output_dir}/prefix.dip.vcf.gz" "${output_dir}/prefix.dip.bed"

###############################################################################
# 5. HiFi-based variant calling
###############################################################################
# The original notes do not provide the HiFi alignment command.
run_deepvariant --model_type="PACBIO" --ref="${reference_fasta}" \
    --reads="${hifi_bam}" --output_vcf="${output_dir}/hifi.vcf.gz" \
    --output_gvcf="${output_dir}/hifi.g.vcf.gz" --num_shards=20
