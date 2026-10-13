# pig-ctcf-atlas

Analysis code for a haplotype-resolved, tissue-wide atlas of CTCF binding in pigs.

## 1. Introduction

This repository contains the analysis code and computational workflows developed to construct a comprehensive, haplotype-resolved atlas of CTCF binding across diverse pig tissues.

In this study, we profiled CTCF binding across 149 anatomical sites representing 57 tissues from two littermate pigs. By integrating CUT&Tag sequencing, trio-based haplotype phasing, Hi-C, and genetic association data, we systematically characterized the genomic distribution, tissue specificity, allelic variation, and regulatory functions of CTCF binding sites (CBSs).

Our analyses provide insights into tissue-specific transcriptional regulation, haplotype-dependent CTCF occupancy, and the potential contributions of CTCF-mediated chromatin organization to complex traits in pigs.

## 2. Main content

### 2.1. Genome assembly and haplotype phasing

- Trio-based haplotype-resolved genome assembly using parental short-read sequencing and offspring long-read sequencing data.
- Identification of paternal and maternal haplotypes using trio-binning approaches.
- Detection of genetic variants and construction of phased variant datasets.
- Assignment of sequencing reads to parental haplotypes for allele-specific analyses.

### 2.2. Construction of the tissue-wide CTCF binding atlas

- CUT&Tag sequencing data processing, read alignment, quality control, and peak calling.
- Integration of CTCF binding peaks across 149 anatomical sites to construct a consensus atlas.
- Quantification and normalization of CTCF binding signals across tissues.
- Characterization of the genomic distribution and functional annotations of CBSs.
- Evaluation of CTCF binding patterns and relationships among tissues through clustering analyses.

### 2.3. Tissue-conserved and tissue-specific CTCF binding

- Identification and classification of tissue-conserved and tissue-specific CBSs.
- Comparison of CTCF binding intensities and genomic distributions across tissues.
- Characterization of canonical and noncanonical CTCF binding motifs.
- Identification of tissue-enriched transcription factor motifs associated with tissue-specific CBSs.
- Investigation of potential cooperative binding between CTCF and tissue-specific transcription factors.

### 2.4. Haplotype-specific CTCF binding analysis

- Assignment of CUT&Tag sequencing reads to paternal and maternal haplotypes.
- Identification of haplotype-specific CTCF binding sites (HS-CBSs).
- Detection of allele-specific binding (ASB) variants using statistical tests of allelic imbalance.
- Characterization of sequence variations associated with haplotype-dependent CTCF occupancy.
- Investigation of the potential regulatory consequences of haplotype-specific CTCF binding.

### 2.5. Integration of CTCF binding with epigenomic and chromatin interaction data

- Integration of CTCF binding profiles with RAD21, H3K4me3, and H3K27ac CUT&Tag data.
- Characterization of CTCF binding at promoters, enhancers, and other regulatory elements.
- Identification of chromatin loops and topologically associating domain (TAD) boundaries using Hi-C data.
- Investigation of the relationships between CTCF binding, chromatin architecture, and tissue-specific gene regulation.
- Analysis of transcription factor-associated regulatory interactions across tissues.

### 2.6. Integration with complex trait-associated genetic variants

- Integration of tissue-specific CBSs with pig genome-wide association study (GWAS) data.
- Enrichment analysis of trait-associated variants within different classes of CBSs.
- Identification of tissue-specific regulatory regions associated with economically important traits.
- Prioritization of candidate regulatory elements and genes potentially involved in complex trait regulation.

## 3. System environment and software

The analyses were primarily performed in a Linux environment using R, Python, and command-line bioinformatics tools.

Major software packages used in this study include:

- BWA-MEM (v0.7.17)
- GATK (v4.1.4.1)
- hifiasm (v0.25.0-r726)
- dipcall (v0.3)
- WhatsHap (v1.4)
- Juicer (v1.6.2)
- Juicer Tools (v1.22.01)
- deepTools (v3.x)
- BEDTools
- HOMER
- FIMO
- JASPAR
- R (v3.5.1)

Additional software dependencies and analysis parameters are described in the corresponding scripts and manuscript Methods.

## 4. Data availability

The sequencing datasets generated in this study are being deposited in the Genome Sequence Archive (GSA) at the China National Center for Bioinformation (CNCB).

Accession numbers and information about the haplotype-resolved reference genome will be provided here.

## 5. Code availability

This repository provides the custom scripts and computational workflows used for the analyses described in this study.

Scripts will be organized according to their corresponding analytical modules to facilitate reproducibility and reuse.

## 6. Citation

Citation information will be updated upon publication of the associated manuscript.
