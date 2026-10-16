# Demo results and quality control

This directory contains selected analysis results and quality-control files for the pig CTCF CUT&Tag workflow in [`scripts/01_cuttag_and_qc_workflow.sh`](../scripts/01_cuttag_and_qc_workflow.sh).

## Sample and reference genome

- **Data source:** [Genome Sequence Archive (GSA)](https://ngdc.cncb.ac.cn/gsa/s/A9393C5E)
- **Run accession:** CRR3293403
- **Sample:** Brain_3
- **Local file prefix:** `HS4M08b1W_CTCF`
- **Reference genome:** Sus scrofa, Sscrofa11.1

## Included files

| File | Description |
| --- | --- |
| `HS4M08b1W_CTCF_bowtie2.txt` | Bowtie2 alignment summary. |
| `HS4M08b1W_CTCF_fragmentLen.txt` | Fragment-length distribution before duplicate removal; columns are fragment length (bp) and count. |
| `HS4M08b1W_CTCF.pbc.qc` | Library complexity metrics; columns are total fragments, distinct fragment positions, positions observed once, positions observed twice, NRF, PBC1, and PBC2. |
| `HS4M08b1W_CTCF_FRIP_score.txt` | Number of reads overlapping peaks, total reads, and FRiP percentage. |
| `HS4M08b1W_CTCF.qual` | Phantompeakqualtools cross-correlation QC output, including NSC and RSC. |
| `HS4M08b1W_CTCF.Rout` | Phantompeakqualtools execution log. |
| `HS4M08b1W_CTCF_peaks.narrowPeak` | MACS2 narrowPeak output. |
| `HS4M08b1W_CTCF_peaks.xls` | MACS2 peak table. |
| `qc_summary.tsv` | Summary of FRiP and library complexity metrics. |

This directory contains a subset of the workflow outputs. Raw FASTQ files, alignment files, and other intermediate files are not included.
