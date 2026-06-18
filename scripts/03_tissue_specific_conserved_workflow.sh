#!/usr/bin/env bash

###############################################################################
# 1. Paths and parameters
###############################################################################
# Tissue-specific and tissue-conserved CTCF binding site analysis.
# Bash sections should be run in bash. Blocks marked "Python script" should be
# run in Python from the project root.
#
# Input:
#   results/cbs/bed/tissue/*.bed
#     One BED file per tissue. The original notes used files ending with
#     _CTCF_peak_sort.bed, but only the BED coordinates are required here.
#
# Output:
#   results/tissue_specific_conserved/specific_peaks/*_CTCF_specific.bed
#   results/tissue_specific_conserved/conserved/shared_peaks_0.5.bed

project_dir="$(pwd)"
tissue_bed_dir="${project_dir}/results/cbs/bed/tissue"
output_dir="${project_dir}/results/tissue_specific_conserved"
threshold=29

mkdir -p "${output_dir}/specific_peaks" "${output_dir}/conserved"

###############################################################################
# 2. Tissue-specific peaks
###############################################################################
# Python script: extract tissue-specific peaks with pybedtools
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

import glob
import os
import pybedtools

project_dir = os.getcwd()
bed_dir = os.path.join(project_dir, "results/cbs/bed/tissue")
output_folder = os.path.join(project_dir, "results/tissue_specific_conserved/specific_peaks")
os.makedirs(output_folder, exist_ok=True)

bed_files = sorted(glob.glob(os.path.join(bed_dir, "*.bed")))

for bed_file in bed_files:
    current = pybedtools.BedTool(bed_file)
    other_files = [f for f in bed_files if f != bed_file]
    other_beds = pybedtools.BedTool.cat(
        *[pybedtools.BedTool(f) for f in other_files],
        postmerge=False
    )

    specific = current.subtract(other_beds)
    tissue = os.path.basename(bed_file).replace("_CTCF_peak_sort.bed", "").replace(".bed", "")
    output_file = os.path.join(output_folder, f"{tissue}_CTCF_specific.bed")
    specific.saveas(output_file)
# End Python script

###############################################################################
# 3. Tissue-conserved peaks present in at least 50% of tissues
###############################################################################
# Python script: extract peaks shared by at least threshold tissue BED files
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

import glob
import os
import pybedtools

project_dir = os.getcwd()
bed_dir = os.path.join(project_dir, "results/cbs/bed/tissue")
output_folder = os.path.join(project_dir, "results/tissue_specific_conserved/conserved")
os.makedirs(output_folder, exist_ok=True)

bed_files = sorted(glob.glob(os.path.join(bed_dir, "*.bed")))
threshold = 29

multi_inter = pybedtools.BedTool().multi_intersect(i=bed_files)
shared_peaks = multi_inter.filter(lambda x: int(x[3]) >= threshold)
shared_peaks.saveas(os.path.join(output_folder, "shared_peaks_0.5.bed"))

print(f"Peaks present in at least {threshold} BED files were written.")
# End Python script
