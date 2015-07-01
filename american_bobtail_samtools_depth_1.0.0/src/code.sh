#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

home_directory=`pwd`

# Fetch reads
dx download "$realigned_bam" -o realigned.bam
dx download "$realigned_bai" -o realigned.bai

# Fetch resources
PROJECT="$DX_PROJECT_CONTEXT_ID"
dx download "$PROJECT:/resources/vcr_b37_100bpBuff.bed" -o targets.bed

# Run all the commands
samtools depth -b targets.bed -q 0 -Q 30 realigned.bam > nonzero_depth.txt

zero_filling.pl targets.bed nonzero_depth.txt > zero_filled_depth.txt

gzip -9 zero_filled_depth.txt

# Upload results
depth_out=`dx upload "zero_filled_depth.txt.gz" --brief --wait -p -o "/${sample}/${sample}.depth.txt.gz"`
dx-jobutil-add-output "files" "$depth_out" --class="array:file"
propagate-user-meta "$realigned_bam" "$depth_out"
propagate-user-meta "$realigned_bai" "$depth_out"
