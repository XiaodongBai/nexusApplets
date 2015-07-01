#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

home_directory=`pwd`

# Fetch reads
dx download "$realigned_bam" -o realigned.bam

# Fetch resources
PROJECT="$DX_PROJECT_CONTEXT_ID"

samtools view -h realigned.bam | perl -e "while (<>) { s/RG:Z:(.+?)\t/RG:Z:$sample\t/; print; }" | samtools view -@ `nproc` -bS -o corrected.bam -
samtools index corrected.bam

# Upload results
bam_id=`dx upload "corrected.bam" --brief --wait -p -o "/$sample/${sample}.corrected.bam"`
dx-jobutil-add-output "corrected_bam" "$bam_id"
propagate-user-meta "$realigned_bam" "$bam_id"
bai_id=`dx upload "corrected.bam.bai" --brief --wait -p -o "/$sample/${sample}.corrected.bai"`
dx-jobutil-add-output "corrected_bai" "$bai_id"
propagate-user-meta "$realigned_bam" "$bai_id"
