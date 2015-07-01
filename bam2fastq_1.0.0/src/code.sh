#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

home_directory=`pwd`

# Fetch reads
dx download "$bam" -o ${sample}.bam

# Fetch resources
PROJECT="$DX_PROJECT_CONTEXT_ID"

# Run all the commands
samtools sort -@ `nproc` -n "${sample}.bam" "${sample}.sortByName"
bam2fastq.pl "${sample}.sortByName.bam" "${sample}.regenerated"
gzip "${sample}.regenerated.R1.fq"
gzip "${sample}.regenerated.R2.fq"

# Upload results
reads1_file_id=`dx upload "${sample}.regenerated.R1.fq.gz" --brief`
dx-jobutil-add-output "reads1_fastqgz" "$reads1_file_id"  --class file
reads2_file_id=`dx upload "${sample}.regenerated.R2.fq.gz" --brief`
dx-jobutil-add-output "reads2_fastqgz" "$reads2_file_id"  --class file
