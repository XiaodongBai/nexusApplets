#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

home_directory=`pwd`

# Fetch reads
dx download "${eve}" -o evefile.gz

# Fetch resources
PROJECT="$DX_PROJECT_CONTEXT_ID"
dx download "$PROJECT:/resources/hs37d5.fasta-index.tar.gz"
tar zxvf hs37d5.fasta-index.tar.gz # => genome.fa, genome.fa.fai, genome.dict

dx download "$PROJECT:/resources/vcrb37_regions.txt.gz"
gunzip vcrb37_regions.txt.gz

generating_master_gvcf.pl "evefile.gz" "${projectname}" "genome.fa" "vcrb37_regions.txt"

bgzip "${projectname}.master.gvcf"
tabix -p vcf "${projectname}.master.gvcf.gz"

# Upload results
selected_vcfgz_file=`dx upload "${projectname}.master.gvcf.gz" --brief`
dx-jobutil-add-output "master_gvcfgz" "$selected_vcfgz_file"
vcfgz_tbi_file=`dx upload "${projectname}.master.gvcf.gz.tbi" --brief`
dx-jobutil-add-output "master_gtbi" "$vcfgz_tbi_file"
