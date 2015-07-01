#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

home_directory=`pwd`

# Fetch reads
dx download "${vcfgz}" -o ${sample}.vcf.gz
dx download "${tbi}" -o ${sample}.vcf.gz.tbi

# Fetch resources
PROJECT="$DX_PROJECT_CONTEXT_ID"
#dx download "$PROJECT:/resources/vcr_b37_6hm_1bpBuffer.bed" -o target.bed

# Run all the commands

#bcftools view --targets-file target.bed --types indels ${sample}.vcf.gz | bgzip > ${sample}.6hm.indels.vcf.gz
#tabix -p vcf ${sample}.6hm.indels.vcf.gz

clean_pvcf.pl ${sample}.vcf.gz | bgzip > ${sample}.cleaned.vcf.gz
tabix -p vcf ${sample}.cleaned.vcf.gz
gzip -9 filter.log

#bcftools isec --complement --prefix compdir --output-type z ${sample}.cleaned.vcf.gz ${sample}.6hm.indels.vcf.gz

# Upload results
#filtered_vcfgz_file=`dx upload compdir/0000.vcf.gz -o "${sample}.filtered.vcf.gz" --brief`
#dx-jobutil-add-output "filtered_vcfgz" "$filtered_vcfgz_file" --class=file

#vcfgz_tbi_file=`dx upload compdir/0000.vcf.gz.tbi -o "${sample}.filtered.vcf.gz.tbi" --brief`
#dx-jobutil-add-output "filtered_tbi" "$vcfgz_tbi_file" --class=file

filtered_vcfgz_file=`dx upload "${sample}.cleaned.vcf.gz" --brief`
dx-jobutil-add-output "filtered_vcfgz" "$filtered_vcfgz_file" --class=file

vcfgz_tbi_file=`dx upload "${sample}.cleaned.vcf.gz.tbi" --brief`
dx-jobutil-add-output "filtered_tbi" "$vcfgz_tbi_file" --class=file

filter_log_file=`dx upload "filter.log.gz" --brief`
dx-jobutil-add-output "logs" "$filter_log_file" --class=array:file

#indel_vcf_file=`dx upload "${sample}.6hm.indels.vcf.gz" --brief`
#dx-jobutil-add-output "logs" "$indel_vcf_file" --class=array:file

#indel_tbi_file=`dx upload "${sample}.6hm.indels.vcf.gz.tbi" --brief`
#dx-jobutil-add-output "logs" "$indel_tbi_file" --class=array:file
