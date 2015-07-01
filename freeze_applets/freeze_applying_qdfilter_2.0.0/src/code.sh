#!/bin/bash

main ()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    home_directory=`pwd`
    PROJECT="$DX_PROJECT_CONTEXT_ID"

    dx download "${pvcf}" -o pvcf.vcf.gz
    dx download "${ptbi}" -o pvcf.vcf.gz.tbi
    dx download "${info}" -o info.file.gz

    # Processing
    pvcf_qdfilter pvcf.vcf.gz info.file.gz ${sample}.QDfiltered.pVCF.vcf.gz
    tabix -p vcf "${sample}.QDfiltered.pVCF.vcf.gz"

    # Output
    filtered_vcfgz_file_id=`dx upload "${sample}.QDfiltered.pVCF.vcf.gz" --brief`
    dx-jobutil-add-output "filtered_vcfgz" "$filtered_vcfgz_file_id"
    filtered_vcfgz_tbi_file_id=`dx upload "${sample}.QDfiltered.pVCF.vcf.gz.tbi" --brief`
    dx-jobutil-add-output "filtered_tbi" "$filtered_vcfgz_tbi_file_id"
}
