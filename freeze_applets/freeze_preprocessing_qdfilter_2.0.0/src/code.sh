#!/bin/bash

main ()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    home_directory=`pwd`
    PROJECT="$DX_PROJECT_CONTEXT_ID"

    dx download "${pheader}" -o pheader.vcf.gz
    dx download "${eve}" -o eve.vcf.gz

    # Processing
    vcf_eve_qdfilter pheader.vcf.gz eve.vcf.gz intermediate_file.gz
    qdfilter_info.pl intermediate_file.gz | bgzip > ${sample}.QDfilter_failed.indices.txt.gz

    # Output
    file_id=`dx upload "${sample}.QDfilter_failed.indices.txt.gz" --brief`
    dx-jobutil-add-output "qdfailed_samples" "$file_id"
}
