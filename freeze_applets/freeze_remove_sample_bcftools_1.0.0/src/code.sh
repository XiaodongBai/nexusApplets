#!/bin/bash

main ()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    home_directory=`pwd`
    PROJECT="$DX_PROJECT_CONTEXT_ID"

    # Input
    dx-download-all-inputs
    vcffolder=$(echo "${vcf_path}" | rev | cut -d/ -f2- | rev)
    mv "${vcf_tbi_path}" "${vcffolder}"

    # Processing
    bcftools view --output-file "$vcf_prefix.selected.vcf.gz" --output-type z --samples-file ^$samplelist_path "${vcf_path}"
    tabix -p vcf "$vcf_prefix.selected.vcf.gz"

    # Output
    mkdir -p $HOME/out/selected_vcfgz
    mkdir -p $HOME/out/selected_vcfgz_tbi
    mv $vcf_prefix.selected.vcf.gz $HOME/out/selected_vcfgz/
    mv $vcf_prefix.selected.vcf.gz.tbi $HOME/out/selected_vcfgz_tbi/
    dx-upload-all-outputs
}
