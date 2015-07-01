#!/bin/bash

main ()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    home_directory=`pwd`
    PROJECT="$DX_PROJECT_CONTEXT_ID"

    dx-download-all-inputs --parallel

    pvcfnum=${#pvcfs[@]}
    let pvcfnum=pvcfnum-1
    for index in $(seq 0 ${pvcfnum})
    do
	echo "${pvcfs_path[$index]}" >> file_list.txt
	folder=$(echo "${pvcfs_path[$index]}" | rev | cut -d/ -f2- | rev)
	mv "${ptbis_path[$index]}" "${folder}"
    done

    # processing
    bcftools concat --output "${sample}.QDfiltered.pVCF.vcf.gz" --output-type z --file-list file_list.txt
    tabix -p vcf "${sample}.QDfiltered.pVCF.vcf.gz"

    mkdir -p $HOME/out/merged_vcfgz
    mkdir -p $HOME/out/merged_vcfgz_tbi
    mv "${sample}.QDfiltered.pVCF.vcf.gz" $HOME/out/merged_vcfgz
    mv "${sample}.QDfiltered.pVCF.vcf.gz.tbi" $HOME/out/merged_vcfgz_tbi
    dx-upload-all-outputs
}
