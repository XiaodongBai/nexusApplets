#!/bin/bash

main ()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    home_directory=`pwd`
    PROJECT="$DX_PROJECT_CONTEXT_ID"

    dx download "${pvcf}" -o pvcf.vcf.gz

    # Processing
    plink --vcf pvcf.vcf.gz --make-bed --out ${sample} --const-fid Geisinger --threads `nproc`
    cat ${sample}.bim | perl -lane 'if ($F[1] !~ /^rs/) { $F[1] = join(":",($F[0],@F[3..5])); } print join("\t",@F);' > temp
    mv temp ${sample}.bim

    # Output
    mkdir -p $HOME/out/plink_files
    mv ${sample}* $HOME/out/plink_files/
    dx-upload-all-outputs --parallel
}
