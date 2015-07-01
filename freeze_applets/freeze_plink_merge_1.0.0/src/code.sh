#!/bin/bash

main ()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    home_directory=`pwd`
    PROJECT="$DX_PROJECT_CONTEXT_ID"

    for chr in $(seq 1 22) X Y
    do
	dx download $(dx find data --project="$PROJECT" --folder="/data/" --name="${sample}.${chr}.${pattern}.bed" --brief)
	dx download $(dx find data --project="$PROJECT" --folder="/data/" --name="${sample}.${chr}.${pattern}.bim" --brief)
	dx download $(dx find data --project="$PROJECT" --folder="/data/" --name="${sample}.${chr}.${pattern}.fam" --brief)
    done
    
    for chr in $(seq 2 22) X Y
    do
	echo "${sample}.${chr}.${pattern}.bed ${sample}.${chr}.${pattern}.bim ${sample}.${chr}.${pattern}.fam" >> plink_merge_list.txt
    done

    # Processing
    plink --merge-list plink_merge_list.txt --make-bed --out ${sample}.${pattern} --bfile ${sample}.1.${pattern}


    # Output
    mkdir -p $HOME/out/plink_files
    mv ${sample}.${pattern}.bed ${sample}.${pattern}.bim ${sample}.${pattern}.fam $HOME/out/plink_files/
    dx-upload-all-outputs --parallel
}
