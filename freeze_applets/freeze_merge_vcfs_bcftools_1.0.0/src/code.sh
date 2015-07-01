#!/bin/bash

main ()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    home_directory=`pwd`
    PROJECT="$DX_PROJECT_CONTEXT_ID"

    # Input
    files=$(dx find data --project="${PROJECT}" --folder="/data/pVCF/" --name="*${pattern}*gz" | tr -s ' ' | cut -d' ' -f6)

    printf "%s\n" "${files}" | while read -r file
    do
	localname=$(echo "${file}" | rev | cut -d/ -f1 | rev)
	echo "${localname}" >> file_list.txt
	dx download "${PROJECT}:${file}"
	dx download "${PROJECT}:${file}.tbi"
    done
    echo >> file_list.txt

    outdir="$HOME/out/merged_files"
    mkdir -p $outdir

    # Processing
    bcftools merge --output "$outdir/$sample.1.pVCF.vcf.gz" --output-type z --file-list file_list.txt 1 &
    bcftools merge --output "$outdir/$sample.2.pVCF.vcf.gz" --output-type z --file-list file_list.txt 2 &
    bcftools merge --output "$outdir/$sample.3.pVCF.vcf.gz" --output-type z --file-list file_list.txt 3 &
    bcftools merge --output "$outdir/$sample.4.pVCF.vcf.gz" --output-type z --file-list file_list.txt 4 &
    bcftools merge --output "$outdir/$sample.5.pVCF.vcf.gz" --output-type z --file-list file_list.txt 5 &
    bcftools merge --output "$outdir/$sample.6.pVCF.vcf.gz" --output-type z --file-list file_list.txt 6 &
    bcftools merge --output "$outdir/$sample.7.pVCF.vcf.gz" --output-type z --file-list file_list.txt 7 &
    bcftools merge --output "$outdir/$sample.8.pVCF.vcf.gz" --output-type z --file-list file_list.txt 8 &
    bcftools merge --output "$outdir/$sample.9.pVCF.vcf.gz" --output-type z --file-list file_list.txt 9 &
    bcftools merge --output "$outdir/$sample.10.pVCF.vcf.gz" --output-type z --file-list file_list.txt 10 &
    bcftools merge --output "$outdir/$sample.11.pVCF.vcf.gz" --output-type z --file-list file_list.txt 11 &
    bcftools merge --output "$outdir/$sample.12.pVCF.vcf.gz" --output-type z --file-list file_list.txt 12 &
    bcftools merge --output "$outdir/$sample.13.pVCF.vcf.gz" --output-type z --file-list file_list.txt 13 &
    bcftools merge --output "$outdir/$sample.14.pVCF.vcf.gz" --output-type z --file-list file_list.txt 14 &
    bcftools merge --output "$outdir/$sample.15.pVCF.vcf.gz" --output-type z --file-list file_list.txt 15 &
    bcftools merge --output "$outdir/$sample.16.pVCF.vcf.gz" --output-type z --file-list file_list.txt 16 &
    bcftools merge --output "$outdir/$sample.17.pVCF.vcf.gz" --output-type z --file-list file_list.txt 17 &
    bcftools merge --output "$outdir/$sample.18.pVCF.vcf.gz" --output-type z --file-list file_list.txt 18 &
    bcftools merge --output "$outdir/$sample.19.pVCF.vcf.gz" --output-type z --file-list file_list.txt 19 &
    bcftools merge --output "$outdir/$sample.20.pVCF.vcf.gz" --output-type z --file-list file_list.txt 20 &
    bcftools merge --output "$outdir/$sample.21.pVCF.vcf.gz" --output-type z --file-list file_list.txt 21 &
    bcftools merge --output "$outdir/$sample.22.pVCF.vcf.gz" --output-type z --file-list file_list.txt 22 &
    bcftools merge --output "$outdir/$sample.X.pVCF.vcf.gz" --output-type z --file-list file_list.txt X &
    bcftools merge --output "$outdir/$sample.Y.pVCF.vcf.gz" --output-type z --file-list file_list.txt Y &
    wait

    cd $outdir
    for i in $(seq 1 22) X Y
    do
	echo "$sample.$i.pVCF.vcf.gz" >> concat_list.txt
    done
    bcftools concat --file-list concat_list.txt --output "$sample.pVCF.vcf.gz" --output-type z &
    tabix -p vcf "$sample.1.pVCF.vcf.gz" &
    tabix -p vcf "$sample.2.pVCF.vcf.gz" &
    tabix -p vcf "$sample.3.pVCF.vcf.gz" &
    tabix -p vcf "$sample.4.pVCF.vcf.gz" &
    tabix -p vcf "$sample.5.pVCF.vcf.gz" &
    tabix -p vcf "$sample.6.pVCF.vcf.gz" &
    tabix -p vcf "$sample.7.pVCF.vcf.gz" &
    tabix -p vcf "$sample.8.pVCF.vcf.gz" &
    tabix -p vcf "$sample.9.pVCF.vcf.gz" &
    tabix -p vcf "$sample.10.pVCF.vcf.gz" &
    tabix -p vcf "$sample.11.pVCF.vcf.gz" &
    tabix -p vcf "$sample.12.pVCF.vcf.gz" &
    tabix -p vcf "$sample.13.pVCF.vcf.gz" &
    tabix -p vcf "$sample.14.pVCF.vcf.gz" &
    tabix -p vcf "$sample.15.pVCF.vcf.gz" &
    tabix -p vcf "$sample.16.pVCF.vcf.gz" &
    tabix -p vcf "$sample.17.pVCF.vcf.gz" &
    tabix -p vcf "$sample.18.pVCF.vcf.gz" &
    tabix -p vcf "$sample.19.pVCF.vcf.gz" &
    tabix -p vcf "$sample.20.pVCF.vcf.gz" &
    tabix -p vcf "$sample.21.pVCF.vcf.gz" &
    tabix -p vcf "$sample.22.pVCF.vcf.gz" &
    tabix -p vcf "$sample.X.pVCF.vcf.gz" &
    tabix -p vcf "$sample.Y.pVCF.vcf.gz" &
    wait

    tabix -p vcf "$sample.pVCF.vcf.gz"

    # Output
    dx-upload-all-outputs
}
