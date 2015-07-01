#!/bin/bash

main() {
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    # Fetch resources
    PROJECT="$DX_PROJECT_CONTEXT_ID"
    dx download "$PROJECT:/resources/hs37d5.fasta-index.tar.gz"
    tar zxvf hs37d5.fasta-index.tar.gz # => genome.fa, genome.fa.fai, genome.dict
    dx download "$PROJECT:/resources/GenomeAnalysisTK-3.2-2.jar" -o /GenomeAnalysisTK.jar

    dx download "${idfile}" -o idfile

    gvcfs=""
    while read -r gvcffile
    do
	gvcf=$(echo ${gvcffile} | cut -d' ' -f1)
	gvcfs="${gvcfs} -V ${gvcf}"
	dx download "${PROJECT}:/data/GVCF/${gvcf}"
	gunzip ${gvcf}
	gvcfbase="${gvcf%.gz}"
	bgzip ${gvcfbase}
	tabix -p vcf ${gvcf}
    done < <(cat "idfile")


    # GATK: 90% of machine memory. As many processors as the machine has.
    gatk_megs=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024))}')

    # Run all the commands

    java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -R genome.fa -T CombineGVCFs ${gvcfs} -o "${sample}.CombineGVCFs.gvcf.gz"

    # Upload results
    merged_gvcfgz_file_id=`dx upload "${sample}.CombineGVCFs.gvcf.gz" --brief`
    dx-jobutil-add-output "merged_gvcfgz" "$merged_gvcfgz_file_id"
    merged_gvcfgz_tbi_file_id=`dx upload "${sample}.CombineGVCFs.gvcf.gz.tbi" --brief`
    dx-jobutil-add-output "merged_gvcfgz_tbi" "$merged_gvcfgz_tbi_file_id"
}
