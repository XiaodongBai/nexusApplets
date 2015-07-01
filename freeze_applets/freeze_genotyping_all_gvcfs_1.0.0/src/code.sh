#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

# Fetch resources
PROJECT="$DX_PROJECT_CONTEXT_ID"
dx download "$PROJECT:/resources/hs37d5.fasta-index.tar.gz"
tar zxvf hs37d5.fasta-index.tar.gz # => genome.fa, genome.fa.fai, genome.dict
dx download "$PROJECT:/resources/GenomeAnalysisTK-3.2-2.jar" -o /GenomeAnalysisTK.jar

gvcffiles=$(dx find data --project="${PROJECT}" --folder="/data/pVCF/" --name="*.CombineGVCFs.gvcf.gz" | tr -s ' ' | cut -d' ' -f6)

gvcfs=""
while read -r gvcffile
do
    gvcf=$(echo ${gvcffile} | rev | cut -d/ -f1 | rev)
    gvcfs="${gvcfs} -V ${gvcf}"
    dx download "${PROJECT}:/data/pVCF/${gvcf}"
    dx download "${PROJECT}:/data/pVCF/${gvcf}.tbi"
done < <(printf '%s\n' "${gvcffiles}")

# GATK: 90% of machine memory. As many processors as the machine has.
gatk_megs=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024))}')
gatk_procs=$(nproc)

# Run all the commands

java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T GenotypeGVCFs -nt ${gatk_procs} -R genome.fa ${gvcfs} -A QualByDepth -A FisherStrand -A DepthPerSampleHC -o ${sample}.pVCF.vcf.gz

# Upload results
genotyped_gvcfgz_file_id=`dx upload "${sample}.pVCF.vcf.gz" --brief`
dx-jobutil-add-output "genotyped_gvcfgz" "$genotyped_gvcfgz_file_id"
genotyped_gvcfgz_tbi_file_id=`dx upload "${sample}.pVCF.vcf.gz.tbi" --brief`
dx-jobutil-add-output "genotyped_gvcfgz_tbi" "$genotyped_gvcfgz_tbi_file_id"
