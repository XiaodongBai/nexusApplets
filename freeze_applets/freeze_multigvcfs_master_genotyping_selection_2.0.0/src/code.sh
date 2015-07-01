#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

# Fetch resources
PROJECT="$DX_PROJECT_CONTEXT_ID"
dx download "$PROJECT:/resources/${master_gvcf}" -o master.gvcf.gz
dx download "$PROJECT:/resources/${master_gvcf_tbi}" -o master.gvcf.gz.tbi

dx download "$PROJECT:/resources/hs37d5.fasta-index.tar.gz"
tar zxvf hs37d5.fasta-index.tar.gz # => genome.fa, genome.fa.fai, genome.dict
dx download "$PROJECT:/resources/GenomeAnalysisTK-3.2-2.jar" -o /GenomeAnalysisTK.jar
dx download "$PROJECT:/resources/gatk_resources-2.8.b37.tar"
tar zxvf gatk_resources-2.8.b37.tar

dx download "$PROJECT:/data/pVCF/${combined_gvcf}" -o gvcf.gz
dx download "$PROJECT:/data/pVCF/${combined_gvcf_tbi}" -o gvcf.gz.tbi

# GATK: 90% of machine memory. As many processors as the machine has.
gatk_megs=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024))}')
gatk_procs=$(nproc)

# CombineGVCFs to combine with master GVCF
java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T GenotypeGVCFs -nt ${gatk_procs} -R genome.fa --variant master.gvcf.gz --variant gvcf.gz -A QualByDepth -A FisherStrand -A DepthPerSampleHC -o ${sample}.master.genotyped.vcf.gz --dbsnp dbsnp_138.b37.vcf.gz -XL MT -XL GL000207.1 -XL GL000226.1 -XL GL000229.1 -XL GL000231.1 -XL GL000210.1 -XL GL000239.1 -XL GL000235.1 -XL GL000201.1 -XL GL000247.1 -XL GL000245.1 -XL GL000197.1 -XL GL000203.1 -XL GL000246.1 -XL GL000249.1 -XL GL000196.1 -XL GL000248.1 -XL GL000244.1 -XL GL000238.1 -XL GL000202.1 -XL GL000234.1 -XL GL000232.1 -XL GL000206.1 -XL GL000240.1 -XL GL000236.1 -XL GL000241.1 -XL GL000243.1 -XL GL000242.1 -XL GL000230.1 -XL GL000237.1 -XL GL000233.1 -XL GL000204.1 -XL GL000198.1 -XL GL000208.1 -XL GL000191.1 -XL GL000227.1 -XL GL000228.1 -XL GL000214.1 -XL GL000221.1 -XL GL000209.1 -XL GL000218.1 -XL GL000220.1 -XL GL000213.1 -XL GL000211.1 -XL GL000199.1 -XL GL000217.1 -XL GL000216.1 -XL GL000215.1 -XL GL000205.1 -XL GL000219.1 -XL GL000224.1 -XL GL000223.1 -XL GL000195.1 -XL GL000212.1 -XL GL000222.1 -XL GL000200.1 -XL GL000193.1 -XL GL000194.1 -XL GL000225.1 -XL GL000192.1 -XL NC_007605 -XL hs37d5

bcftools view --output-file ${sample}.master.genotyped.selected.vcf.gz --output-type z --samples ^template ${sample}.master.genotyped.vcf.gz
tabix -p vcf ${sample}.master.genotyped.selected.vcf.gz

# Upload results
genotyped_vcfgz_file=`dx upload "${sample}.master.genotyped.selected.vcf.gz" --brief`
dx-jobutil-add-output "genotyped_vcfgz" "$genotyped_vcfgz_file"
vcfgz_tbi_file=`dx upload "${sample}.master.genotyped.selected.vcf.gz.tbi" --brief`
dx-jobutil-add-output "genotyped_tbi" "$vcfgz_tbi_file"
