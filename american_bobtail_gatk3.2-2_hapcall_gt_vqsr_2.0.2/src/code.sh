#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

home_directory=`pwd`

# Fetch reads
dx download "$realigned_bam" -o realigned.bam
dx download "$realigned_bai" -o realigned.bai

# Fetch resources
PROJECT="$DX_PROJECT_CONTEXT_ID"
dx download "$PROJECT:/resources/hs37d5.fasta-index.tar.gz"
tar zxvf hs37d5.fasta-index.tar.gz # => genome.fa, genome.fa.fai, genome.dict
dx download "$PROJECT:/resources/gatk_resources-2.8.b37.tar"
tar xvf gatk_resources-2.8.b37.tar # => dbsnp_138.b37.vcf.gz, dbsnp_138.b37.vcf.gz.tbi, ...
dx download "$PROJECT:/resources/GenomeAnalysisTK-3.2-2.jar" -o /GenomeAnalysisTK.jar
dx download "$PROJECT:/resources/vcr_b37_100bpBuff.interval_list" -o intList.interval_list
dx download "$PROJECT:/resources/GHS50sample.CombineGVCFs.gvcf.gz" -o combined.gvcf.gz
dx download "$PROJECT:/resources/GHS50sample.CombineGVCFs.gvcf.gz.tbi" -o combined.gvcf.gz.tbi

# GATK: 90% of machine memory. As many processors as the machine has.
gatk_megs=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024))}')
gatk_procs=$(nproc)

# Run all the commands

java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T HaplotypeCaller \
    -R genome.fa \
    -stand_call_conf 30.0 \
    -stand_emit_conf 10.0 \
    --dbsnp dbsnp_138.b37.vcf.gz \
    -o ${sample}.gvcf.gz \
    -I realigned.bam \
    -L intList.interval_list \
    -ERC GVCF \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    -pairHMM VECTOR_LOGLESS_CACHING 

java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T GenotypeGVCFs -nt $gatk_procs \
    -R genome.fa \
    --variant ${sample}.gvcf.gz \
    --variant combined.gvcf.gz \
    -A QualByDepth \
    -A FisherStrand \
    -A DepthPerSampleHC \
    -o ${sample}.annGT.vcf.gz \
    --dbsnp dbsnp_138.b37.vcf.gz 

java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T SelectVariants -R genome.fa --variant ${sample}.annGT.vcf.gz -o ${sample}.select.vcf.gz -se "${sample}" -env

java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T VariantRecalibrator -nt $gatk_procs \
    -R genome.fa \
    -input ${sample}.select.vcf.gz \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf.gz \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 1000G_omni2.5.b37.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz \
    -an DP \
    -an QD \
    -an FS \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile snp.output.recal \
    -tranchesFile snp.output.tranches \
    -rscriptFile snp.output.R \
    -XL NC_007605 -XL hs37d5

java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T ApplyRecalibration -nt $gatk_procs \
    -R genome.fa \
    -input ${sample}.select.vcf.gz \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile snp.output.recal \
    -tranchesFile snp.output.tranches \
    -o ${sample}.vcf.gz \
    -XL NC_007605 -XL hs37d5
#    -o ${sample}.recalibrated_snps_raw_indels.vcf.gz


#java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T VariantRecalibrator -nt $gatk_procs \
#    -R genome.fa \
#    -input ${sample}.recalibrated_snps_raw_indels.vcf.gz \
#    -resource:mills,known=true,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
#    -an DP \
#    -an FS \
#    -an MQRankSum \
#    -an ReadPosRankSum \
#    -mode INDEL \
#    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
#    --maxGaussians 4 \
#    -recalFile indel.output.recal \
#    -tranchesFile indel.output.tranches \
#    -rscriptFile indel.output.R \
#    -XL NC_007605 -XL hs37d5

#java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T ApplyRecalibration -nt $gatk_procs \
#    -R genome.fa \
#    -input ${sample}.recalibrated_snps_raw_indels.vcf.gz \
#    -mode INDEL \
#    --ts_filter_level 99.0 \
#    -recalFile indel.output.recal \
#    -tranchesFile indel.output.tranches \
#    -o ${sample}.vcf.gz \
#    -XL NC_007605 -XL hs37d5

# Upload results
gvcfgz_file_id=`dx upload "${sample}.gvcf.gz" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "gvcfgz" "$gvcfgz_file_id"
propagate-user-meta "$realigned_bam" "$gvcfgz_file_id"
propagate-user-meta "$realigned_bai" "$gvcfgz_file_id"
gtbi_file_id=`dx upload "${sample}.gvcf.gz.tbi" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "gtbi" "$gtbi_file_id"
propagate-user-meta "$realigned_bam" "$gtbi_file_id"
propagate-user-meta "$realigned_bai" "$gtbi_file_id"
genotyped_vcfgz_file=`dx upload "${sample}.vcf.gz" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "genotyped_vcfgz" "$genotyped_vcfgz_file"
propagate-user-meta "$realigned_bam" "$genotyped_vcfgz_file"
propagate-user-meta "$realigned_bai" "$genotyped_vcfgz_file"
vcfgz_tbi_file=`dx upload "${sample}.vcf.gz.tbi" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "genotyped_tbi" "$vcfgz_tbi_file"
propagate-user-meta "$realigned_bam" "$vcfgz_tbi_file"
propagate-user-meta "$realigned_bai" "$vcfgz_tbi_file"
