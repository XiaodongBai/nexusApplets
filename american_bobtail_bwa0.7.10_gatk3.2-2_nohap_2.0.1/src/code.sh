#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

home_directory=`pwd`

# Fetch reads
dx download "$reads_fastqgz" -o reads.fq.gz
dx download "$reads2_fastqgz" -o reads2.fq.gz

# Fetch resources
PROJECT="$DX_PROJECT_CONTEXT_ID"
dx download "$PROJECT:/resources/hs37d5.bwa-index.tar.gz"
tar zxvf hs37d5.bwa-index.tar.gz # => genome.bwt, genome.pac, ...
dx download "$PROJECT:/resources/hs37d5.fasta-index.tar.gz"
tar zxvf hs37d5.fasta-index.tar.gz # => genome.fa, genome.fa.fai, genome.dict
dx download "$PROJECT:/resources/gatk_resources-2.8.b37.tar"
tar xvf gatk_resources-2.8.b37.tar # => dbsnp_138.b37.vcf.gz, dbsnp_138.b37.vcf.gz.tbi, ...
dx download "$PROJECT:/resources/GenomeAnalysisTK-3.2-2.jar" -o /GenomeAnalysisTK.jar
dx download "$PROJECT:/resources/vcr_b37_100bpBuff.interval_list" -o intList.interval_list
dx download "$PROJECT:/resources/vcr_b37.bed" -o vcr.bed

# Parameters specific to instance type
# BWA: As many processors as the machine has.
bwa_procs=$(nproc)

# SortSam: 90% of machine memory, after removing 6GB for BWA
# Also calculate maximum records in RAM (rule of thumb taken from
# http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page )
sort_megs=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024 - 6000))}')
sort_recs=$(echo $sort_megs | awk '{print int(250000*$1/1024)}')

# MarkDuplicates: 90% of machine memory.
markdup_megs=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024))}')
markdup_recs=$(echo $markdup_megs | awk '{print int(250000*$1/1024)}')

# GATK: 90% of machine memory. As many processors as the machine has.
gatk_megs=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024))}')
gatk_procs=$(nproc)

# Run all the commands
bwa mem -t $bwa_procs -M -R "@RG\tID:${sample}\tPL:ILLUMINA\tSM:${sample}" genome.fa reads.fq.gz reads2.fq.gz | java -Xmx${sort_megs}m -jar /SortSam.jar MAX_RECORDS_IN_RAM=$sort_recs I=/dev/stdin O=${sample}.bwa.mapped.sorted.bam VALIDATION_STRINGENCY=LENIENT SO=coordinate CREATE_INDEX=true

java -Xmx${markdup_megs}m -jar /MarkDuplicates.jar MAX_RECORDS_IN_RAM=$markdup_recs I=${sample}.bwa.mapped.sorted.bam O=${sample}.bwa.markdup.bam VALIDATION_STRINGENCY=LENIENT AS=T CREATE_INDEX=true M=${sample}.bwa.markdup.log

java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -nt $gatk_procs -T RealignerTargetCreator -R genome.fa -I ${sample}.bwa.markdup.bam -o ${sample}.realign.intervals -known 1000G_phase1.indels.b37.vcf.gz -known Mills_and_1000G_gold_standard.indels.b37.vcf.gz

java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T IndelRealigner -R genome.fa -I ${sample}.bwa.markdup.bam -targetIntervals ${sample}.realign.intervals -known 1000G_phase1.indels.b37.vcf.gz -known Mills_and_1000G_gold_standard.indels.b37.vcf.gz -o ${sample}.realign.bam


# Upload results
bam_file_id=`dx upload "${sample}.realign.bam" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "realigned_bam" "$bam_file_id"  --class file --array
propagate-user-meta "$reads_fastqgz" "$bam_file_id"
propagate-user-meta "$reads2_fastqgz" "$bam_file_id"
bai_file_id=`dx upload "${sample}.realign.bai" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "realigned_bai" "$bai_file_id"  --class file --array
propagate-user-meta "$reads_fastqgz" "$bai_file_id"
propagate-user-meta "$reads2_fastqgz" "$bai_file_id"
