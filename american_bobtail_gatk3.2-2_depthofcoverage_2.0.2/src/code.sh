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
dx download "$PROJECT:/resources/GenomeAnalysisTK-3.2-2.jar" -o /GenomeAnalysisTK.jar
dx download "$PROJECT:/resources/vcr_b37_100bpBuff.interval_list" -o intList.interval_list

# GATK: 90% of machine memory. As many processors as the machine has.
gatk_megs=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024))}')
gatk_procs=$(nproc)

# Run all the commands
java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T DepthOfCoverage -R genome.fa -o ${sample}.gatk_readDepth_1x_q30.out -I realigned.bam -L intList.interval_list -mbq 0 -mmq 30 

# Upload results
rd1x_q30_out=`dx upload "${sample}.gatk_readDepth_1x_q30.out" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "rd1x_q30" "$rd1x_q30_out"
propagate-user-meta "$realigned_bam" "$rd1x_q30_out"
propagate-user-meta "$realigned_bai" "$rd1x_q30_out"
rd1x_q30_out_summary=`dx upload "${sample}.gatk_readDepth_1x_q30.out.sample_summary" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "rd1x_q30_summary" "$rd1x_q30_out_summary"
propagate-user-meta "$realigned_bam" "$rd1x_q30_out_summary"
propagate-user-meta "$realigned_bai" "$rd1x_q30_out_summary"
rd1x_q30_out_statistics=`dx upload "${sample}.gatk_readDepth_1x_q30.out.sample_statistics" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "rd1x_q30_statistics" "$rd1x_q30_out_statistics"
propagate-user-meta "$realigned_bam" "$rd1x_q30_out_statistics"
propagate-user-meta "$realigned_bai" "$rd1x_q30_out_statistics"
rd1x_q30_out_interval_summary=`dx upload "${sample}.gatk_readDepth_1x_q30.out.sample_interval_summary" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "rd1x_q30_interval_summary" "$rd1x_q30_out_interval_summary"
propagate-user-meta "$realigned_bam" "$rd1x_q30_out_interval_summary"
propagate-user-meta "$realigned_bai" "$rd1x_q30_out_interval_summary"
rd1x_q30_out_interval_statistics=`dx upload "${sample}.gatk_readDepth_1x_q30.out.sample_interval_statistics" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "rd1x_q30_interval_statistics" "$rd1x_q30_out_interval_statistics"
propagate-user-meta "$realigned_bam" "$rd1x_q30_out_interval_statistics"
propagate-user-meta "$realigned_bai" "$rd1x_q30_out_interval_statistics"
rd1x_q30_out_cumulative_coverage_counts=`dx upload "${sample}.gatk_readDepth_1x_q30.out.sample_cumulative_coverage_counts" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "rd1x_q30_cumulative_coverage_counts" "$rd1x_q30_out_cumulative_coverage_counts"
propagate-user-meta "$realigned_bam" "$rd1x_q30_out_cumulative_coverage_counts"
propagate-user-meta "$realigned_bai" "$rd1x_q30_out_cumulative_coverage_counts"
rd1x_q30_out_cumulative_coverage_proportions=`dx upload "${sample}.gatk_readDepth_1x_q30.out.sample_cumulative_coverage_proportions" --brief --wait -p -o "/$sample/"`
dx-jobutil-add-output "rd1x_q30_cumulative_coverage_proportions" "$rd1x_q30_out_cumulative_coverage_proportions"
propagate-user-meta "$realigned_bam" "$rd1x_q30_out_cumulative_coverage_proportions"
propagate-user-meta "$realigned_bai" "$rd1x_q30_out_cumulative_coverage_proportions"
