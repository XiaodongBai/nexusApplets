#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
#
# Fetch files
#
#

# Download bam file inside a directory called "bam" while keeping the original filename
# (same for bed/vcf files). This is done because the stat tools use the filename in the report.
PROJECT="$DX_PROJECT_CONTEXT_ID"
mkdir bam bed vcf
dx download "$bam" -o ./bam/       # => downloads ./bam/<whatever>.bam
local_bam=(bam/*)                  # => sets $local_bam to "bam/whatever.bam"
bam_prefix="${local_bam#bam/}"     # => sets $bam_prefix to "whatever.bam"
bam_prefix="${bam_prefix%.bam}"    # => sets $bam_prefix to "whatever"

# Similar for the bed file
dx download "$PROJECT:/resources/vcr_b37.bed" -o bed/         # => downloads ./bed/<whatever>.bed
local_bed=(bed/*)                  # => sets $local_bed to "bed/whatever.bed"

dx download "${PROJECT}:/resources/hs37d5.fasta-index.tar.gz"
tar zxvf hs37d5.fasta-index.tar.gz # => genome.fa, genome.fa.fai, genome.dict

#
# Calculate 90% of memory size, for java
#
mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.9/1024)}'`
java="java -Xmx${mem_in_mb}m"

#
# Run all the commands and upload results
#

# But first, define a bash function to help us upload each result into an
# array of files. Usage: up 'local_filename' 'remote_filename'
#
function up {
  id=`dx upload "$1" -o "$2" --brief`
  dx-jobutil-add-output stats "$id" --class=array:file
#  echo "$2 $id" >> file_ids.txt
}

# Picard CollectAlignmentSummaryMetrics
$java -jar /picard/CollectAlignmentSummaryMetrics.jar I="$local_bam" O=summarymetrics.txt R=genome.fa VALIDATION_STRINGENCY=SILENT
up summarymetrics.txt "$bam_prefix".stats-alignment-summary.txt

# Picard CollectInsertSizeMetrics
$java -jar /picard/CollectInsertSizeMetrics.jar H=insert.pdf I="$local_bam" O=paired.txt R=genome.fa VALIDATION_STRINGENCY=SILENT MINIMUM_PCT=0.5
up insert.pdf "$bam_prefix".stats-insertsize.pdf
up paired.txt "$bam_prefix".stats-insertsize.txt

# Picard CalculateHsMetrics -- unfortunately doesn't work with bed files directly,
# so first we have to convert the bed into a "plist"
mkdir -p plist/bed
local_plist="plist/$local_bed"
samtools view -H "$local_bam" | grep '^@SQ' > "$local_plist"
awk '{print $1 "\t" $2+1 "\t" $3 "\t+\t" $1 ":" $2+1 "-" $3}' < "$local_bed" >> "$local_plist"

$java -jar /picard/CalculateHsMetrics.jar TI="$local_plist" BI="$local_plist" I="$local_bam" PER_TARGET_COVERAGE=pertarget.txt O=hsmetrics.txt R=genome.fa VALIDATION_STRINGENCY=SILENT
up pertarget.txt "$bam_prefix".stats-pertarget-coverage.txt
up hsmetrics.txt "$bam_prefix".stats-capture-metrics.txt

# FastQC
# Modified on 6/20/2014 to use fastqc version 0.11.2, which doesn't have a fastqc jar file
mkdir fastqc_out
/fastqc/fastqc -o fastqc_out -t `nproc` -f bam --nogroup --extract --contaminants /fastqc/Configuration/contaminant_list.txt $local_bam
rm -f fastqc_out/*zip
up fastqc_out/*/fastqc_data.txt "$bam_prefix".stats-fastqc.txt
dx-build-report-html --local stats-fastqc.html fastqc_out/*/*html
up stats-fastqc.html "$bam_prefix".stats-fastqc.html
