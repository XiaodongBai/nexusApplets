#!/bin/bash

main()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    # (same for bed/vcf files). This is done because the stat tools use the filename in the report.
    PROJECT="$DX_PROJECT_CONTEXT_ID"
    mkdir bed vcf
    dx download "$PROJECT:/resources/vcr_b37.bed" -o bed/         # => downloads ./bed/<whatever>.bed
    local_bed=(bed/*)                  # => sets $local_bed to "bed/whatever.bed"

    # Similar for vcf file
    dx download "$vcfgz_file" -o vcf/       # => downloads ./vcf/<whatever>.vcf.gz
    local_vcf=(vcf/*)                  # => sets $local_vcf to "vcf/whatever.vcf.gz"
    vcf_file="${local_vcf#vcf/}"     # => sets $vcf_prefix to "whatever.vcf.gz"

    # Modified on 6/20/2014 by adding * so that it is able to handle both vcf and gvcf files
    vcf_prefix="${vcf_file%.*vcf.gz}" # => sets $vcf_prefix to "whatever"

    remote_runs_folder="/data/runs"
    lane="${lane#00}"
    folder="${remote_runs_folder}/${run}/lane.${lane}"

    # But first, define a bash function to help us upload each result into an
    # array of files. Usage: up 'local_filename' 'remote_filename'
    #
    function up {
	id=`dx upload "$1" -o "$2" --brief`
	dx-jobutil-add-output stats "$id" --class=array:file
	echo "$2 $id" >> file_ids.txt
    }

    # BCFtools
    mkdir bcf_out
    dx download "${PROJECT}:${folder}/${vcf_prefix}/${vcf_file}.tbi" -o vcf/
    #tabix -p vcf "$local_vcf"
    tabix -h -B "$local_vcf" "$local_bed" >target.vcf
    bcftools stats target.vcf >vcf.all.stats
    indelline=$(grep "number of indels:" vcf.all.stats)
    bcftools stats --apply-filters PASS target.vcf >vcf.stats
#    novelvarnum=$(bcftools view --apply-filters PASS --novel target.vcf | wc -l)
    novelvarnum=$(bcftools view --novel target.vcf | wc -l)
    sed -i "s/.*number of indels:.*/$indelline/" vcf.stats
    up vcf.stats "$vcf_prefix".stats-vcf.txt
    plot-vcfstats -p bcf_out/ vcf.stats
    up bcf_out/summary.pdf "$vcf_prefix".stats-vcf.pdf

    #
    # Generate JSON file
    #
    gvcfname="${vcf_prefix}.gvcf.gz"
    dx download "${PROJECT}:${folder}/${vcf_prefix}/${gvcfname}"
    dx download "${PROJECT}:${folder}/${vcf_prefix}/${gvcfname}.tbi"
    ygvcf="${vcf_prefix}.Y.gvcf"
    bcftools view -H -t Y -O v  -o "${ygvcf}" "${gvcfname}"
    gender_obs=$(y_basecov.pl "${ygvcf}")

    # Capture metrics stats file is local
    metricsfilename="${vcf_prefix}*.stats-capture-metrics.txt"
    dx download "${PROJECT}:${folder}/${vcf_prefix}*.stats-capture-metrics.txt"

    # Insert size stats file is local
    insertsizefilename="${vcf_prefix}*.stats-insertsize.txt"
    dx download "${PROJECT}:${folder}/${vcf_prefix}*.stats-insertsize.txt"

    # Get startDate and endDate information                                                                                                        
    casavajobfile="run.${run}.lane.${lane}.launched-casava-jobid.txt"
    dxcasavafile=$(dx find data --project "$PROJECT" --folder="$folder" --name "$casavajobfile" --brief)
    jobid=$(dx cat $dxcasavafile)

    creationtime=$(dx describe $dxcasavafile | grep --perl-regex "Created\s\s" | tr -s ' ')

    startmonth=$(echo "${creationtime}" | cut -d' ' -f3)
    startdate=$(echo "${creationtime}" | cut -d' ' -f4)
    starttime=$(echo "${creationtime}" | cut -d' ' -f5)
    startyear=$(echo "${creationtime}" | cut -d' ' -f6)
    startdateAlt="$startdate $startmonth $startyear"
    startdateCorrect=`date -d"$startdateAlt" +%Y/%m/%d`
    startDate="$startdateCorrect $starttime"
    echo "The run processing started: $startDate"

    enddate=`ls --full-time bcf_out/summary.pdf | cut -d' ' -f6`
    endtime=`ls --full-time bcf_out/summary.pdf | cut -d' ' -f7 | cut -d. -f1`
    enddateCorrect=`date -d"$enddate" +%Y/%m/%d`
    endDate="$enddateCorrect $endtime"
    echo "The run processing ended: $endDate"

    stats2json_single_sample.pl $PROJECT $metricsfilename ${insertsizefilename} "vcf.stats" ${run} $flowcell "lane.${lane}" $sample "$startDate" "$endDate" "file_ids.txt" ${novelvarnum} "vcf.all.stats" ${gender_obs}

    time=$(date "+%m%d%y")
    up outputJSON $run.$lane.$sample.stats.v2.0.2.$time.json
}
