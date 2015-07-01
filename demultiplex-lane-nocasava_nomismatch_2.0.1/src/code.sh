#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

PROJECT="$DX_PROJECT_CONTEXT_ID"

# Get startDate information
startDate=`date +"%x %T"`

mkdir -p input output

# Download the sample sheet
dx download "$sample_sheet" -o sample_sheet.csv

# Remove extras from sample sheet for CASAVA
cut -f1-10 -d, < sample_sheet.csv > casava_sample_sheet.csv

# Download and unpack run data file(s)
for i in "${!run_data[@]}"
do
  dx download "${run_data[$i]}" -o run-data.tar
  tar xf run-data.tar -C ./input/
  rm -f run-data.tar
done

# Locate the "BaseCalls" directory or fail otherwise.
input_dir=$(find input -type d -name BaseCalls)
if [ "$input_dir" == "" ]; then
  # Couldn't find a BaseCalls directory. Fail.
  dx-jobutil-report-error "The BaseCalls/ directory could not be found in the run directory tarball." AppError
fi

#
# Run CASAVA
# ----------
#
configureBclToFastq.pl --input-dir "$input_dir" --output-dir output/Unaligned --sample-sheet casava_sample_sheet.csv --fastq-cluster-count 0 $advanced_opts
make -C output/Unaligned -j `nproc`

#
# Additional metadata extraction
# ------------------------------
#
# Extract instrument, flowcell, and run folder from the Bustard XML config
#
instrument=$(xmllint output/Unaligned/Basecall_Stats_*/config.xml --xpath '/BaseCallAnalysis/Run/RunParameters/Instrument/text()')
run_folder=$(xmllint output/Unaligned/Basecall_Stats_*/config.xml --xpath '/BaseCallAnalysis/Run/RunParameters/RunFolder/text()')
flowcell=$(grep '^FLOWCELL:=' output/Unaligned/Makefile | cut -f2 -d=)

#
# Upload results, add useful meta
# -------------------------------
#

declare -A lanes

i=0
for project_dir in output/Unaligned/Project_*
do
  project=$(echo "$project_dir" | cut -f3 -d/ | cut -f2- -d_)
  for sample_dir in "$project_dir"/Sample_*
  do
    let i=i+1
    sample=$(echo "$sample_dir" | cut -f4 -d/ | cut -f2- -d_)
    for file in "$sample_dir"/*.fastq.gz
    do
      #
      # Typical file: output/Unaligned/Project_<project>/Sample_<sample>/<sample>_<index>_L<lane>_R<mate>_<chunk>.fastq.gz 
      #
      index=$(echo "$file" | rev | cut -f4 -d_ | rev)
      lane=$(echo "$file" | rev | cut -f3 -d_ | rev | cut -c2-)
      mate=$(echo "$file" | rev | cut -f2 -d_ | rev | cut -c2-)
      chunk=$(echo "$file" | rev | cut -f1 -d_ | rev | cut -f1 -d.)
      clean_lane=$(echo "$lane" | sed 's/^0*//')
      lanes[$clean_lane]=true

      # Extract contributor from expanded sample sheet
      contributor=$(awk -F, '{if ($3 == "'$sample'") print $11}' < sample_sheet.csv)

      file_id=$(dx upload $file -p --brief -o "/$sample/" \
        --property instrument="$instrument" \
        --property flowcell="$flowcell" \
        --property run_folder="$run_folder" \
        --property project="$project" \
        --property sample="$sample" \
        --property index="$index" \
        --property lane="$lane" \
        --property mate="$mate" \
        --property chunk="$chunk" \
        --property contributor="$contributor" )
      if [[ "$mate" == "1" ]]; then
	  read1="-ireads_fastqgz=$file_id"
        dx-jobutil-add-output reads "$file_id" --class file --array
      else
	  read2="-ireads2_fastqgz=$file_id"
        dx-jobutil-add-output reads2 "$file_id" --class file --array
      fi
    done
    cmd="dx run rgc_pipeline_bwa_gatk3.1_nohap_1.0.6 $read1 $read2 -isample=$sample -y --brief"
    echo $cmd >> concat_cmd
  done
done

which_lanes=${!lanes[@]}
which_lanes=${which_lanes// /,}

file_id=$(dx upload output/Unaligned/Basecall_Stats_*/BustardSummary.xml --brief \
  --property instrument="$instrument" \
  --property flowcell="$flowcell" \
  --property run_folder="$run_folder" \
  -o "run.$run_folder.lane.$which_lanes.BustardSummary.xml")
dx-jobutil-add-output stats "$file_id" --class file --array
echo "run.$run_folder.lane.$which_lanes.BustardSummary.xml $file_id" >> file_ids.txt

file_id=$(dx upload output/Unaligned/Basecall_Stats_*/Flowcell_demux_summary.xml --brief \
  --property instrument="$instrument" \
  --property flowcell="$flowcell" \
  --property run_folder="$run_folder" \
  -o "run.$run_folder.lane.$which_lanes.Flowcell_demux_summary.xml")
dx-jobutil-add-output stats "$file_id" --class file --array
echo "run.$run_folder.lane.$which_lanes.Flowcell_demux_summary.xml $file_id" >> file_ids.txt

file_id=$(dx upload output/Unaligned/Basecall_Stats_*/Demultiplex_Stats.htm --brief \
  --property instrument="$instrument" \
  --property flowcell="$flowcell" \
  --property run_folder="$run_folder" \
  -o "run.$run_folder.lane.$which_lanes.Demultiplex_Stats.htm")
dx-jobutil-add-output stats "$file_id" --class file --array
echo "run.$run_folder.lane.$which_lanes.Demultiplex_Stats.htm $file_id" >> file_ids.txt

file_id=$(dx upload concat_cmd --brief \
  -o "run.$run_folder.lane.$which_lanes.bwa.commands.txt")
dx-jobutil-add-output stats "$file_id" --class file --array

# Enddata is defined as the time when all files have been uploaded
endData=`date +"%x %T"`

TMPPROJECT=${DX_WORKSPACE_ID}
mv output/Unaligned/Basecall_Stats_*/Demultiplex_Stats.htm Demultiplex_Stats.htm
stats2json_single_lane.pl $TMPPROJECT "Demultiplex_Stats.htm" $run_folder $flowcell $which_lanes "$startDate" "$endData" "file_ids.txt"

time=$(date "+%m%d%y")
file_id=$(dx upload outputJSON -o $run_folder.$which_lanes.stats.v2.0.0.$time.json --brief)
dx-jobutil-add-output stats "$file_id" --class file --array
