#!/bin/bash

main()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    PROJECT="$DX_PROJECT_CONTEXT_ID"
    function up {
	id=`dx upload "$1" -o "$2" --brief`
	dx-jobutil-add-output stats "$id" --class=array:file
    }

    tag=$(dx ls ${PROJECT}:/data/runs/${run}/${lane}/${sample}/${sample}*bam | cut -d. -f2)
    filehandle_update_single_sample.pl ${PROJECT} ${run} ${lane} ${sample} ${tag}

    time=$(date "+%m%d%y")
    up outputJSON ${run}.${lane}.${sample}.file_handle.update.$time.json
}
