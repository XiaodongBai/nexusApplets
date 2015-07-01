#!/bin/bash



main() 
{
    echo "Value of dx_project: '$dx_project'"
    echo "Value of sample_id: '$sample_id'"
    
    echo "Downloading coverage data..."
    file_id=$(dx find data --project "$dx_project" --folder /data/runs --name "${sample_id}.gatk_readDepth_1x_q30.out" --brief)
    echo "File ID of coverage file: $file_id"
    dx download "$file_id"

    echo "Converting coverage data to BED file..."
    /gatk_coverage_to_bed "${sample_id}.gatk_readDepth_1x_q30.out" /windows.bed > "${sample_id}.coverage.bed"
    echo "Normalizing coverage data..."
    /normalize_coverage "${sample_id}.coverage.bed" /windows.bed > "${sample_id}.coverage.normalized.bed"
    SEX=`grep "^Y" "${sample_id}.coverage.normalized.bed" | awk '$5 == 0 { x += $4; n++ } END { if (x/n > 0.1) print "M"; else print "F" }'`
     echo "Calling CNVs..."
    /call_cnv "${sample_id}.coverage.normalized.bed" /models.out --sex "$SEX" > "${sample_id}.cnv.bed"
    
    echo "Uploading results..."
    file_id=$(dx upload "${sample_id}.cnv.bed" --brief)
    dx-jobutil-add-output cnv_bed "$file_id" --class=file

    file_id=$(dx upload "${sample_id}.coverage.bed" --brief)
    dx-jobutil-add-output coverage_bed "$file_id" --class=file
}
