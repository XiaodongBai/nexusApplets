#!/bin/bash


main() {

    #
    # Fetch VCF file
    #

    echo "Downloading VCF file..."
    dx download "$vcf_file" -o file.vcf.gz
    vcf_file_name=$(dx describe "$vcf_file" --name)
    prefix="${vcf_file_name%.vcf.gz}"

    #
    # Fetch and unpack annotation resources
    #

    echo "Downloading resources..."
    PROJECT="$DX_PROJECT_CONTEXT_ID"
    dx download "$PROJECT:/resources/rgc_annotation_pipeline_v2.2.0_AmericanBobtail_resources.tar.gz"
    gunzip rgc_annotation_pipeline_v2.2.0_AmericanBobtail_resources.tar.gz
    tar -xvf rgc_annotation_pipeline_v2.2.0_AmericanBobtail_resources.tar
    
    #
    # Run commands
    #
    
    #  /bgzip -cd file.vcf.gz | \
    #  /rgc_annotation_pipeline_version rgc_annotation_pipeline_v2.1.0_AmericanBobtail_sources.txt | \
    #  java -Xmx4g -jar /snpEff-3.6/snpEff.jar eff -noStats -v GRCh37.75 - | /parseSnpEff /Homo_sapiens.GRCh37.75.whiteList.transcripts.txt - | /bgzip -cd | \
    #  java -Xmx4g -jar /snpEff-3.6/snpEff.jar eff -noStats -v GRCh37.75 - | sed 's/=EFF/=eff_ensembl_75/' | sed 's/EFF=/eff_ensembl_75=/g' | \
    #  java -Xmx4g -jar /snpEff-3.6/snpEff.jar eff -noStats -v hg19 - | sed 's/=EFF/=eff_refseq/' | sed 's/EFF=/eff_refseq=/g' | \
    #  java -Xmx4g -jar /snpEff-3.6/snpEff.jar eff -noStats -v hg19kg - | sed 's/=EFF/=eff_knownGenes/' | sed 's/EFF=/eff_knownGenes=/g' | \
    #  /vcfAnnotateBatch /annotationBatch.txt | /bgzip -c > "${prefix}.annotated.vcf.gz"
    #  /tabix -p vcf "${prefix}.annotated.vcf.gz"

    /bgzip -cd file.vcf.gz | \
    /rgc_annotation_pipeline_version rgc_annotation_pipeline_v2.2.0_AmericanBobtail_sources.txt > a.vcf
    java -Xmx4g -jar /snpEff-3.6/snpEff.jar eff -noStats -v GRCh37.75 a.vcf | /parseSnpEff /Homo_sapiens.GRCh37.75.whiteList.transcripts.txt - | /bgzip -cd > b.vcf
    java -Xmx4g -jar /snpEff-3.6/snpEff.jar eff -noStats -v GRCh37.75 b.vcf | sed 's/=EFF/=eff_ensembl_75/' | sed 's/EFF=/eff_ensembl_75=/g' > c.vcf
    java -Xmx4g -jar /snpEff-3.6/snpEff.jar eff -noStats -v hg19 c.vcf | sed 's/=EFF/=eff_refseq/' | sed 's/EFF=/eff_refseq=/g' > d.vcf
    java -Xmx4g -jar /snpEff-3.6/snpEff.jar eff -noStats -v hg19kg d.vcf | sed 's/=EFF/=eff_knownGenes/' | sed 's/EFF=/eff_knownGenes=/g' | \
    /vcfAnnotateBatch /annotationBatch.txt | /bgzip -c > "${prefix}.annotated.vcf.gz"
    /tabix -p vcf "${prefix}.annotated.vcf.gz"

    #
    # Upload files
    #

    file_id=$(dx upload "${prefix}.annotated.vcf.gz" --brief)
    dx-jobutil-add-output vcf_annotated_file "$file_id" --class=file
    file_id=$(dx upload "${prefix}.annotated.vcf.gz.tbi" --brief)
    dx-jobutil-add-output vcf_annotated_file_index "$file_id" --class=file
}

