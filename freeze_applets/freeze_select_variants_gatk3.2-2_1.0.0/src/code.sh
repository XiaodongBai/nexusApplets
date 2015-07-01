#!/bin/bash

main()
{
    # The following line causes bash to exit at any point if there is any error
    # and to output each line as it is executed -- useful for debugging
    set -e -x -o pipefail

    home_directory=`pwd`

    dx-download-all-inputs
    vcfgzfolder=$(echo $vcfgz_path | rev | cut -d/ -f2- | rev)
    mv $gtbi_path $vcfgzfolder

    # Fetch resources
    PROJECT="$DX_PROJECT_CONTEXT_ID"
    dx download "$PROJECT:/resources/hs37d5.fasta-index.tar.gz"
    tar zxvf hs37d5.fasta-index.tar.gz # => genome.fa, genome.fa.fai, genome.dict
    dx download "$PROJECT:/resources/GenomeAnalysisTK-3.2-2.jar" -o /GenomeAnalysisTK.jar

    # GATK: 90% of machine memory. As many processors as the machine has.
    gatk_megs=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024))}')
    gatk_procs=$(nproc)


    java -Xmx${gatk_megs}m -jar /GenomeAnalysisTK.jar -T SelectVariants -nt ${gatk_procs} -R genome.fa --variant $vcfgz_path -o $vcfgz_prefix.filtered.vcf.gz --exclude_sample_file $samplelist_path

    mkdir -p $HOME/out/selected_vcfgz
    mkdir -p $HOME/out/selected_tbi

    mv $vcfgz_prefix.filtered.vcf.gz $HOME/out/selected_vcfgz
    mv $vcfgz_prefix.filtered.vcf.gz.tbi $HOME/out/selected_tbi

    dx-upload-all-outputs
}
