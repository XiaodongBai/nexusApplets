#!/bin/bash

main() {

    echo "Value of signal_file: '$signal_file'"
    echo "Downloading signal file..."
    dx download "$signal_file" -o signal_file
    fn=$(dx describe "$signal_file" --name)
    echo $fn
    prefix=${fn%.signal.txt}
    echo $prefix

    echo "Downloading and extracting resources..."
    dx download "$DX_PROJECT_CONTEXT_ID:/resources/penncnv_resources.tar.gz" 
    gunzip penncnv_resources.tar.gz
    tar -xvf penncnv_resources.tar

    mkdir -p out/raw_cnv_output out/clean_cnv_output out/log

    perl /detect_cnv.pl -test -hmm hhall.hmm -pfb ghs_500_hg19.pfb -log out/log/${prefix}.log.txt -out out/raw_cnv_output/${prefix}.rawCNVs.txt -confidence -gcmodel HumanOmniExpressExome-8v1-2_B_StrandReport_FP.gcmodel signal_file

    dx-upload-all-outputs --parallel 

}
 