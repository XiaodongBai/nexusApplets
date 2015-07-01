#$ -S /bin/bash
#$ -l h_vmem=20G -l mem_requested=20G
#$ -e /net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/cluster_output.e
#$ -q nick-grad.q
#$ -o /net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/cluster_output.o

cd /nfs/home/grapas2/projects/2011/reconstruct_pedigrees/PRIMUS_dev/bin/
module load R/3.0.0
#pwd
GSITSHOST=`/bin/hostname`
GSITSPWD=`/bin/pwd` 
GSITSDATE=`/bin/date` 

THRESHOLD='0.375'
THRESHOLD='0.1875'
#THRESHOLD='0.1680'
#THRESHOLD='0.09375'
#THRESHOLD='0.11'

MITO='0.1'

#module unload perl
#perl -v

#perl ../PRIMUS_dev/bin/run_PRIMUS.pl -p ../data/7_9_14_hhu3/full.genome -t $THRESHOLD -o ../data/7_9_14_hhu3/full.genome_$THRESHOLD --cluster

perl run_PRIMUS.pl -p ../../data/JHS/mergedJHS_ARIC_PRIMUS/mergedJHS_ARIC_prePRIMUS/mergedJHS_ARIC_cleaned.genome -t $THRESHOLD -o ../../data/JHS/mergedJHS_ARIC_PRIMUS_$THRESHOLD --cluster --sexes FILE=../../data/JHS/mergedJHS_ARIC_PRIMUS/mergedJHS_ARIC_prePRIMUS/mergedJHS_ARIC_cleaned.sexcheck SEX=4 --age_file ../../data/JHS/JHS_age.txt --affection_file ../../data/JHS/mergedJHS_ARIC.affections.txt

perl run_PRIMUS.pl -p ../../data/JHS/mergedJHS_ARIC_PRIMUS/mergedJHS_ARIC_prePRIMUS/mergedJHS_ARIC_cleaned.genome -t $THRESHOLD -o ../../data/JHS/mergedJHS_ARIC_PRIMUS_$THRESHOLD\_mito$MITO --cluster --sexes FILE=../../data/JHS/mergedJHS_ARIC_PRIMUS/mergedJHS_ARIC_prePRIMUS/mergedJHS_ARIC_cleaned.sexcheck SEX=4 --age_file ../../data/JHS/JHS_age.txt --affection_file ../../data/JHS/mergedJHS_ARIC.affections.txt --mito_matches FILE=../../data/JHS/mergedJHS_ARIC_PRIMUS/mergedJHS_ARIC_prePRIMUS/mergedJHS_ARIC_MT_estimates.txt --y_matches FILE=../../data/JHS/mergedJHS_ARIC_PRIMUS/mergedJHS_ARIC_prePRIMUS/mergedJHS_ARIC_Y_estimates.txt --MT_error_rate $MITO



exit;

#perl run_PRIMUS.pl -p ../data/Silverman/silverman_uwcmg_eocopd_1.batch_1.all.polymorphic.filtered_PRIMUS/silverman_uwcmg_eocopd_1.batch_1.all.polymorphic.filtered_PLINK_IBD/results/silverman_uwcmg_eocopd_1.batch_1.all.polymorphic.filtered.genome --sexes FILE=../data/Silverman/silverman_uwcmg_eocopd_1.batch_1.all.polymorphic.filtered_PRIMUS/silverman_uwcmg_eocopd_1.batch_1.all.polymorphic.filtered_PLINK_IBD/Sex/silverman_uwcmg_eocopd_1.batch_1.all.polymorphic.filtered.sexcheck SEX=4 --output_dir ../data/Silverman/silverman_uwcmg_eocopd_1.batch_1.all.polymorphic.filtered_PRIMUS_cluster_noisy_0.102 --cluster --threshold 0.102


## Santhi Pedigrees
#perl run_PRIMUS.pl -p ../data/Santhi/sample_filter_IBD_geno.poly.indep_1000000_1000_0.2.prune.in.genome_no_complex_or_distant_peds.genome -o ../data/Santhi2/ --sexes FILE=/net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/data/Santhi/select_sample.fam SEX=5 --cluster --threshold $THRESHOLD


## PC_pedigrees
#perl run_PRIMUS.pl -p ../data/Stanford/Stanford_Clean_Final_IID_names_no_complex_or_distant_rels_PRIMUS/Stanford_Clean_Final_IID_names_no_complex_or_distant_rels_PLINK_IBD/results/Stanford_Clean_Final_IID_names_no_complex_or_distant_rels.genome -o ../data/Stanford/Stanford_Clean_Final_IID_names_no_complex_or_distant_rels_PRIMUS_$THRESHOLD --sexes FILE=/net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/data/PCfamilies.ped SEX=5 --affections FILE=/net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/data/PCfamilies.ped AFFECTION=6 --cluster --threshold $THRESHOLD

#perl run_PRIMUS.pl -p /net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/data/PC_pedigree_2.genome -o ../data/Stanford/ --sexes FILE=../data/PCfamilies.ped SEX=5 --affections FILE=../data/PCfamilies.ped AFFECTION=6 --cluster --threshold $THRESHOLD

#./run_dataset_cluster.pl PC_pedigrees_v1_$THRESHOLD ../data/PCfirstpass.genome ../data/PCfamilies.ped $THRESHOLD 4 5


## SJK
#perl run_PRIMUS.pl -p ../data/SJK/all_sjk_final_PRIMUS/all_sjk_final_PLINK_IBD/results/all_sjk_final.genome --cluster --threshold 0.375 --sexes FILE=../data/SJK/all_sjk_final.sexcheck.edited SEX=3 --affections FILE=../data/SJK/all_sjk_final.pedind_edited AFFECTION=6 AFFECTION_VALUE=2 --output ../data/SJK/all_sjk_final_PRIMUS5


## Najma
#perl run_PRIMUS.pl -p ../data/GENNID_no_GP/GENNID_all10012012_noGP_real_samples_PRIMUS/GENNID_all10012012_noGP_real_samples_PLINK_IBD/results/GENNID_all10012012_noGP_real_samples.genome --cluster --threshold 0.09375 --sexes FILE=../data/GENNID/GENNID_all10012012_real_samples_PLINK_IBD/Sex/GENNID_all10012012_real_samples.sexcheck SEX=4 --output ../data/GENNID_no_GP/GENNID_all10012012_noGP_real_samples_PRIMUS5/ --ages FILE=../data/GENNID/pedin-GENNID.txt 

#### GENNID
#perl run_PRIMUS.pl -p ../data/GENNID_no_GP/GENNID_all10012012_noGP_real_samples_PRIMUS/GENNID_all10012012_noGP_real_samples_PLINK_IBD/results/GENNID_all10012012_noGP_real_samples.genome --cluster --threshold 0.09375 --sexes FILE=../data/GENNID/GENNID_all10012012_real_samples_PLINK_IBD/Sex/GENNID_all10012012_real_samples.sexcheck SEX=4 --output ../data/GENNID_no_GP/GENNID_all10012012_noGP_real_samples_PRIMUS5/ --ages FILE=../data/GENNID/pedin-GENNID.txt 

## no AIMs
#perl run_PRIMUS.pl -p ../data/GENNID/AIMs_removed_1/GENNID_all10012012_noGP_real_samples_PLINK_IBD/results/GENNID_all10012012_noGP_real_samples.genome --cluster --threshold 0.09375 --sexes FILE=../data/GENNID/GENNID_all10012012_real_samples_PLINK_IBD/Sex/GENNID_all10012012_real_samples.sexcheck SEX=4 --output ../data/GENNID/AIMs_removed_3 --affections FILE=../data/GENNID/pedin-GENNID.txt AFFECTION=6 --ages FILE=../data/GENNID/pedin-GENNID.txt AGE=8

## with AIMs
#perl run_PRIMUS.pl -p ../data/GENNID/AIMs_not_removed_1/GENNID_all10012012_noGP_real_samples_PLINK_IBD/results/GENNID_all10012012_noGP_real_samples.genome --cluster --threshold 0.09375 --sexes FILE=../data/GENNID/GENNID_all10012012_real_samples_PLINK_IBD/Sex/GENNID_all10012012_real_samples.sexcheck SEX=4 --output ../data/GENNID/AIMs_not_removed_3 --affections FILE=../data/GENNID/pedin-GENNID.txt AFFECTION=6 --ages FILE=../data/GENNID/pedin-GENNID.txt AGE=8

#### MEXAM
## no AIMs
#./run_PRIMUS.pl --cluster --plink ../data/mexam/AIMs_removed_1/justStarrCo_CC_PLINK_IBD/results/justStarrCo_CC.genome --sexes FILE=../data/mexam/justStarrCo_CC_0.09375/justStarrCo_CC_PLINK_IBD/Sex/justStarrCo_CC.sexcheck SEX=4 --output ../data/mexam/AIMs_removed_2/recontruction_$THRESHOLD --threshold $THRESHOLD --ages FILE=../data/mexam/MexAmPheno_041711.txt AGE=5 --affections FILE=../data/mexam/MexAmPheno_041711.txt AFFECTION=4
#./run_PRIMUS.pl --cluster --plink ../data/mexam/AIMs_removed_3/reconstruction_0.375/justStarrCo_CC.genome_no_dups.genome --sexes FILE=../data/mexam/justStarrCo_CC_0.09375/justStarrCo_CC_PLINK_IBD/Sex/justStarrCo_CC.sexcheck SEX=4 --output ../data/mexam/AIMs_removed_4/reconstruction_$THRESHOLD --threshold $THRESHOLD --ages FILE=../data/mexam/MexAmPheno_041711.txt AGE=5 --affections FILE=../data/mexam/MexAmPheno_041711.txt AFFECTION=4

exit;

## AIMs
#./run_PRIMUS.pl --cluster --plink ../data/mexam/AIMs_not_removed_1/justStarrCo_CC_PLINK_IBD/results/justStarrCo_CC.genome --sexes FILE=../data/mexam/justStarrCo_CC_0.09375/justStarrCo_CC_PLINK_IBD/Sex/justStarrCo_CC.sexcheck SEX=4 --output ../data/mexam/AIMs_not_removed_2/recontruction_$THRESHOLD --threshold $THRESHOLD --ages FILE=../data/mexam/MexAmPheno_041711.txt AGE=5 --affections FILE=../data/mexam/MexAmPheno_041711.txt AFFECTION=4

#### JHS
## no AIMs
#./run_PRIMUS.pl --cluster --plink ../data/JHS/AIMs_removed_1/mergedJHS_ARIC_PLINK_IBD/results/mergedJHS_ARIC.genome --sexes FILE=../data/JHS/mergedJHS_ARIC.sexcheck SEX=4 --output ../data/JHS/AIMs_removed_1/recontruction_$THRESHOLD --threshold $THRESHOLD --ages FILE=../data/JHS/JHS_age.txt AGE=3 
## AIMs
#./run_PRIMUS.pl --cluster --plink ../data/JHS/AIMs_not_removed_1/mergedJHS_ARIC_PLINK_IBD/results/mergedJHS_ARIC.genome --sexes FILE=../data/JHS/mergedJHS_ARIC.sexcheck SEX=4 --output ../data/JHS/AIMs_not_removed_1/recontruction_$THRESHOLD --threshold $THRESHOLD --ages FILE=../data/JHS/JHS_age.txt AGE=3 


## Gibbons
#THRESHOLD='0.1875'
#./run_PRIMUS.pl -p ../data/gibbons/justgibbons_030513_nodups.genome_PRIMUS/justgibbons_030513_nodups.genome_network18/justgibbons_030513_nodups.genome_network18.genome --cluster --threshold $THRESHOLD --sexes FILE=../data/gibbons/gibbons_030513_nodups.sexcheck SEX=4 --output ../data/gibbons/justgibbons_030513_nodups.genome_network18.genome_$THRESHOLD\_PRIMUS

#./run_PRIMUS.pl -p ../data/gibbons/justgibbons_030513_nodups.genome_PRIMUS/justgibbons_030513_nodups.genome_network1/justgibbons_030513_nodups.genome_network1.genome --cluster --threshold 0.09375 --sexes FILE=../data/gibbons/gibbons_030513_nodups.sexcheck SEX=4 --output ../data/gibbons/justgibbons_030513_nodups.genome_network1.genome_0.09375\_PRIMUS6

#./run_PRIMUS.pl -p ../data/gibbons/justgibbons_030513_nodups.genome --cluster --threshold $THRESHOLD --sexes FILE=../data/gibbons/gibbons_030513_nodups.sexcheck SEX=4


## CACHE
#./run_PRIMUS.pl -p ../data/cache/cacheibdfinal.genome -o ../data/cache/cacheibdfinal.genome_PRIMUS_$THRESHOLD\_gens2 --cluster --threshold $THRESHOLD --max_gens 2

## GENEVA
#./run_PRIMUS.pl -i FILE=../data/geneva_data/geneva.genome FID2=1 IID2=3 IBD0=4 IBD1=5 IBD2=6 RELATEDNESS=7 -o ../data/geneva_data/geneva.genome_PRIMUS4 --cluster --threshold $THRESHOLD
#./run_PRIMUS.pl -i FILE=../data/geneva_data/geneva_real_ks.txt FID1=1 IID1=1 FID2=2 IID2=2 IBD0=3 IBD1=4 IBD2=5 RELATEDNESS=6 -o ../data/geneva_data/geneva_real_ks.txt_PRIMUS --cluster --threshold 0.05

## Luquetti
#./run_PRIMUS.pl -i FILE=../../../../public_html/luquetti/luquetti_grc_ear_1.genome --no_IMUS --threshold=0.12 -o ../../../../public_html/luquetti/luquetti_grc_ear_1.genome_PRIMUS --sexes FILE=../../../../public_html/luquetti/lookup_luquetti.txt FID=1 IID=1 SEX=3 MALE=Male FEMALE=Female
#./run_PRIMUS.pl -i FILE=../../../../public_html/luquetti/luquetti_grc_ear_1.kinship.homo.kin0 IBD0=6 IBD1=7 IBD2=8 RELATEDNESS=9 --no_IMUS --threshold=0.05 -o ../../../../public_html/luquetti/luquetti_grc_ear_1.kinship.homo.kin0_PRIMUS --sexes FILE=../../../../public_html/luquetti/lookup_luquetti.txt FID=1 IID=1 SEX=3 MALE=Male FEMALE=Female

## SCA
#./run_dataset_cluster.pl sca ../data/sca/justsca.genome ../data/sca/sca.sexcheck 0.1 3 -1

## Oshima
#./run_dataset_cluster.pl oshima2 ../data/oshima/justoshima.genome ../data/oshima/oshima.sexcheck 0.1 3 -1
#./run_PRIMUS.pl --cluster -p ../data/oshima/justoshima.genome --verbose=1 --no_IMUS --output ../data/oshima/justoshima.genome_PRIMUS7 --sexes FILE=../data/oshima/oshima.sexcheck SEX=4 --threshold=0.19

## Reichenberger
#./run_dataset_cluster.pl Reichenberger4 ../data/Reichenberger/justReichenberger.genome ../data/Reichenberger/Reichenberger.sexcheck 0.1 3 -1

## YALE
#./run_dataset_cluster.pl yale/justYale011412_v2_$THRESHOLD ../data/yale/justYale011412.genome ../data/yale/justYale011412.sexcheck $THRESHOLD

## JHS
#THRESHOLD='0.3'
THRESHOLD='0.375'
#THRESHOLD='0.1875'
#THRESHOLD='0.1680'
#THRESHOLD='0.09375'
#./run_PRIMUS.pl --cluster --plink ../data/JHS/justJHSARICnoAIMS.genome --sexes FILE=../data/JHS/mergedJHS_ARIC_IID__IID.sexcheck SEX=4 --output ../data/JHS/justJHSARICnoAIMS_v6_$THRESHOLD --threshold $THRESHOLD --ages FILE=../data/JHS/JHS_age.txt AGE=3

#./run_dataset_cluster.pl JHS/justJHSARICnoAIMS_v5_$THRESHOLD ../data/JHS/justJHSARICnoAIMS.genome ../data/JHS/mergedJHS_ARIC_IID__IID.sexcheck $THRESHOLD 3 -1 ../data/JHS/JHS_age.txt

## MEXAM
#THRESHOLD='0.3'
THRESHOLD='0.375'
#THRESHOLD='0.1875'
#THRESHOLD='0.1680'
#THRESHOLD='0.09375'
#./run_PRIMUS.pl --cluster --plink ../data/mexam/justStarrCo_CC_0.09375/justStarrCo_CC_PLINK_IBD/results/justStarrCo_CC.genome_no_dups.genome --sexes FILE=../data/mexam/justStarrCo_CC_0.09375/justStarrCo_CC_PLINK_IBD/Sex/justStarrCo_CC.sexcheck SEX=4 --output ../data/mexam/justStarrCo_CC_v7_$THRESHOLD\_test --threshold $THRESHOLD --ages FILE=../data/mexam/MexAmPheno_041711.txt AGE=5
#./run_PRIMUS.pl --cluster --plink ../data/mexam/justStarrCo_CC_PRIMUS2/justStarrCo_CC_PLINK_IBD/results/justStarrCo_CC.genome_no_dups.genome --sexes FILE=../data/mexam/justStarrCo_CC_PRIMUS2/justStarrCo_CC_PLINK_IBD/Sex/justStarrCo_CC.sexcheck SEX=4 --output ../data/mexam/justStarrCo_CC_v8_$THRESHOLD --threshold $THRESHOLD --ages FILE=../data/mexam/MexAmPheno_041711.txt AGE=5
#./run_PRIMUS.pl --cluster --plink /nfs/home/below/MexAm/IBDs_ForJEFFtalk/final_noaims_jeanmeth.genome --sexes FILE=../data/mexam/justStarrCo_CC_PRIMUS2/justStarrCo_CC_PLINK_IBD/Sex/justStarrCo_CC.sexcheck SEX=4 --output ../data/mexam/justStarrCo_CC_v9_$THRESHOLD --threshold $THRESHOLD --ages FILE=../data/mexam/MexAmPheno_041711.txt AGE=5


## Jarvick
#./run_dataset_cluster.pl jarvik/justJarvick_v2_$THRESHOLD ../data/jarvik/justJarvick.genome ../data/mexam/consensus_calls.sexcheck $THRESHOLD


## PC_pedigrees
#./run_PRIMUS.pl -p ../data/PCfirstpass.genome -o ../data/PC_pedigrees_v4_$THRESHOLD --sexes FILE=../data/PCfamilies.ped SEX=5 --affections FILE=../data/PCfamilies.ped AFFECTION=6 --cluster --threshold $THRESHOLD
#./run_dataset_cluster.pl PC_pedigrees_v1_$THRESHOLD ../data/PCfirstpass.genome ../data/PCfamilies.ped $THRESHOLD 4 5

## Hapmap
#if [ "1" == "1" ]; then
#	declare -a hapmap_pops=('ASW' 'CEU' 'CHB' 'CHD' 'GIH' 'JPT' 'LWK' 'MEX' 'MKK' 'TSI' 'YRI');
#	for i in "${hapmap_pops[@]}"
#	do:
#		echo "./run_dataset_cluster.pl hapmap3_v4_$THRESHOLD/$i ../data/hapmap3/1kLDpruned_hm3pop$i.genome ../data/hapmap3/hapmap3_r2_b36_fwd.$i.qc.poly.ped.ids $THRESHOLD 4 5"
#		mkdir ../data/hapmap3_v4_$THRESHOLD
#		./run_dataset_cluster.pl hapmap3_v4_$THRESHOLD/$i ../data/hapmap3/1kLDpruned_hm3pop$i.genome ../data/hapmap3/hapmap3_r2_b36_fwd.$i.qc.poly.ped.ids $THRESHOLD 4 5
#	done;
#fi

if [ "1" == "2" ]; then
	#declare -a hapmap_pops=('ASW' 'CEU' 'CHB' 'CHD' 'GIH' 'JPT' 'LWK' 'MEX' 'TSI' 'YRI' 'MKK');
	#declare -a hapmap_pops=('MKK');
	declare -a hapmap_pops=('ASW');
	for i in "${hapmap_pops[@]}"
	do
	:
		#THRESHOLD='0.375'
		#THRESHOLD='0.1875'
		#THRESHOLD='0.1680'
		THRESHOLD='0.09375'
		VERSION='8'
		mkdir ../data/hapmap3_v$VERSION\_$THRESHOLD
		#echo "./run_PRIMUS.pl --ped /net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/data/hapmap3/hapmap3_r2_b36_fwd.$i.qc.poly.ped --map /net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/data/hapmap3/hapmap3_r2_b36_fwd.$i.qc.poly.map --genome --cluster --threshold $THRESHOLD --output ../data/hapmap3_v$VERSION\_$THRESHOLD/$i --ref_pops $i"
		#./run_PRIMUS.pl --ped /net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/data/hapmap3/hapmap3_r2_b36_fwd.$i.qc.poly.ped --map /net/grc/vol1/grapas2/projects/2011/reconstruct_pedigrees/data/hapmap3/hapmap3_r2_b36_fwd.$i.qc.poly.map --genome --cluster --threshold $THRESHOLD --output ../data/hapmap3_v$VERSION\_$THRESHOLD/$i --ref_pops $i --plink_ex ./plink --smartpca_ex ./smartpca
		echo "./run_PRIMUS.pl -p ../data/hapmap3_v8_0.09375/$i/hapmap3_r2_b36_fwd.$i.qc.poly_PLINK_IBD/results/hapmap3_r2_b36_fwd.$i.qc.poly.genome --cluster --threshold $THRESHOLD --output ../data/hapmap3_v$VERSION\_$THRESHOLD/$i --sexes FILE=../data/hapmap3/hapmap3_r2_b36_fwd.$i.qc.poly.ped.ids SEX=5 --no_IMUS"
		./run_PRIMUS.pl -p ../data/hapmap3_v8_0.09375/$i/hapmap3_r2_b36_fwd.$i.qc.poly_PLINK_IBD/results/hapmap3_r2_b36_fwd.$i.qc.poly.genome --cluster --threshold $THRESHOLD --output ../data/hapmap3_v$VERSION\_$THRESHOLD/$i --sexes FILE=../data/hapmap3/hapmap3_r2_b36_fwd.$i.qc.poly.ped.ids SEX=5 --no_IMUS
	done;
fi




