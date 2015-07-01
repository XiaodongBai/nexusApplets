#! /usr/bin/perl

use strict;

my $version = "v1.8.0";

## List all necessary files, lite is a subset of full.
my @bin_files_lite = qw(run_PRIMUS.pl make_dataset_summary.pl make_dataset_pairwise_summary.pl cranefoot cranefoot_mac primus_kickoff7.pl find_expected_pedigree.pl);
my @bin_files;
push(@bin_files,@bin_files_lite);

my @lib_files_lite = qw(KDE_data);
my @lib_files = ();
push(@lib_files, @lib_files_lite);

my @perl_mods_lite = qw(Getopt File Statistics);
my @perl_mods = ();
push(@perl_mods,@perl_mods_lite);

my @PRIMUS_mod_files_lite = qw(compare_fam_files.pm get_age_flags.pm node_v7.pm predict_relationships_2D.pm reconstruct_pedigree_v7.pm IMUS.pm PRIMUS_plus_ERSA.pm);
my @PRIMUS_mod_files = @PRIMUS_mod_files_lite;
push(@PRIMUS_mod_files,"prePRIMUS_pipeline_v7.pm");

my @data_files_lite = qw(1K_genomes_MEX_family.genome 1K_genomes_MEX_family.features complete.genome incomplete.genome complete.fam complete.sex incomplete.sex mt_and_y_halfsib3_size20_sim46-7_v2.genome mt_and_y_halfsib3_size20_sim46-7_v2.sex mt_and_y_halfsib3_size20_sim46-7_v2_MT_estimates.txt mt_and_y_halfsib3_size20_sim46-7_v2_Y_estimates.txt);
my @data_files = qw(MEX_pop.bim MEX_pop.bed MEX_pop.fam MEX_pop.bim mt_and_y_halfsib3_size20_sim46-7_v2.map mt_and_y_halfsib3_size20_sim46-7_v2.ped);
push(@data_files,@data_files_lite);

my @hapmap3_lite = ();
my @hapmap3 = qw(allhapmapUNREL_r2_b36_fwd.qc.poly.fam allhapmapUNREL_r2_b36_fwd.qc.poly.bim allhapmapUNREL_r2_b36_fwd.qc.poly.bed big_snp_annot.txt hapmap3_r2_b36_fwd.*.qc.poly.genome_maximum_independent_set allhapmapUNREL_r2_b36_fwd.qc.poly_SNP_alleles.txt hapmap3_sample_populations.txt);
push(@hapmap3,@hapmap3_lite);


#make_release($version);
make_release("$version\_lite",1);


########################################################################################################################
sub make_release
{

	my $version = shift;
	my $lite = shift;
	my $root = "/nfs/home/grapas2/projects/2011/reconstruct_pedigrees/PRIMUS_$version";
	my $reconstruct_dir = "/nfs/home/grapas2/projects/2011/reconstruct_pedigrees";
	my $bin = "$root/bin";
	my $lib = "$root/lib/";
	my $perl_mod = "$lib/perl_modules";
	my $PRIMUS_mod = "$perl_mod/PRIMUS";
	my $data = "$root/example_data";
	my $hapmap3 = "$lib/hapmap3";
	my $README = "$root/README";	

	if(!-d $root){system("mkdir $root");}
	else{system("rm -r $data/*")}
	if(!-d $bin){system("mkdir $bin");}
	if(!-d $lib){system("mkdir $lib");}
	if(!-d $perl_mod){system("mkdir $perl_mod");}
	if(!-d $hapmap3){system("mkdir $hapmap3");}
	if(!-d $PRIMUS_mod){system("mkdir $PRIMUS_mod");}
	if(!-d $data){system("mkdir $data");}

	if($lite)
	{
		foreach(@bin_files_lite){system("cp -vru ../bin/$_ $bin/")}
		foreach(@lib_files_lite){system("cp -vru ../lib/$_ $lib/")}
		foreach(@perl_mods_lite){system("cp -vru ../lib/perl_modules/$_ $perl_mod/")}
		foreach(@hapmap3_lite){system("cp -vru ../lib/hapmap3/$_ $hapmap3/")}
		foreach(@PRIMUS_mod_files_lite){system("cp -vru ../lib/perl_modules/PRIMUS/$_ $PRIMUS_mod/")}
		foreach(@data_files_lite){system("cp -vru ../example_data/$_ $data/")}

		system("perl -i -pe 's/^use PRIMUS::prePRIMUS_pipeline_v7;//' $bin/primus_kickoff7.pl");
	}
	else
	{
		foreach(@bin_files){system("cp -vru ../bin/$_ $bin/")}
		foreach(@lib_files){system("cp -vru ../lib/$_ $lib/")}
		foreach(@perl_mods){system("cp -vru ../lib/perl_modules/$_ $perl_mod/")}
		foreach(@hapmap3){system("cp -vru ../lib/hapmap3/$_ $hapmap3/")}
		foreach(@PRIMUS_mod_files){system("cp -vru ../lib/perl_modules/PRIMUS/$_ $PRIMUS_mod/")}
		foreach(@data_files){system("cp -vru ../example_data/$_ $data/")}
	}
	system("cp -vru ../bin/README $README");

	chdir $reconstruct_dir;
	system("tar -cvzf PRIMUS_$version.tgz PRIMUS_$version");
	chdir "$reconstruct_dir/bin";
}

