package PRIMUS::find_expected_pedigree;

use strict;
use lib "../PRIMUS_dev/lib/perl_modules";
use PRIMUS::compare_fam_files;
use PRIMUS::node_v7;

sub find_expected_pedigree
{
	my $summary_file = shift;
	my $ref_ped = shift;
	my $verbose = shift;

	if($verbose eq ""){$verbose = 0}
	$summary_file =~ /^(.*)Summary_([\w\.-]+).txt$/;
	my $dataset_dir = $1;
	my $dataset_name = $2;

	my $USAGE = "\n\n$0 [dataset_summary_file] [reference .fam file]\n\n\n"; 

	my $PASS = 0;

	if(!-e $summary_file){die "Summary file $_ does not exist:$USAGE";}
	if(!-e $ref_ped){die "Your reference pedigree $_ does not exist. USAGE:$USAGE";}
	if($verbose > 0)
	{
		print "dataset_dir = $dataset_dir\n";
		print "Ref_ped file = $ref_ped\n";
	}

	my %ref_hash = load_ref_ped_file($ref_ped);
	my $unrelateds = load_unrelated_samples("$dataset_dir/$dataset_name\_unrelated_samples.txt");
	
	my @ref_names = keys (%ref_hash);
	
	## Find the relevant family network
	my ($networks_ref,$samples_ref) = get_relevant_network($summary_file,@ref_names);
	my @networks = @$networks_ref;
	my %network_names = %$samples_ref;


	## Check that the reference pedigree is larger than one sample and does not have an individual in the unrelated_samples.txt file
	foreach(@ref_names)
	{
		if(@ref_names > 1)
		{
			if(exists $$unrelateds{$_})
			{
				if($verbose > 0)
				{
					print "sample $_ from your reference pedigree is unrelated to all other sample in this dataset; therefore, the reference pedigree is not among those reconstructed\n";
					print "NO MATCH found\n";
				}
				if(@networks < 1){push(@networks,"NA")}
				return (0,join(',',@networks),0,"NA","$_\_predicted_unrelated_to_all_other_samples");
			}
		}
		elsif(@ref_names eq 1)
		{
			if(!exists $$unrelateds{$_})
			{
				if($verbose > 0)
				{
					print "sample $_ from your reference pedigree is related to at least one other sample in this dataset; therefore, the reference pedigree is not among those reconstructed\n";
					print "NO MATCH found\n";
				}
				return (0,join(',',@networks),0,"NA","$_\_is_related_to_other_samples");
			}
			else
			{
				if($verbose > 0)
				{
					print "Your pedigree of one person contains an individual ($_) that is unrelated to everyone else in this dataset. \nPedigree confirmed\n";
					print "MATCH found\n";
				}
				return (1,join(',',@networks),0,"NA","NA");
			}
		}
	}

	# TEST A: ref_ped only matches to one family network
	if(@networks > 1)
	{
		if($verbose > 0)
		{
			print "Multiple family networks match you samples: @networks. Your pedigree is not one of the possible.\n";
		}
		return (0,join(',',@networks),0,"NA","Samples_from_expected_PED_not_all_in_same_network");
	}
	if(@networks < 1)
	{
		if($verbose > 0)
		{
			print "No network found in this dataset summary file with samples that match you pedigree\n";
		}
		return (0,"NA",0,'NA',"No_network_found_in_this_dataset_summary_file_with_samples_that_match_you_pedigree");
	}
	if($verbose > 0)
	{
		print "Family network of interest: @networks\n";
	}

	## Find the intersection of non-missing samples between ref-ped and the network
	my @intersection = get_intersection(\%ref_hash,$network_names{'last'});
	my %intersection = map{$_,1}@intersection;
	if(@ref_names < @intersection){die "ERROR!!! ref_names (@ref_names) has more samples than intersection (@intersection).\n";}
	#print "\n\nintersection: @intersection\n";


	## if necessary, convert !intersection samples to missing in ref-ped; REMOVE extra missing samples
	convert_non_intersection_to_missing(\%intersection,\%ref_hash);
	#print "\n";
	#foreach(keys %ref_hash){print "$_; $ref_hash{$_}{'PID'}; $ref_hash{$_}{'MID'}\n";}
	remove_unnecessary_missing(\%ref_hash);
	@ref_names = keys %ref_hash;
	#print "\nPruned Reference ped\n";
	open(PP,">$ref_ped\_pruned.fam");
	foreach(keys %ref_hash)
	{
		#print PP "$ref_hash{$_}{'FID'}\t$_\t$ref_hash{$_}{'PID'}\t$ref_hash{$_}{'MID'}\t$ref_hash{$_}{'SEX'}\t$ref_hash{$_}{'AFF'}\n";
		#print "$_; $ref_hash{$_}{'PID'}; $ref_hash{$_}{'MID'}\n";
	}
	close(PP);

	## foreach possible pedigree in family network, 
	##   if necessary, convert !intersection samples to missing and remove extra missing samples
	##   Compare pedigree to new Ref_ped
	my $ped_num = 1;

	my $ped_file = "$dataset_dir/$dataset_name\_network@networks[0]/$dataset_name\_network@networks[0]_$ped_num.fam";
	print "ped_file: $ped_file\n";
	if(!-e $ped_file)
	{
		if($verbose > 0)
		{
			print "Not valid .ped file for network@networks[0]: $ped_file\n";
		}
		return (0,join(',',@networks),0,"NA","Reconstruction_failed");
	}
	while(-e $ped_file)
	{
		if($verbose > 0)
		{
			#print "ped file: $ped_file\n";
		}
		my %ped = load_reconstructed_ped_file($ped_file);
		my $converted = convert_non_intersection_to_missing(\%intersection,\%ped);
		#print "\n";
		#foreach(keys %ped){print "$_; $ped{$_}{'PID'}; $ped{$_}{'MID'}\n";}
		remove_unnecessary_missing(\%ped);
		#print "\nReconstructed pruned ped\n";
		#foreach(keys %ped){print "$_; $ped{$_}{'PID'}; $ped{$_}{'MID'}\n";}

		if(PRIMUS::compare_ped_files::are_peds_same(\%ped,\%ref_hash))
		{
			if($verbose > 0)
			{
				print "MATCH!!! network@networks[0] ped $ped_num: $ped_file\n";
			
			}
			if($converted)
			{
				print "Your reference pedigree was merged with other samples in the dataset\n";
			}
			#return $ped_file;
			return (1,join(',',@networks),$converted,$ped_file,$ped_num);
		}	

		$ped_num++;
		$ped_file = "$dataset_dir/$dataset_name\_network@networks[0]/$dataset_name\_network@networks[0]_$ped_num.fam";
	}

	if($verbose > 0)
	{
		print "NO MATCHING PEDIGREE FOUND\n";
	}
	return (0,join(',',@networks),0,"NA","NO_MATCHING_PEDIGREE_FOUND");
	#return 0;
}

#########################################################
#########################################################
sub remove_unnecessary_missing
{
	my $hash = shift;

	my $network_ref = make_network_from_ped($hash);
	my $remove = 1;
	while($remove eq 1)
	{
		$remove = 0;
		foreach my $node_name (keys %$network_ref)
		{
			#if(!exists $$hash{$node_name}){next;}
			#my @parents = ($$hash{$node_name}{'PID'},$$hash{$node_name}{'MID'});
			#if(@parents eq 0 || @parents[0] !~ /Missing/ || @parents[1] !~ /Missing/){next;}

			#my @parents0 = ($$hash{@parents[0]}{'PID'},$$hash{@parents[0]}{'MID'});
			#my @parents1 = ($$hash{@parents[1]}{'PID'},$$hash{@parents[1]}{'MID'});
			#if(@parents0 > 0 || @parents1 > 0){next;}
			
			#print "\nnode: $node_name\n";

			if(!exists $$network_ref{$node_name}){next;}
			
			my @children = $$network_ref{$node_name}->children();
			my @parents = $$network_ref{$node_name}->parents();
			#print "parents: @parents\n";
			#print "children: @children\n";

			## Remove is node is Missing and has no children
			if($node_name =~ /Missing/i && @children < 1)
			{
				#print "here\n";
				if(@parents > 0)
				{
					$$network_ref{@parents[0]}->remove_child($node_name);
					$$network_ref{@parents[1]}->remove_child($node_name);
				}

				#my @children0 = $$network_ref{@parents[0]}->children();
				#print "children0: @children0\n";
				#my @children1 = $$network_ref{@parents[1]}->children();
				#print "children1: @children1\n";
				
				delete $$network_ref{$node_name};
				delete $$hash{$node_name};
				$remove = 1;
				next;
			}

			#print "parents: @parents\n";
			if(($node_name !~ /Missing/i && @parents eq 0) || @parents[0] !~ /Missing/i || @parents[1] !~ /Missing/i){next;}
			
			my @parents0 = $$network_ref{@parents[0]}->parents();
			my @parents1 = $$network_ref{@parents[1]}->parents();
			#print "parents0: @parents0\n";
			#print "parents1: @parents1\n";
			if(@parents0 > 0 || @parents1 > 0){next;}
			
			my @children0 = $$network_ref{@parents[0]}->children();
			my @children1 = $$network_ref{@parents[1]}->children();
			#print "children0: @children0\n";
			#print "children1: @children1\n";
			if(@children0 > 1 || @children1 > 1){next;}
		
			## Passed all criteria, remove dummy parents
			$remove = 1;
			#print "**REMOVE** parents: @parents\n";
			$$network_ref{$node_name}->delete_parents();
			$$hash{$node_name}{'PID'} = 0;
			$$hash{$node_name}{'MID'} = 0;
			#$$network_ref{@parents[0]}->delete_children();
			#$$network_ref{@parents[1]}->delete_children();
			delete $$network_ref{@parents[0]};
			delete $$network_ref{@parents[1]};
			delete $$hash{@parents[0]};
			delete $$hash{@parents[1]};
		}
	}
}

sub make_network_from_ped
{
	my $hash = shift;
	my %network;
	my $network_ref = \%network;
	
	#print "here\n";

	foreach my $child (keys %$hash)
	{
		#print "child: $child\n";
		my $PID = $$hash{$child}{'PID'};
		my $MID = $$hash{$child}{'MID'};
		my $dad = $$hash{$child}{'PID'};
		my $mom = $$hash{$child}{'MID'};
		if(!exists $$network_ref{$child})
		{
			my $node = new PRIMUS::node_v7($child);
			$$network_ref{$child} = $node;
		}
		if(!exists $$network_ref{$dad} && $PID ne 0)
		{
			my $node = new PRIMUS::node_v7($dad);
			$$network_ref{$dad} = $node;
		}
		if(!exists $$network_ref{$mom} && $MID ne 0)
		{
			my $node = new PRIMUS::node_v7($mom);
			$$network_ref{$mom} = $node;
		}
		if($PID ne 0)
		{
			$$network_ref{$child}->add_parent($dad);
			$$network_ref{$dad}->add_child($child);
		}
		
		if($MID ne 0)
		{
			$$network_ref{$child}->add_parent($mom);
			$$network_ref{$mom}->add_child($child);
		}
	}
	return $network_ref;
}

sub convert_non_intersection_to_missing
{
	my $intersection = shift;
	my $ped = shift;
	my $ref_dummy_ctr = 15000;
	my $converted = 0;
	foreach(keys %$ped)
	{
		#print "name: $_\n";
		if(!exists $$intersection{$_} && $_ !~ /Missing/i)
		{
			#print "make Missing: $_ to Missing$ref_dummy_ctr\n";
			change_name($_,"Missing$ref_dummy_ctr",$ped);
			$ref_dummy_ctr++;
			$converted = 1;
		}
	}
	return $converted;
}

sub change_name
{
	my $old = shift;
	my $new = shift;
	my $ped = shift;

	if($old =~ /Missing/i){return};
	#print "Converting $old to $new\n";
	
	foreach my $name (keys %$ped)
	{
		if($$ped{$name}{'PID'} eq $old){$$ped{$name}{'PID'} = $new}
		if($$ped{$name}{'MID'} eq $old){$$ped{$name}{'MID'} = $new}
		if($name eq $old)
		{
			$$ped{$new} = $$ped{$name};
			delete $$ped{$name};
		}
	}
}

sub load_unrelated_samples
{
	my $file = shift;
	my %unrelateds;
	
	open(IN,$file) or die "can't open $file; $!\n";
	<IN>;
	while(my $line = <IN>)
	{
		my ($FID,$IID) = split(/\s+/,$line);
		$unrelateds{$IID} = 1;
	}
	
	return \%unrelateds;
}	

sub load_ref_ped_file
{
	my $ctr = 1000;
	my $file = shift;
	open(IN,$file);
	my %ped;
	my %all_samples;
	while(my $line = <IN>)
	{
		while($line =~ /\s0\s/)
		{
			$line =~ s/\s0\s/\tMissing$ctr\t/;
			$ctr++;
		}
		$line =~ s/\s0\n/\tMissing$ctr\n/;
		chomp($line);
		$ctr++;
		my ($FID,$IID,$PID,$MID,$sex,$aff,@rest) = split(/\s+/,$line);
		#print "$FID,$IID,$PID,$MID\n";
		$ped{$IID}{"FID"} = $FID;
		$ped{$IID}{"PID"} = $PID;
		$ped{$IID}{"MID"} = $MID;
		$ped{$IID}{"SEX"} = $sex;
		$ped{$IID}{"AFF"} = $aff;
		$all_samples{$IID} = 1;
		$all_samples{$PID} = 1;
		$all_samples{$MID} = 1;
	}
	close(IN);
	
	## If a parent is not and IID, add it
	foreach(keys %all_samples)
	{
		if(!exists $ped{$_})
		{
			$ped{$_}{'FID'} = "?";
			$ped{$_}{'PID'} = "Missing$ctr";
			$ctr++;
			$ped{$_}{'MID'} = "Missing$ctr";
			$ctr++;
		}
	}
	return %ped;
}

sub load_reconstructed_ped_file
{
	my $file = shift;
	open(IN,$file);
	my %ped;
	my $ctr = 5000;
	while(my $line = <IN>)
	{
		while($line =~ /\s0\s/)
		{
			$line =~ s/\s0\s/\tMissing$ctr\t/;
			$ctr++;
		}
		if($line =~ /\s0\n/)
		{
			$line =~ s/\s0\n/\tMissing$ctr\n/;
			$ctr++;
		}
		chomp($line);
		#print "line: $line\n";
		my ($network,$name,$PID_long,$MID_long,@rest) = split(/\s+/,$line);
		my ($FID,$IID); 
		
		if($name !~ /Missing/i && $name !~ /Missing/i)
		{
			($FID,$IID) = split(/__/,$name);
		}
		else
		{
			$IID = $name;
			$FID = "?";
		}
		my ($PID,$MID,$junk);
		if($PID_long !~ /Missing/i && $PID_long !~ /Missing/i)
		{
			($junk,$PID) = split(/__/,$PID_long);
		}
		else
		{
			$PID = $PID_long;
		}
		if($MID_long !~ /Missing/i && $MID_long !~ /Missing/i)
		{
			($junk,$MID) = split(/__/,$MID_long);
		}
		else
		{
			$MID = $MID_long;
		}
		#print "$FID,$IID,$PID,$MID\n";
		$ped{$IID}{"FID"} = $FID;
		$ped{$IID}{"PID"} = $PID;
		$ped{$IID}{"MID"} = $MID;
	}
	close(IN);
	return %ped;
}

sub get_intersection
{
	my $ref1 = shift;
	my $ref2 = shift;

	my @intersection;
	foreach(keys %$ref1)
	{
		if(exists $$ref2{$_})
		{
			push(@intersection,$_);
		}
	}
	return @intersection;
}


sub get_relevant_network
{
	my $summary_file = shift;
	my @ref_names = @_;
	my @matching_networks;
	my %samples;
	
	## Check the summary file
	if(!-e $summary_file){die "Cannot find $summary_file; provide full path to a PRIMUS dataset summary file (NOT a network summary file)\n";}

	#print "Dataset name: $dataset_name\n";
	#print "REF_NAMES: @ref_names\n";

	open(IN, $summary_file);
	while(my $line = <IN>)
	{
		#print "line: $line";
		my $summary_sample_column = 6; ## zero based
		chomp($line);
		if($line !~ /^\d+\s/){next;}
		my @temp = split(/\s+/,$line);
		my @long_names = split(/,/,@temp[$summary_sample_column]);
		my @names;
		foreach(@long_names)
		{
			my($FID,$IID) = split(/__/,$_);
			push(@names,$IID);
		}
		#print "names: @names\n\n";

		my %names = map{$_,1}@names;

		my $any_match = 0;
		my $full_match = 1;
		foreach(@ref_names)
		{
			## If at least one name matches, add it as a network
			if(exists $names{$_})
			{
				$any_match = 1; 
				push(@matching_networks,@temp[0]);
				$samples{@temp[0]} = \%names;
				$samples{'last'} = \%names;
				last;
			}
		}
	}
	close(IN);

	return (\@matching_networks,\%samples);
}

sub do_arrays_match
{
	my $arr1_ref = shift;
	my $arr2_ref = shift;

	if(@$arr1_ref ne @$arr2_ref)
	{
		return 0;
	}
	foreach (@$arr1_ref)
	{
		my $val = $_;
		if(!grep($_ eq $val,@$arr2_ref))
		{
			return 0;
		}
	}
	foreach (@$arr2_ref)
	{
		my $val = $_;
		if(!grep($_ eq $val,@$arr1_ref))
		{
			return 0;
		}
	}
	return 1;
}

return 1;
