package Gerald::Jerboa;

BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT = qw(
      &info 
      &param_file_is_legacy_style 
      &readFileXML 
      &key_matched
      &calc_total_by_field 
      &getBustardPath 
      &extract_run_info 
      &calc_percent_by_field 
      &derive_uniqueness_rescue_stats 
      &init_lane_results
      &read_config_file 
      &process_pair_info_file 
      &sum_vals_under_node 
      &derive_mispairing_stats
      &derive_insert_freq_stats 
      &derive_pair_stats &read_tiles_file 
      &process_score_file
      &adjust_measured_phasings
      &store_applied_phasings 
      &store_measured_phasings 
      &read_sample_file
      &legacy_read_applied_phasing_file 
      &getIndexRead
      &deleteIndexingRead
      &get_cluster_counts
      &get_tile_cluster_counts
      &merge_lanes
      &merge_tiles
      $maxNumLanes
    );
    @EXPORT_OK = qw();
}

#------------------------------------------------------------------------------

# PROJECT: GERALD
# MODULE:  Jerboa.pm
# AUTHOR:  M. Zerara
# functions taken from the old jerboa.pl
# written by L. J. Davies and extended by by R. J. Shaw 
#
# Copyright (c) 2008-2010 Illumina Inc.
#
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).
#
# This file is part of the Consensus Assessment of Sequence And VAriation
# (CASAVA) software package.
#
# Functions common to jerboa.pl component

use Data::Dumper;
use warnings FATAL => 'all';
use strict;
use POSIX;
use XML::Simple;
use File::Basename;
use constant DEBUG => 0; # set to 1 to get debug info

use Casava::Common::Utils;


sub info($);
sub param_file_is_legacy_style($);
sub readFileXML($);
sub key_matched($;$);
sub calc_total_by_field($;$;$;$);
sub getBustardPath($;$);
sub extract_run_info($$$$$$);
sub calc_percent_by_field($;$;$;$);
sub derive_uniqueness_rescue_stats($;$;$);
sub init_lane_results($;$;$;$;$;$);
sub read_config_file($;$;$;$;$;$;$;$;$);
sub process_pair_info_file($;$);
sub sum_vals_under_node($);
sub derive_mispairing_stats($;$);
sub derive_insert_freq_stats($);
sub derive_pair_stats($);
sub read_tiles_file($;$);
sub process_score_file($;$;$;$;$;$;$;$;$;$;$;$);
sub adjust_measured_phasings($$$$$);
sub store_applied_phasings($$$$$$$$);
sub store_measured_phasings($$$$$$);
sub read_sample_file($$$$);
sub legacy_read_applied_phasing_file($$$$$);
sub getIndexRead($);
sub deleteIndexingRead($$$);
sub get_cluster_counts($;\%);
sub get_tile_cluster_counts(\%;$;$;$;\$;\$);
sub merge_lanes($$$$$$$$$);
sub merge_tiles($$$$$$$);


my $paired_not_done_key = 'NoPairedAlignmentDone';
my $Insert_Size_str = 'InsertSize';
my @insert_stat_display_map = (['Median', 'Median'],
                               ['LowSD', 'Below-median SD'],
                               ['HighSD', 'Above-median SD'],
                               ['Min', 'Low thresh.'],
                               ['Max', 'High thresh.']);
my $not_avail_str = 'N/A';
my $Total_str = 'Total';

my $Orientation_str = 'Orientation';
my @orient_display_map = (['Fm', 'F-: &gt; R2 R1 &gt;'],
                          ['Fp', 'F+: &gt; R1 R2 &gt;'],
                          ['Rm', 'R-: &lt; R2 R1 &gt;'],
                          ['Rp', 'R+: &gt; R1 R2 &lt;']);

my $nominal_orient_key = 'Nominal';

my $num_orient_ok_small_insert_key = 'NominalOrientationButSmallInsert';
my $num_orient_ok_large_insert_key = 'NominalOrientationButLargeInsert';
my $num_pairs_ok_key = 'PairsOK';

my @inserts_display_map = ([$num_orient_ok_small_insert_key, 'Too small'],
                           [$num_orient_ok_large_insert_key, 'Too large'],
                           [$num_pairs_ok_key, 'Orientation and size OK']);


my $insert_pc_format = "%0.1f";

my $pc_suffix = 'Percent';
my $pc_prefix = 'pc';
my $pc_orient_ok_small_insert_key
    = $num_orient_ok_small_insert_key . $pc_suffix;
my $pc_orient_ok_large_insert_key
    = $num_orient_ok_large_insert_key . $pc_suffix;
my $pc_pairs_ok_key = $num_pairs_ok_key . $pc_suffix;


my %insert_num_pc_key_map = ($num_orient_ok_small_insert_key
                             => $pc_orient_ok_small_insert_key,
                             $num_orient_ok_large_insert_key
                             => $pc_orient_ok_large_insert_key,
                             $num_pairs_ok_key
                             => $pc_pairs_ok_key);


my $Unique_key = 'SingleAlignmentFound';
my $Unique_str = 'Unique';
my $Rescuable_key = 'ManyAlignmentsFound';
my $Rescuable_str = 'Rescuable';
my $Repeat_key = 'Repeat';
my $Repeat_str = $Repeat_key;
my $Repeat_Masked_key = 'RM';
my $Repeat_Masked_str = 'Repeat Masked';
my $Not_Matched_key = 'NM';
my $Not_Matched_str = 'Not Matched';
my $Low_Quality_key = 'QC';
my $Low_Quality_str = 'Low Quality';
my $Read_str = 'Read';

my $unique_pair_align_key = 'UniquePairedAlignment';
my $multi_pair_align_key = 'ManyPairedAlignments';

my @align_outcome_display_map = ([$Unique_key, $Unique_str,],
                                 [$Rescuable_key, $Rescuable_str],
                                 [$Repeat_key, $Repeat_str],
                                 [$Repeat_Masked_key, $Repeat_Masked_str],
                                 [$Not_Matched_key, $Not_Matched_str],
                                 [$Low_Quality_key, $Low_Quality_str]);


my $Lost_str = 'Lost';
my $pair_key_separ_str = ':';

my $num_paired_reads = 2;
our $maxNumLanes = 8;

# BustardSummary.xml keys
my $cluster_count_raw_key = 'clusterCountRaw';
my $cluster_count_PF_key = 'clusterCountPF';

sub info($)
{
    my $message = shift @_;
    print STDERR "Info: ", basename($0), $message, "\n";
}

#------------------------------------------------------------------------------

sub store_measured_phasings($$$$$$)
{
    my $lane_results_ref = shift;
    my $bustard_path = shift;
    my $num_lanes = shift;
    my $num_reads = shift;
    my $phasing_options_ref = shift;
    my $param_is_legacy = shift;

    my $phasing_path = "$bustard_path/Phasing";
    my $xml_data;

    my $measured_phasing_file;

    for (my $lane_num = 1; $lane_num <= $num_lanes; ++$lane_num) {
        for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
            my $lane_read_ref = $lane_results_ref->[$read_num][$lane_num];

            if (!$param_is_legacy) {
                my $curr_options_ref = $phasing_options_ref->{$read_num};

                # FIXME : error?
                next if (!defined($curr_options_ref));

                # Take a copy so AutoFlag can be overridden.
                my %curr_options = %{$curr_options_ref};
                $curr_options{'AutoFlag'} = 2; # `use measured per lane'

                $measured_phasing_file = getPhasingFilename($phasing_path,
                                                            \%curr_options,
                                                            $lane_num);

            } else {
                # Look for lane-specific measured phasing values 
                # (pipeline ver >= 0.1.7).
                $measured_phasing_file = sprintf("%s/s_%d_phasing.xml",
                                                 $phasing_path, $lane_num);
            }

            if (defined($xml_data = readFileXML($measured_phasing_file))) {
                if (defined($xml_data->{Phasing})) {
                    $lane_read_ref->{phasingmeasured} = $xml_data->{Phasing};
                }

                if (defined($xml_data->{Prephasing})) {
                    $lane_read_ref->{prephasingmeasured}
                    = $xml_data->{Prephasing};
                }
            }
        }
    }
}

#------------------------------------------------------------------------------

sub adjust_measured_phasings($$$$$)
{
    my $results_ref = shift;
    my $bustard_path = shift;
    my $num_lanes = shift;
    my $phasing_options_src = shift;
    my $phasing_options_dest = shift;

    # last read from original into last read of demux bin
    my $from = (sort (keys %{$phasing_options_src}))[-1];
    my $to   = (sort (keys %{$phasing_options_dest}))[-1];
    my $phasing_path = "$bustard_path/Phasing";

    # Take a copy so AutoFlag can be overridden.
    my %curr_options_src  = %{$phasing_options_src->{$from}};
    my %curr_options_dest = %{$phasing_options_dest->{$to}};
    # `use measured per lane'
    $curr_options_src{'AutoFlag'} = $curr_options_dest{'AutoFlag'} = 2;

    for (my $lane_num = 1; $lane_num <= $num_lanes; ++$lane_num)
    {
        my $phasing_lane_src  = getPhasingFilename($phasing_path,\%curr_options_src,$lane_num);
        my $phasing_lane_dest = getPhasingFilename($phasing_path,\%curr_options_dest,$lane_num);
        $results_ref->{$phasing_lane_src} = $phasing_lane_dest;
    }
    $curr_options_src{'AutoFlag'} = $curr_options_dest{'AutoFlag'} = 1;
    my $phasing_total_src  = getPhasingFilename($phasing_path,\%curr_options_src);
    my $phasing_total_dest = getPhasingFilename($phasing_path,\%curr_options_dest);
    $results_ref->{$phasing_total_src} = $phasing_total_dest;
}

#------------------------------------------------------------------------------

sub legacy_read_applied_phasing_file($$$$$)
{
    my $app_name = shift;
    my $phasings_applied_file = shift;
    my $phasing_ref = shift;
    my $prephasing_ref = shift;
    my $bustard_dir = shift;

    my $xmlData;

    if (defined($xmlData = readFileXML($phasings_applied_file))) {
        if (defined($xmlData->{Run})) {
            # FIXME : Use forcearray and keyattr in readFileXML to avoid 
            #         the ARRAY / HASH division. 

            my $run_ref = undef;

            if (ref($xmlData->{Run}) eq 'ARRAY') {
                for (my $i = 0; $i < scalar(@{$xmlData->{Run}}); $i++) {
                    if ($xmlData->{Run}[$i]{Name} eq $bustard_dir) {
                        $run_ref = $xmlData->{Run}[$i];
                        last;
                    }
                }
            } elsif (ref($xmlData->{Run}) eq 'HASH') {
                $run_ref = $xmlData->{Run};
            }

            $$phasing_ref = $run_ref->{Parameters}[0]
            if (defined($run_ref->{Parameters}[0]));

            $$prephasing_ref = $run_ref->{Parameters}[1]
            if (defined($run_ref->{Parameters}[1]));

            printf STDERR "%s: found phase=%5.3f prephase=%5.3f\n",
            $app_name, $$phasing_ref, $$prephasing_ref;
        }
    }
}

#------------------------------------------------------------------------------

sub store_applied_phasings($$$$$$$$)
{
    my $app_name = shift;
    my $lane_results_ref = shift;
    my $num_lanes = shift;
    my $num_reads = shift;
    my $bustard_path = shift;
    my $configXMLFile = shift;
    my $phasing_options_ref = shift;
    my $param_is_legacy = shift;

    my $phasing_applied = undef;
    my $prephasing_applied = undef;
    my $phasings_applied_file = $configXMLFile;
    my $bustard_dir = basename($bustard_path);

    if ($param_is_legacy) {
        legacy_read_applied_phasing_file($app_name, $phasings_applied_file,
                                         \$phasing_applied,
                                         \$prephasing_applied,
                                         $bustard_dir);
    }

    my $phasing_dir = 'Phasing';
    my $xml_data;
    my $curr_options_ref;
    my %legacy_options;

    if ($param_is_legacy) {
        $legacy_options{'AutoFlag'} = 0;
        $legacy_options{'Path'} = '';
        $curr_options_ref = \%legacy_options;
    }

    for (my $lane_num = 1; $lane_num <= $num_lanes; ++$lane_num) {
        for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
            my $lane_read_ref = $lane_results_ref->[$read_num][$lane_num];

            if ($param_is_legacy) {
                if ($phasing_applied > 0.) {
                    $lane_read_ref->{phasingApplied} = $phasing_applied;
                }

                if ($prephasing_applied > 0.) {
                    $lane_read_ref->{prephasingApplied} = $prephasing_applied;
                }
            }

            # Applied phasing files now provided by getPhasingFilename.

            if (!$param_is_legacy) {
                $curr_options_ref = $phasing_options_ref->{$read_num};

                # FIXME : error?
                next if (!defined($curr_options_ref));
            }

            my $applied_phasing_file
                = join('/', ($bustard_path,
                             getPhasingFilename($phasing_dir,
                                                $curr_options_ref,
                                                $lane_num)));

            if (defined($xml_data = readFileXML($applied_phasing_file))) {
                if (defined($xml_data->{Phasing})) {
                    $lane_read_ref->{phasingApplied} = $xml_data->{Phasing};
                }

                if (defined($xml_data->{Prephasing})) {
                    $lane_read_ref->{prephasingApplied}
                    = $xml_data->{Prephasing};
                }
            }
        }
    }
}

#------------------------------------------------------------------------------

sub read_sample_file($$$$)
{
    my $appName = shift;
    my $sampleFile = shift;
    my $laneResultsRef = shift;
    my $num_reads = shift;

    my $xmlData;
    my $unknown='unknown';

    if (defined($xmlData = readFileXML($sampleFile))) {

        for (my $i = 0; $i < scalar(@{$xmlData->{Lane}}); $i++) {

            for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
                $laneResultsRef->
                    [$read_num][$xmlData->{Lane}[$i]{laneNumber}]{sample}
                = $xmlData->{Lane}[$i]{SampleID};
            }
        }
    } else {
        warn "Unable to find file $sampleFile";
        for (my $lane = 1; $lane <= $maxNumLanes; $lane++) {
            for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
                $laneResultsRef->[$read_num][$lane]{sample} = $unknown;
            }
        }
    }
}

#------------------------------------------------------------------------------


sub param_file_is_legacy_style($)
{
    my $bustard_dir = shift;

    # This was designed to test the format of the phasing  files used before
    # 0.3, therefore it is now irrelevant
    return 0;

    my %dir_info;

    if (!parseAnalysisDirName($bustard_dir, \%dir_info)) {
        warn "Failed to parse bustard directory name : $bustard_dir\n";
        exit(1);
    }

    if (!(defined($dir_info{major_ver}) && defined($dir_info{minor_ver})
          && defined($dir_info{release}))) {
        warn "Failed to extract bustard version from : $bustard_dir\n";
        exit(1);
    }

    my @bustard_ver = ($dir_info{major_ver}, $dir_info{minor_ver},
                       $dir_info{release});
    my @new_param_bustard_ver = (1, 9, 0);

    my $style_is_legacy
        = (compareVersions(\@new_param_bustard_ver, \@bustard_ver) < 0);

    return $style_is_legacy;
}

# Read XML file into a hash. Return 1 if successfully completed
# Check for zero padding at end of file
sub readFileXML($)
{
    my ($fileName)=@_;
    return undef unless (-f $fileName);
    unless(open (FILE, $fileName))
    {
        warn "Failed to open $fileName: $!"; return undef;
    }
    my $text="";
    while (<FILE>)
    {
        $text.=$_;
    } # while
    close (FILE);

    my $l=length($text);

    while (($l!=0)&&(substr($text,$l-1,1) eq "\0"))
    {
        $l--;
    } # while

    if ($l==0)
    {
        warn "file $fileName seems to be completely blank";
        return undef;
    } # if   

    $text=substr($text,0,$l) unless ($l==length($text));

    my $xml = XML::Simple->new();
    my $data = $xml->XMLin($text);

    return $data;
} # sub readFileXML

sub key_matched($;$)
{
    my $match_key_fields_ref = shift;
    my $key_fields_ref = shift;

    my @match_key_fields = @{$match_key_fields_ref};
    my $num_match_key_fields = scalar(@match_key_fields);

    my @key_fields = @{$key_fields_ref};

    if (scalar(@key_fields) != $num_match_key_fields) {
        return 0;
    }

    my $matched = 1;

    for (my $field_ind = 0; $field_ind < $num_match_key_fields;
         ++$field_ind) {
        if (!($key_fields[$field_ind]
              =~ /$match_key_fields[$field_ind]/)) {
            $matched = 0;
            last;
        }
    }

    return $matched;
}

sub calc_total_by_field($;$;$;$)
{
    my $hash_ref = shift;
    my $match_key_fields_ref = shift;
    my $wild_field_ind = shift;
    my $separ_str = shift;

    foreach my $key_str (keys %{$hash_ref}) {
        my @key_fields = split(/$separ_str/, $key_str);
        next if (!key_matched($match_key_fields_ref, \@key_fields));

        my @total_key_fields = @key_fields;
        splice(@total_key_fields, $wild_field_ind, 1, $Total_str);
        my $total_key = join($separ_str, @total_key_fields);

        if (!defined($hash_ref->{$total_key})) {
            $hash_ref->{$total_key} = 0;
        }

        my $val = $hash_ref->{$key_str};

        if (!defined($val) || ($val eq $not_avail_str)) {
            $val = 0;
        }

        $hash_ref->{$total_key} += $val;
    }
}

sub getBustardPath($;$) {
   my $bustardPathRef = shift;
   my $configFile = shift;
   my $xml = XML::Simple->new();
   my $xmlData = $xml->XMLin($configFile);

   if (!defined($xmlData)) {
     warn "Failed to parse $configFile as XML\n";
     exit(1);
   }
     $$bustardPathRef = $xmlData->{'Defaults'}{"EXPT_DIR"} or die "Can't find EXPT_DIR in xml file!";
}

sub extract_run_info($$$$$$)
{
    my $appName = shift;
    my $bustardPath = shift;
    my $runFolderRef = shift;
    my $machineNameRef = shift;
    my $bustardDirRef = shift;
    my $tileAreaRef = shift;

    my $paramsfile    = getConfigurationFilePath($bustardPath);
    my @path          = grep{ $_ } File::Spec->splitdir($bustardPath);
    my @emptyArray    = ();
    my %emptyHashMap  = ();
    my $runParameters =  getBustardRunInfo($paramsfile, $path[-1], 'RunParameters', @emptyArray, %emptyHashMap);

    my $runFolder= $runParameters->{'RunFolder'} if $runParameters and exists $runParameters->{'RunFolder'};
    $runFolder= 'unknown' unless $runFolder;
    my $machineName = $runParameters->{'Instrument'} if $runParameters and exists $runParameters->{'Instrument'};
    $machineName = $1 if not $machineName and ($runFolder =~ /^\d{6}_([^_]*)_\d+/);

    unless (defined($machineName)) {
        my $runFolderPath = File::Spec->catfile($bustardPath, '..', '..', '..');
        my $paramsFile = "$runFolderPath/$runFolder.params";

        my $xmlData;
        if (defined($xmlData = readFileXML($paramsFile))) {

            if (defined($xmlData->{run}{instrument})) {
                $machineName = $xmlData->{run}{instrument};
            } else {
                warn "Failed to find instrument name in $paramsFile\n";
            }
        }
    }

    $machineName = 'unknown' unless $machineName;

    $$runFolderRef = $runFolder;
    $$machineNameRef = $machineName;
    $$tileAreaRef = $runParameters->{'TileArea'} if $runParameters and exists $runParameters->{'TileArea'};
}

sub calc_percent_by_field($;$;$;$)
{
    my $hash_ref = shift;
    my $match_key_fields_ref = shift;
    my $total_key = shift;
    my $separ_str = shift;

    my $divisor = $hash_ref->{$total_key};

    if (!defined($divisor)) {
        warn("Divisor ($total_key) not found when calculating percentage!");
        $divisor = 0;
    }

    $divisor /= 100.0;

    foreach my $key_str (keys %{$hash_ref}) {
        my @key_fields = split(/$separ_str/, $key_str);
        next if (!key_matched($match_key_fields_ref, \@key_fields));

        my $pc_key_str = join($separ_str, ($pc_prefix, $key_str));
        $hash_ref->{$pc_key_str} = (($divisor != 0)
                                    ? ($hash_ref->{$key_str} / $divisor)
                                    : 0);
    }
}

# Break down the unique paired alignments according to which of the individual
# read alignments were non-unique but `rescuable' - read1 only, read2 only
# or both reads.
#
# Calculate the percentages relative to the total in each rescuable category.

sub derive_uniqueness_rescue_stats($;$;$)
{
    my $lane_pair_ref = shift;
    my $num_paired_reads = shift;
    my $separ_str = shift;

    my $total = 0;
    my $total_divisor = 0;

    my $Read_str = 'Read';
    my $All_Reads_key = 'AllReads';
    my $Rescued_Unique_key = 'RescuedUnique';
    my $Rescuable_key = 'ManyAlignmentsFound';
    my $Unique_key = 'SingleAlignmentFound';
    my $unique_pair_align_key = 'UniquePairedAlignment';
    my $Total_str = 'Total';
    my $pc_prefix = 'pc';

    # The $bitmask is just a device to represent each read separately and
    # then both reads.
    for (my $bitmask = 1; $bitmask < (1 << $num_paired_reads); ++$bitmask) {
        my $all_rescued = 1;
        my @key_parts = ();
        my $new_key = '';

        for (my $read_num = 1; $read_num <= $num_paired_reads; ++$read_num) {
            my $rescued = ($bitmask & (1 << ($read_num - 1)));

            if ($rescued) {
                $new_key = $Read_str . $read_num;
            } else {
                $all_rescued = 0;
            }

            push(@key_parts, $Read_str . $read_num . ($rescued
                                                      ? $Rescuable_key
                                                      : $Unique_key));
        }

        if ($all_rescued) {
            $new_key = $All_Reads_key;
        }

        $new_key = join($separ_str, ($new_key, $Rescued_Unique_key));
        my $new_pc_key = join($separ_str, ($pc_prefix, $new_key));
        my $base_key = join($separ_str, @key_parts);
        my $val = $lane_pair_ref->{join($separ_str, ($base_key,
                                                     $unique_pair_align_key))};

        my $divisor = $lane_pair_ref->{$base_key};
        my $pc_val;

        if (defined($val) && defined($divisor)) {
            $pc_val = (($divisor > 0)
                       ? ($val / $divisor * 100.0)
                       : 0);
        } else {
            $pc_val = 0;
            $val = 0 if !defined($val);
        }

        $lane_pair_ref->{$new_key} = $val;
        $lane_pair_ref->{$new_pc_key} = $pc_val;
        $total += $val;
        $total_divisor += (defined($divisor) ? $divisor : 0);
    }

    my $total_key = join($separ_str, ($Total_str, $Rescued_Unique_key));
    $lane_pair_ref->{$total_key} = $total;
    $lane_pair_ref->{join($separ_str, ($pc_prefix, $total_key))}
    = (($total_divisor > 0)
       ? ($total / $total_divisor * 100.0)
       : 0);
}

sub init_lane_results($;$;$;$;$;$)
{
    my $lane_results_ref = shift;
    my $lane_parameters_ref = shift;
    my $lane_read_list_ref = shift;
    my $lane_num_tiles_ref = shift;
    my $maxLanes = shift;
    my $num_available_reads = shift;

  
    my %undef_hash;
    my $unknown = 'unknown';

    for (my $lane_num = 1; $lane_num <= $maxLanes; ++$lane_num) {
        # The number of reads for a lane is initialised to the number
        # available but may be reduced later according to analysis type.
        # $lane_read_list_ref->{$lane_num} = [1..$num_available_reads];
        $lane_num_tiles_ref->[$lane_num] = 0;

        for (my $read_num = 1; $read_num <= $num_available_reads;
             ++$read_num) {
            # Force $lane_results_ref->[$read_num][$lane] to exist before use.
            %{$lane_parameters_ref->[$lane_num]} = %undef_hash;
            %{$lane_results_ref->[$read_num][$lane_num]} = %undef_hash;

            my $parameters = $lane_parameters_ref->[$lane_num];
            $parameters->{tileCountRaw} = 0;
            #$parameters->{tileCountPF} = 0;
            # set to 1 if any _(re)score.txt successfully opened for 
            # tile in that lane.
            $parameters->{alignedFlag} = 0;

            foreach my $thisStat
                (qw/template sample type length lengthsList purity eland/) {
                $parameters->{$thisStat} = $unknown;
            }
        }
    }
}

sub read_config_file($;$;$;$;$;$;$;$;$)
{
    my $configFile = shift;
    my $lane_results_ref = shift;
    my $lane_parameters_ref = shift;
    my $maxNumLanes = shift;
    my $num_available_reads = shift;
    my $lane_read_list_ref = shift;
    my $max_reads_used_ref = shift;
    my $ref_read_num = shift;

    my $xmlData = readFileXML($configFile);

    if (!defined($xmlData)) {
        warn "Failed to parse $configFile as XML\n";
        exit(1);
    }

    my $chip_wide = $xmlData->{Defaults};
    die "ChipWideRunParameters missing in $configFile" unless $chip_wide;
    my $lane_specific = $xmlData->{Lane};
    die "LaneSpecificRunParameters missing in $configFile" unless $lane_specific;

    # get the analysis for each lane
    my %analyses;
    my %prefixes;
    foreach my $lane_prefix (sort keys %{$lane_specific->{ANALYSIS}}) {
        die "Invalid lane prefix in $configFile: $lane_prefix" unless $lane_prefix =~ /^[^_]_(\d+)$/;
        my $lane_num = $1;
        $prefixes{$lane_num} = $lane_prefix;
        my $analysis = $lane_specific->{ANALYSIS}->{$lane_prefix};
        $analyses{$lane_num} = $analysis if $analysis;
    }
    # get the read list for each lane
    foreach my $lane_num (sort keys %analyses) {
        my $analysis = $analyses{$lane_num};
        # Use the first non-empty use base masks
        my $numberOfReads = 1; # default is single ended
        $numberOfReads = 0 if $analysis eq "none";
        $numberOfReads = 2 if $analysis eq "eland_pair" or $analysis eq "sequence_pair"; 
        my $readsFound = 0;
        my @readList = ();
        my $currentRead = 0; # incremented in the test of the while
        # Arbitrarily limits the number of reads to 10
        while (scalar(@readList) <= $numberOfReads and ++$currentRead < 10)
        {
            my $useBasesKey = "USE_BASES$currentRead";
            my $useBases;
            $useBases = $chip_wide->{$useBasesKey} if exists $chip_wide->{$useBasesKey};
            $useBases = $lane_specific->{$useBasesKey}->{"s_$lane_num"} if exists $lane_specific->{$useBasesKey}->{"s_$lane_num"};
            push @readList, $currentRead if $useBases and $useBases =~ /[yY]/;
        }
        die "Found only ", scalar(@readList), " used reads un the USE_BASES masks. Expected $numberOfReads for analysis $analysis" if @readList < $numberOfReads;
        $lane_read_list_ref->{$lane_num} = \@readList;
    }

    my $unknown = 'unknown';

    my %configKeywords = ( "ANALYSIS"    => "type",
                           "GENOME_FILE" => "template",
                           "ELAND_GENOME" => "eland",
                           "READ_LENGTH" => "originalReadLength");

    foreach my $lane_num (sort keys %analyses) {
        #foreach my $read_num (@{$lane_read_list_ref->{$lane_num}}) {
            foreach my $keyword (keys %configKeywords) {
                my $value;
                my $chip_value = $chip_wide->{$keyword};
                my $lane_value = $lane_specific->{$keyword}->{$prefixes{$lane_num}};
                $value = $chip_value if $chip_value;
                $value = $lane_value if $lane_value;
                #$lane_results_ref->[$read_num][$lane_num]{$configKeywords{$keyword}} = $value;
                $lane_parameters_ref->[$lane_num]{$configKeywords{$keyword}} = $value;
            }
        #}
    }

    # no output for ANALYSIS none, use ELAND_GENOME as Eland target
    foreach my $lane_num (sort keys %analyses) {
        my $analysis_type = uc($analyses{$lane_num});
        #= uc($lane_results_ref->[$ref_read_num][$lane_num]{type});

        my $num_reads_used = @{$lane_read_list_ref->{$lane_num}};

        if ($num_reads_used > $$max_reads_used_ref) {
            $$max_reads_used_ref = $num_reads_used;
        }

        #for (my $read_num = 1; $read_num <= $num_available_reads;
        #     ++$read_num) {
        #foreach my $read_num (@{$lane_read_list_ref->{$lane_num}}) {
        #    my $lane_read_ref = $lane_results_ref->[$read_num][$lane_num];
        #    $lane_read_ref->{type} = $analysis_type;
        my $parameters = $lane_parameters_ref->[$lane_num];
        $parameters->{type} = $analysis_type;
            if ($analysis_type eq "NONE") {
                foreach my $value (qw/length template/) {
                    #$lane_read_ref->{$value} = $unknown;
                    $parameters->{$value} = $unknown;
                }
            } elsif (($analysis_type eq "SEQUENCE")
                     || ($analysis_type eq "SEQUENCE_PAIR")) {
                foreach my $value (qw/template/) {
                    #$lane_read_ref->{$value} = $unknown;
                    $parameters->{$value} = $unknown;
                }
            } elsif (($analysis_type eq "ELAND")
                     || ($analysis_type eq "ELAND_EXTENDED")
                     || ($analysis_type eq "ELAND_PAIR")
		     || ($analysis_type eq "ELAND_RNA")) {
                #$lane_read_ref->{template} = basename($lane_read_ref->{eland});
                $parameters->{template} = basename($parameters->{eland});
            }
        #}
    }

    # sets lengthsList and length values

    my $read_length_key = 'READ_LENGTH';

    # Initialise from chip-wide lengths in case lane-specific lengths are
    # not available.

    foreach my $lane_num (sort keys %analyses) {
        # Use num_read_used to determine whether to use READ_LENGTH
        # or to look for READ_LENGTHn.
        my $num_reads_used = @{$lane_read_list_ref->{$lane_num}};
        my @read_lengths;

        #for (my $read_num = 1; $read_num <= $num_reads_used; ++$read_num) {
        foreach my $read_num (@{$lane_read_list_ref->{$lane_num}}) {
            my $param_key = (($num_reads_used > 1)
                             ? $read_length_key . $read_num
                             : $read_length_key);
            my $read_length = $chip_wide->{$param_key};

            if (!defined($read_length)) {
                $read_length = 0;
            }

            push(@read_lengths, $read_length);
        }

        $lane_parameters_ref->[$lane_num]{lengthsList} = join(', ', @read_lengths);
    }

    my @lengths;

    foreach my $param_key (sort keys %{$lane_specific}) {
        next unless $param_key =~ m/^READ_LENGTH(.*)/;
        my $read_num = $1;
        my $read_lengths_ref = $lane_specific->{$param_key};
        next unless $read_lengths_ref;
        foreach my $lane_key (sort keys %{$read_lengths_ref}) {
            if ($lane_key =~ m/s_(.*)/) {
                my $lane_num = $1;
                my $analysis = uc($analyses{$lane_num});
                if ($analysis =~ /^(.*)_PAIR$/) {
                    push(@{$lengths[$lane_num]}, $read_lengths_ref->{$lane_key}) if $read_num;
                } else {
                    push(@{$lengths[$lane_num]}, $read_lengths_ref->{$lane_key}) unless $read_num;
                }
                $lane_parameters_ref->[$lane_num]{lengthsList} = join(', ', @{$lengths[$lane_num]});
            }
        } # for each lane
    } # for each read
}

sub process_pair_info_file($;$)
{
    my $pair_results_ref = shift;
    my $lane_num = shift;

    my $pair_info_fn = "s_".$lane_num."_pair.xml";
    my $xml_data = readFileXML($pair_info_fn);

    if (defined($xml_data)) {
        my %undef_hash;
        %{$pair_results_ref->{$lane_num}} = %undef_hash;
        my $lane_pair_ref = $pair_results_ref->{$lane_num};
        $lane_pair_ref->{$paired_not_done_key} = 0; # force to exist

        my $inserts_data = $xml_data->{$Insert_Size_str};

        if (defined($inserts_data)) {
            my @insert_stat_strs = map { $_->[0] } @insert_stat_display_map;

            foreach my $stat_str (@insert_stat_strs) {
                $lane_pair_ref->{"$Insert_Size_str:$stat_str"}
                = (defined($inserts_data->{$stat_str})
                   ? $inserts_data->{$stat_str}
                   : $not_avail_str);
            }
        }

        my $orient_data = $xml_data->{$Orientation_str};

        if (defined($orient_data)) {
            # Insert size stats

            my @orientation_keys = map { $_->[0] } @orient_display_map;

            foreach my $stat_key (@orientation_keys, $nominal_orient_key) {
                $lane_pair_ref->{"$Orientation_str:$stat_key"}
                = (defined($orient_data->{$stat_key})
                   ? $orient_data->{$stat_key}
                   : 0);
            }

            # Insert frequency stats

            my @insert_freq_keys = map { $_->[0] } @inserts_display_map;

            foreach my $stat_key (@insert_freq_keys) {
                # not provided so derived instead
                next if ($stat_key eq $num_pairs_ok_key);

                $lane_pair_ref->{"$Insert_Size_str:$stat_key"}
                = (defined($orient_data->{$stat_key})
                   ? $orient_data->{$stat_key}
                   : 0);

                my $pc_key = $insert_num_pc_key_map{$stat_key};

                $lane_pair_ref->{"$Insert_Size_str:$pc_key"}
                = (defined($orient_data->{$pc_key})
                   ? sprintf($insert_pc_format, $orient_data->{$pc_key})
                   : '0.0');
            }
        }

        my $reads_data = $xml_data->{'Reads'};

        if (defined($reads_data)) {
            my @align_outcome_keys
                = map { $_->[0] } @align_outcome_display_map;

            foreach my $read1_outcome_key (@align_outcome_keys) {
                my $key1_str = $Read_str . "1" . $read1_outcome_key;
                my $read1_data = $reads_data->{$key1_str};

                if (defined($read1_data)) {
                    foreach my $read2_outcome_key (@align_outcome_keys) {
                        my $key2_str = $Read_str . "2" . $read2_outcome_key;
                        my $read2_data = $read1_data->{$key2_str};

                        if (defined($read2_data)) {
                            # Find all leaf values at arbitrary levels of 
                            # nesting under this node and sum them to produce 
                            # the main entry of interest.
                            $lane_pair_ref->{"$key1_str:$key2_str"}
                            = sum_vals_under_node($read2_data);

                            # A few special cases at this level.
                            my @inner_keys = keys (%{$read2_data});

                            foreach my $inner_key (@inner_keys) {
                                if (($inner_key eq $unique_pair_align_key)
                                    || ($inner_key eq $multi_pair_align_key)) {
                                    $lane_pair_ref->
                                    {"$key1_str:$key2_str:$inner_key"}
                                    = sum_vals_under_node($read2_data->
                                                          {$inner_key});
                                } elsif ($inner_key eq $paired_not_done_key) {
                                    $lane_pair_ref->{$paired_not_done_key} = 1;
                                }
                            }
                        } else {
                            $lane_pair_ref->{"$key1_str:$key2_str"} = 0;
                        }
                    }
                } else {
                    foreach my $read2_outcome_key (@align_outcome_keys) {
                        my $key2_str = $Read_str . "2" . $read2_outcome_key;
                        $lane_pair_ref->{"$key1_str:$key2_str"} = 0;
                    }
                }
            }
        }
    }
}

#------------------------------------------------------------------------------
# Sum all leaf values under an XML node handling arbitrary levels of
# nesting using recursion.

# forward declaration apparently needed for recursion in Perl.

sub sum_vals_under_node($)
{
    my $xml_ref = shift;

    if (!ref($xml_ref)) {
        # There is no nesting; this is already a leaf value.
        return $xml_ref;
    }

    my $sum = 0;
    my @sub_tags = keys (%{$xml_ref});

    foreach my $sub_tag (@sub_tags) {
        my $val = $xml_ref->{$sub_tag};

        if (!ref($val)) {
            # a leaf value
            $sum += $val;
        } elsif (ref($val) eq 'HASH') {
            # more nesting
            $sum += sum_vals_under_node($val);
        } else {
            die "Found unexpected reference type in XML structure\n";
        }
    }

    return $sum;
}


# whether uniquely or otherwise) but its partner has not (whether because
# it was low quality - too many unknown bases - or simply could not be matched
# to the reference).
#
# The percentages calculated are relative to total number of aligned partner
# reads (not total number of read pairs).  

sub derive_mispairing_stats($;$)
{
    my $lane_pair_ref = shift;
    my $separ_str = shift;

    for (my $ref_read_num = 1; $ref_read_num <= 2; ++$ref_read_num) {
        my $lost_read_num = ($ref_read_num == 1) ? 2 : 1;

        my $count = 0;
        my $divisor = 0;

        foreach my $ref_key ($Unique_key, $Rescuable_key, $Repeat_key) {
            foreach my $lost_key($Not_Matched_key, $Low_Quality_key) {

                my $read1_key = $Read_str."1" . (($ref_read_num == 1)
                                                  ? $ref_key
                                                  : $lost_key);
                my $read2_key = $Read_str."2" . (($ref_read_num == 2)
                                                  ? $ref_key
                                                  : $lost_key);

                my $val = $lane_pair_ref->{join($separ_str,
                                                ($read1_key, $read2_key))};

                if (defined($val)) {
                    $count += $val;
                }
            }

            my @key_strs = ($Read_str . $ref_read_num . $ref_key, $Total_str);

            if ($ref_read_num != 1) {
                @key_strs = reverse @key_strs;
            }

            my $val = $lane_pair_ref->{join($separ_str, @key_strs)};
            if (defined($val)) {
                $divisor += $val;
            }
        }

        my $key = $Read_str . $lost_read_num . $Lost_str;
        my $pc_key = join($separ_str, ($pc_prefix, $key));

        $lane_pair_ref->{$key} = $count;
        $lane_pair_ref->{$pc_key} = (($divisor != 0)
                                     ? ($count / $divisor * 100.0)
                                     : 0.0);
    }
}

sub derive_insert_freq_stats($)
{
    my $lane_pair_ref = shift;

    my $nominal_orient =
        $lane_pair_ref->{"$Orientation_str:$nominal_orient_key"};

    my $ok_freq = (defined($nominal_orient)
                   ? $lane_pair_ref->{"$Orientation_str:$nominal_orient"}
                   : 0);

    foreach my $bad_freq_key ($num_orient_ok_small_insert_key,
                              $num_orient_ok_large_insert_key) {
        my $key = "$Insert_Size_str:$bad_freq_key";
        my $bad_freq = (defined($lane_pair_ref->{$key})
                        ? $lane_pair_ref->{$key}
                        : 0);
        $ok_freq -= $bad_freq;
    }

    $lane_pair_ref->{"$Insert_Size_str:$num_pairs_ok_key"} = $ok_freq;


    # Calc total num inserts and use as denominator in percentage.

    my $total_inserts = 0;
    my @orientation_keys = map { $_->[0] } @orient_display_map;

    foreach my $stat_key (@orientation_keys) {
        my $num_inserts = $lane_pair_ref->{"$Orientation_str:$stat_key"};
        $total_inserts += (defined($num_inserts) ? $num_inserts : 0);
    }

    my $pc_key = $insert_num_pc_key_map{$num_pairs_ok_key};

    $lane_pair_ref->{"$Insert_Size_str:$pc_key"}
    = sprintf($insert_pc_format, (($total_inserts != 0)
                                  ? (100.0 * $ok_freq / $total_inserts)
                                  : 0.0));
}

sub derive_pair_stats($)
{
    my $lane_pair_ref = shift;
    my $paired_align_done = !$lane_pair_ref->{$paired_not_done_key};

    # Individual Alignment Outcomes
    calc_total_by_field($lane_pair_ref,
                        ["^".$Read_str."1", "^".$Read_str."2"],
                        1, $pair_key_separ_str);
    calc_total_by_field($lane_pair_ref,
                        ["^".$Read_str."1", "^(".$Read_str."2|".$Total_str.")"],
                        0, $pair_key_separ_str);
   # calc_percent_by_field($lane_pair_ref,
   #                       ["^(${Read_str}1|${Total_str})",
   #                        "^(${Read_str}2|${Total_str})"],
   #                       "${Total_str}:${Total_str}",
   #                       $pair_key_separ_str);


    # Breakdown of UniquePairedAlignment, ManyPairedAlignments
    if ($paired_align_done) {
        foreach my $pair_align_key ($unique_pair_align_key,
                                    $multi_pair_align_key) {
            calc_total_by_field($lane_pair_ref,
                                ["^".$Read_str."1", "^".$Read_str."2",
                                 "^$pair_align_key\$"],
                                1, $pair_key_separ_str);
            calc_total_by_field($lane_pair_ref,
                                ["^".$Read_str."1",
                                 "^(".$Read_str."2|".$Total_str.")",
                                 "^$pair_align_key\$"],
                                0, $pair_key_separ_str);

           # calc_percent_by_field($lane_pair_ref,
           #                      ["^(${Read_str}1|${Total_str})",
           #                       "^(${Read_str}2|${Total_str})",
           #                       "^${pair_align_key}\$"],
           #                     "${Total_str}:${Total_str}:${pair_align_key}",
           #                      $pair_key_separ_str);
        }
    }
    # Uniqueness rescue
    if ($paired_align_done) {
        derive_uniqueness_rescue_stats($lane_pair_ref, $num_paired_reads, $pair_key_separ_str);
    }

    # Mispairing
    derive_mispairing_stats($lane_pair_ref, $pair_key_separ_str);

    # Orientations
    my @orientation_strs = map { $_->[0] } @orient_display_map;
    my $orient_regex = '^(' . join('|', @orientation_strs) . ')$';
    calc_total_by_field($lane_pair_ref,
                        ["^$Orientation_str\$", $orient_regex],
                        1, $pair_key_separ_str);
   # calc_percent_by_field($lane_pair_ref,
   #                       ["^${Orientation_str}\$", $orient_regex],
   #                       "${Orientation_str}:${Total_str}",
   #                       $pair_key_separ_str);

    # Insert size freqs
    if ($paired_align_done) {
        derive_insert_freq_stats($lane_pair_ref);
    }
}


sub read_tiles_file($;$)
{
    my $tiles_file = shift;
    my $tiles_ref = shift;

    if (!open(TILES, $tiles_file)) {
        warn "Failed to open tile list file `", $tiles_file, "'.\n";
        exit(0);
    }

    while (my $line = <TILES>) {
        next if ($line =~ m/^\#/);

        my @curr_tiles = split(" ", $line);
        push @$tiles_ref, @curr_tiles;
    }

    close(TILES);
}


sub process_score_file($;$;$;$;$;$;$;$;$;$;$;$)
{
    my $file_name = shift;
    my $infoContentKey = shift;
    my $errorKey = shift;
    my $pc_align_key = shift;
    my $av_align_score_key = shift;
    my $tile_num = shift;
    my $lane_num = shift;
    my $read_list_ref = shift;
    my $results_ref = shift;
    my $lane_results_ref = shift;
    my $lane_parameters_ref = shift;
    my $cluster_count = shift;

    my $info_content_line = undef;
    my $num_unique_aligns = undef;
    my $total_unique_align_score = undef;

    #for (my $read_num = 0; $read_num <= $num_reads; ++$read_num) { # 0 -> All
    my @read_list;
    push @read_list, 0 if 1 < @{$read_list_ref}; # 0 -> All
    push @read_list, @{$read_list_ref};

    foreach my $read_num (@read_list) {
        $results_ref->{$read_num}->{$lane_num}->{$tile_num} = {} unless exists $results_ref->{$read_num}->{$lane_num}->{$tile_num};
        my $tile_read_ref = $results_ref->{$read_num}->{$lane_num}->{$tile_num};
        # Should these all be set to undef rather than to 0,
        # i.e. should attempts to divide 0 by 0 be reported as 0 or unknown?
        $tile_read_ref->{$errorKey} = 0;
        $tile_read_ref->{$infoContentKey} = 0;
        $tile_read_ref->{$pc_align_key} = 0;
        $tile_read_ref->{$av_align_score_key} = 0;

    }

    my $read_num = $read_list_ref->[0];
    #my $highest_read_num = $read_num;
    my $tile_read_ref = $results_ref->{$read_num}->{$lane_num}->{$tile_num};

    open IN, $file_name or die "Failed to open $file_name: $!";

    #for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
    foreach my $read_num (@{$read_list_ref}) {
        $lane_parameters_ref->[$lane_num]{alignedFlag} = 1;
    }
    
    my $in_info_content = 0;
    
    while (<IN>) {
        #           print $_; die;
        
        if ($_ =~ /^\# Lane.*Read/) {
            if ($_ =~ /^\# Lane.*Read (\d+)/) {
                $read_num = $1;
                grep(/^$read_num$/, @{$read_list_ref}) or die ""
                . "Unexpected read: $read_num not in [" 
                . join(', ', @{$read_list_ref}) . "]";
                #if ($read_num > $num_reads) {
                #    die("Expected $num_reads reads but found "
                #        . "Read $read_num in $file_name\n");
                #}
                
                #if ($read_num > $highest_read_num) {
                #    $highest_read_num = $read_num;
                #}
            } else {
                $read_num = 0; # All Reads
            }
            
            $tile_read_ref = $results_ref->{$read_num}->{$lane_num}->{$tile_num};
        } elsif ($_ =~ m/ percent error rate/) {
            my @array = split(' ',$_);
            $tile_read_ref->{$errorKey} = $array[0];
        } elsif ($_ =~ /unique alignments : (\d+) \(total score (\d+)\)/) {
            $num_unique_aligns = $1;
            $total_unique_align_score = $2;
            
            if ($cluster_count != 0) {
                $tile_read_ref->{$pc_align_key}
                = ($num_unique_aligns / $cluster_count) * 100;
                
                $tile_read_ref->{$av_align_score_key}
                = $total_unique_align_score / $cluster_count;
                
            }
        } elsif ($_ =~ /^\~/) {
            $info_content_line = $_;
            $in_info_content = 1;
        } elsif ($in_info_content) {
            # Now at the line after the last info content for a stats set.
            # Process that last info content line.
            $in_info_content = 0;
            
            chomp($info_content_line);
            $info_content_line =~ s/^\~//;
            my @array = split('\t', $info_content_line);
            
            $tile_read_ref->{$infoContentKey} = $array[4] / $array[0]
            if ((@array == 7) && ($array[0] != 0));
        }
        
    } # while <IN>
    
    close(IN);

    #if ($highest_read_num < $num_reads) {
    #    die("Expected $num_reads but highest Read found in "
    #        . " $file_name was Read $highest_read_num\n");
    #}
}


sub getIndexRead($) {
   my $configFile = shift;

   my $xml = XML::Simple->new();
   my $xmlData = $xml->XMLin($configFile, suppressempty=>undef);

   if (!defined($xmlData)) {
     warn "Failed to parse $configFile as XML\n";
     exit(1);
   }
   
   my $indexRead = $xmlData->{'ChipWideRunParameters'}{'INDEXING_READ_POSITION'};
   return $indexRead;
}

#deletes indexing reads from hash
sub deleteIndexingRead($$$) {
   my $bxmlData_ref=shift;
   my $index=shift;
   my $max_reads_used=shift;

   #decrement since index read in config.xml is 1-offset
   $index-- if (defined($index) && ($index>0));

   if (defined($bxmlData_ref->{ExpandedLaneSummary}{Read}[$index])) {
      splice(@{$bxmlData_ref->{ExpandedLaneSummary}{Read}},$index,1);
      #renumbering reads
      foreach (@{$bxmlData_ref->{ExpandedLaneSummary}{Read}}) {
         $_->{readNumber}-- if ($_->{readNumber}>$index);
      }
   }

   if (defined($bxmlData_ref->{LaneResultsSummary}{Read}[$index])) {
      splice(@{$bxmlData_ref->{LaneResultsSummary}{Read}},$index,1);
      #renumbering reads
      foreach (@{$bxmlData_ref->{LaneResultsSummary}{Read}}) {
         $_->{readNumber}-- if ($_->{readNumber}>$index);
      }
   }

   for (my $lane_num=0; $lane_num<$maxNumLanes; ++$lane_num) {
      my $lane_read=$bxmlData_ref->{TileResultsByLane}{Lane};
      my $lane_ref=@{$lane_read}[$lane_num];

      if (defined(@{$lane_ref->{Read}}[$index])) {
         splice(@{$lane_ref->{Read}},$index,1);
         foreach (@{$lane_ref->{Read}}) {
            $_->{readNumber}-- if ($_->{readNumber}>$index);
         }
      }
   }
}

#------------------------------------------------------------------------------

sub make_lrt_key($;$;$)
{
    my ($lane_num, $read_num, $tile_num) = @_;
    return join(':', $lane_num, $read_num, $tile_num);
}

#------------------------------------------------------------------------------

sub get_cluster_counts($;\%)
{
    my ($bustard_xml_ref, $cluster_counts_ref) = @_;
    my %undef_hash;

    my $lane_hash_ref = $bustard_xml_ref->{TileResultsByLane}{Lane};
    die "Undefined TileResultsByLane in Bustard Summary" unless $lane_hash_ref;

    #for (my $lane_ind = 0; $lane_ind <scalar(@{$lane_arr_ref});
    #     ++$lane_ind) {
    foreach my $lane_num (sort keys %{$lane_hash_ref}) {
        my $lane_ref = $lane_hash_ref->{$lane_num};
        my $read_hash_ref = $lane_ref->{Read};
        next unless $read_hash_ref;
        foreach my $read_num (sort keys %{$read_hash_ref}) {
            my $read_ref = $read_hash_ref->{$read_num};
            my $tile_hash_ref = $read_ref->{Tile};
            next unless $tile_hash_ref;
            foreach my $tile_num (sort keys %{$tile_hash_ref}) {
                my $tile_ref = $tile_hash_ref->{$tile_num};
                my $count_hash = {};
                foreach my $count_key ($cluster_count_raw_key, $cluster_count_PF_key) {
                    $count_hash->{$count_key} = $tile_ref->{$count_key};
                }
                $cluster_counts_ref->{$lane_num}->{$read_num}->{$tile_num} = $count_hash;
            }
        }
    }
}

#------------------------------------------------------------------------------

sub get_tile_cluster_counts(\%;$;$;$;\$;\$) {
    my ($cluster_counts_ref, $lane_num, $read_num, $tile_num,
	$cluster_count_raw_ref, $cluster_count_PF_ref) = @_;

    $$cluster_count_raw_ref = 0;
    $$cluster_count_PF_ref = 0;
    return unless exists $cluster_counts_ref->{$lane_num};
    return unless exists $cluster_counts_ref->{$lane_num}->{$read_num};
    return unless exists $cluster_counts_ref->{$lane_num}->{$read_num}->{$tile_num};
    my $count_hash = $cluster_counts_ref->{$lane_num}->{$read_num}->{$tile_num};
    $$cluster_count_raw_ref = $count_hash->{$cluster_count_raw_key};
    $$cluster_count_PF_ref = $count_hash->{$cluster_count_PF_key};
}

#------------------------------------------------------------------------------

sub merge_lanes($$$$$$$$$) {
   my $laneResults_ref = shift;
   my $laneParameters_ref = shift;
   my $results_ref = shift;
   my $max_reads = shift;
   my $num_tiles = shift;
   my $lane_type_ref = shift;
   my $lane_num_reads_ref = shift;
   my $bxmlData_ref = shift;
   my $stats = shift;

   my $thisVal;
   my %formats = ('percentUniquelyAlignedPF' => '%3.2f',
                  'averageAlignScorePF' => '%3.2f',
                  'errorPF' => '%3.2f',
                  'errorRaw' => '%3.2f',
                  'errorPF' => '%3.2f',
                  'infoContentRaw' => '%d',
                  'infoContentPF' => '%d');
   
   #for (my $lane_num = 1; $lane_num <= $maxNumLanes; ++$lane_num) {
   #foreach my $lane_num (@all_lanes) {
   #   if (exists $lane_type_ref->{$lane_num}) {
   #      my $num_reads = $lane_num_reads_ref->[$lane_num];
   
   #      for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
   foreach my $read_num (sort keys %{$bxmlData_ref->{Read}}) {
       foreach my $lane_num (sort keys %{$bxmlData_ref->{Read}->{$read_num}->{Lane}}) {
           if (exists $lane_type_ref->{$lane_num}) {
               my $laneResults_read_ref = $laneResults_ref->[$read_num][$lane_num];
               my $parameters = $laneParameters_ref->[$lane_num];
               my $analysis_type = $parameters->{type};
               next if ($analysis_type eq "SEQUENCE") or ($analysis_type eq "SEQUENCE_PAIR");
               
               #my $lane_read = $bxmlData_ref->{Read};
               #my $read_ref = $lane_read->{$read_num};
               
               #my $lane_ref = $read_ref->{Lane};
               #my $lane_read_ref = $lane_ref->[$lane_num-1];
               my $lane_read_ref = $bxmlData_ref->{Read}->{$read_num}->{Lane}->{$lane_num};
               foreach my $tile_num (sort keys %{$results_ref->{$read_num}{$lane_num}}) {
                   my $tile_read_ref = $results_ref->{$read_num}{$lane_num}{$tile_num};
                   die "Undefined lane for s_${lane_num}_${tile_num}" unless defined $tile_read_ref->{lane};
                   die "Unexpected lane for s_${lane_num}_${tile_num}:" . $tile_read_ref->{lane} unless $tile_read_ref->{lane} == $lane_num;
                   foreach (@{$stats}) {
                       $thisVal=0;
                       if (defined($tile_read_ref->{$_})) {
                           $thisVal = $tile_read_ref->{$_};
                       }
                       $lane_read_ref->{$_}->{mean} += $thisVal;
                       $lane_read_ref->{$_}->{sumsq} += ($thisVal
                                                       * $thisVal);
                   }
               }

               foreach (@{$stats}) {
                   my $numSamples = $parameters->{tileCountPF};
                   
                   my $stat_ref = $lane_read_ref->{$_};
                   
                   if (defined($stat_ref)) {
                       
                       my $format = $formats{$_};
                       $stat_ref->{mean} = sprintf($format, $stat_ref->{mean});
                       $stat_ref->{sumsq} = sprintf($format, $stat_ref->{sumsq});
                       
                       if (defined($numSamples) && ($numSamples >= 2)) {
                           $stat_ref->{stdev} = ($stat_ref->{sumsq}
                                                 - (($stat_ref->{mean} * $stat_ref->{mean})
                                                    / $numSamples));
                           
                           $stat_ref->{stdev} /= ($numSamples - 1);
                           if ($stat_ref->{stdev} < 0) {
                               $stat_ref->{stdev} = 0;
                           }
                           
                           $stat_ref->{stdev} = sprintf($format, sqrt($stat_ref->{stdev}));
                       } else {
                           $stat_ref->{stdev} = 0;
                           #$stat_ref->{stdev} = undef;
                       }
                       
                       if (defined($numSamples) && ($numSamples > 0)) {
                           $stat_ref->{mean} /= $numSamples;
                           $stat_ref->{mean} = sprintf($format, $stat_ref->{mean});
                       } else {
                           $stat_ref->{mean} = undef;
                       }
                   }
               }
           }#for read
       }#if
   }#for lane
}

# FIXME : would seem cleaner to process results of using tileNumber as KeyAttr

sub merge_tiles($$$$$$$) {
   my $laneResults_ref = shift;
   my $results_ref = shift;
   my $max_reads = shift;
   my $num_aligned_tiles = shift;
   my $lane_type_ref = shift;
   my $lane_read_list_ref = shift;
   my $bxmlData_ref = shift;

   for (my $lane_num = 1; $lane_num <= $maxNumLanes; ++$lane_num) {
      if (exists $lane_type_ref->{$lane_num}) {
          #my $num_reads = @{$lane_read_list_ref->{$lane_num}};

         #for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
          foreach my $read_num (@{$lane_read_list_ref->{$lane_num}}) {

             #my $lane_read = $bxmlData_ref->{Lane};
             #my $lane_ref = $lane_read->{$lane_num};

             #my $read_ref = $lane_ref->{Read};

             #my $lane_read_ref = $read_ref->[$read_num-1];

             my $lane_read_ref = $bxmlData_ref->{Lane}->{$lane_num}->{Read}->{$read_num};
             my $lane_num_basecalled_tiles = scalar(keys %{$lane_read_ref->{Tile}});

             my $laneResults_read_ref
               = $laneResults_ref->[$read_num][$lane_num];

             my $basecalled_tile_index = 0;

             #my $num_aligned_tiles_for_read
             #  = scalar(@{$results_ref->[$read_num]});

             #for (my $aligned_tile_ind = 0;
             #     $aligned_tile_ind < $num_aligned_tiles_for_read;
             #     ++$aligned_tile_ind) {
             foreach my $tile_num (sort keys %{$results_ref->{$read_num}->{$lane_num}}) {
                 my $aligned_tile_read_ref
                   = $results_ref->{$read_num}->{$lane_num}->{$tile_num};
                 my $aligned_lane_num = $aligned_tile_read_ref->{lane};

                 if (!defined($aligned_lane_num)) {
                     warn("'lane' not specified for aligned tile");
                     next;
                 }

                 #next if ($aligned_tile_read_ref->{lane} != $lane_num);
                 #if ($aligned_lane_num != $lane_num) {
                 #    $basecalled_tile_index = 0;
                 #    next;
                 #}

                 #my $aligned_tile_num = $aligned_tile_read_ref->{tileNumber};

                 #if (!defined($aligned_tile_num)) {
                 #    warn("'tileNumber' not specified for aligned tile");
                 #    next;
                 #}

                 my $lane_read_tile_ref
                   = $lane_read_ref->{Tile}->{$tile_num};
                 #my $basecalled_tile_num = $lane_read_tile_ref->{tileNumber};

                 #if (!defined($basecalled_tile_num)) {
                 #    warn("'tileNumber' not specified for basecalled tile");
                 #    $basecalled_tile_num = 0;
                 #}

                 # If tiles have been selected at the alignment stage,
                 # may need to skip past tiles filtered out.
                 #my $tile_num_matched = 1;

                 #while ($basecalled_tile_num != $aligned_tile_num) {
                 #    ++$basecalled_tile_index;

                 #    if ($basecalled_tile_index
                 #        >= $lane_num_basecalled_tiles) {
                 #        warn("Ran out of basecalled tiles looking for tile "
                 #             . "num match for aligned tile "
                 #             . $aligned_tile_num
                 #             . " in lane $aligned_lane_num");
                 #        $tile_num_matched = 0;
                 #        last;
                 #    }

                 #    $lane_read_tile_ref
                 #      = $lane_read_ref->{Tile}[$basecalled_tile_index];
                 #    $basecalled_tile_num = $lane_read_tile_ref->{tileNumber};

                 #    if ($basecalled_tile_num > $aligned_tile_num) {
                 #        warn("Failed to find tile num match for aligned tile "
                 #             . $aligned_tile_num
                 #             . " in lane $aligned_lane_num");
                 #        $tile_num_matched = 0;
                 #        last;
                 #    }
                 #}

                 #next if (!$tile_num_matched);

                 $lane_read_tile_ref->{averageAlignScorePF}
                   = sprintf("%3.2f",
                             $aligned_tile_read_ref->{averageAlignScorePF})
                     if $aligned_tile_read_ref->{averageAlignScorePF};

                 $lane_read_tile_ref->{percentUniquelyAlignedPF}
                   = sprintf("%3.2f",
                             $aligned_tile_read_ref->{percentUniquelyAlignedPF})
                     if $aligned_tile_read_ref->{percentUniquelyAlignedPF};

                 $lane_read_tile_ref->{errorPF}
                   = sprintf("%3.2f",
                             $aligned_tile_read_ref->{errorPF})
                     if $aligned_tile_read_ref->{errorPF};

                 ++$basecalled_tile_index;
             }
         }
     }
  }
}


1;    # says use was ok
