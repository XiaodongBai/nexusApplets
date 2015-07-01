#!/usr/bin/env perl

# Copyright (c) 2008 Illumina Cambridge Ltd
# Author M. Zerara
# Several functions are taken from jerboa.pl
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).

# produceIntensityStats.pl 
# Extracts intensity data from _all.txt files 
# and generates a xml file BustardSummary.xml


use warnings;
use strict;
use POSIX qw(strftime);

use XML::Simple;
use File::Basename;
use File::Spec;
use Carp;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Common::Utils;
use Gerald::Jerboa;
use Casava::Demultiplex::SampleSheet;


# some useful constants
my $unknown='unknown';
my $maxNumLanes = 8;
my $numChannels=4;
my $ref_read_num = 1; # Use this to index values that are read-independent.
my $allSuffix='_all.txt';
my $timeStampFormat = "%Y-%m-%d %H:%M:%S";

my $not_avail_str = 'N/A';

# Intensities less than minIntensity get ignored in dropoff calculations
my $minIntensity=5;

#------------------------------------------------------------------------------

sub init_lane_stats($$$$$)
{
    my $lane_results_ref = shift;
    my $lane_num_reads_ref = shift;
    my $lane_num_tiles_ref = shift;
    my $maxNumLanes = shift;
    my $num_available_reads = shift;


    my %undef_hash;

    for (my $lane_num = 1; $lane_num <= $maxNumLanes; ++$lane_num) {
        # The number of reads for a lane is initialised to the number
        # available but may be reduced later according to analysis type.
        $lane_num_reads_ref->[$lane_num] = $num_available_reads;
        $lane_num_tiles_ref->[$lane_num] = 0;

        for (my $read_num = 1; $read_num <= $num_available_reads;
             ++$read_num) {
            # Force $lane_results_ref->[$read_num][$lane] to exist before use.
            %{$lane_results_ref->[$read_num][$lane_num]} = %undef_hash;
            my $lane_read_ref = $lane_results_ref->[$read_num][$lane_num];

            $lane_read_ref->{phasingApplied} = undef;
            $lane_read_ref->{prephasingApplied} = undef;
            $lane_read_ref->{phasingmeasured} = undef;
            $lane_read_ref->{prephasingmeasured} = undef;
            $lane_read_ref->{tileCountRaw} = 0;
            $lane_read_ref->{tileCountPF} = 0;

            # set to 1 if any _(re)score.txt successfully opened for 
            # tile in that lane.
            $lane_read_ref->{alignedFlag} = 0;

            $lane_read_ref->{type} = $unknown;

        }
    }
}

#------------------------------------------------------------------------------

sub logLinearFit
{
    my ($x, $y)=@_;
    die unless (scalar(@{$x})==scalar(@{$y}));
    my ($x_sum, $y_sum, $xx_sum, $xy_sum, $i, $ylog, $a, $b)
    =(0,0,0,0,0,0,0,0);
    my $numDataPoints=0;

    for ($i=0; $i<@{$x}; $i++)
    {
#        print "$i ".$x->[$i]." ".$y->[$i]."\n";
        next unless ($y->[$i]>=$minIntensity);
        ++$numDataPoints;

        $ylog=log($y->[$i]);
        $x_sum+=$x->[$i];
        $y_sum+=$ylog;
        $xx_sum+=($x->[$i]*$x->[$i]);
        $xy_sum+=($x->[$i]*$ylog);
    } # for
#    print "$i\n";

    unless ($numDataPoints>=2)
    {
        warn "Not enough valid data points to do fit";
        return 100.0;
    } # unless

    # Work out gradient of best fit line
    $a=($xy_sum-($x_sum*$y_sum/$i))/($xx_sum-($x_sum*$x_sum/$i));

    # Convert to a percent loss per cycle
    return (100*(1-exp($a)));
} # sub logLinearFit

#------------------------------------------------------------------------------

sub process_all_file($;$;$;$;$;\@;$)
{
    my $file_name = shift;
    my $tile_ind = shift;
    my $curr_tile_lane = shift;
    my $lane_results_ref = shift;
    my $results_ref = shift;
    my $read_startcycles = shift;
    my $num_reads = shift;

    my @array;

    my $thisIntensity;
    my @cyc2to10num;
    my @cyc2to10inten;
    my @cyc10to20num;
    my @cyc10to20inten;

    my ($combined_cycle_num, $read_num, $cycle_num);

    # Initialise tile statistics.
    # (These are also the final values if the tile _all file cannot be opened.)
    for ($read_num = 1; $read_num <= $num_reads; ++$read_num) {
        my $tile_read_ref = $results_ref->[$read_num][$tile_ind];
        $tile_read_ref->{oneSig} = 0;
        $tile_read_ref->{signalAverage1to3} = 0;
        $tile_read_ref->{signalAverage2to4} = 0;
        $tile_read_ref->{signal20AsPctOf1} = 0;
        $tile_read_ref->{signalLoss2to10} = undef;
        $tile_read_ref->{signalLoss10to20} = undef;
    }

    my $raw_cluster_count = 0;
    my $PF_cluster_count = 0;

    # Read _all file (one line per combined cycle after the comments)
    # and accumulate statistics.

    if (defined(open(IN, $file_name))) {
        @cyc2to10num = ();
        @cyc2to10inten = ();
        @cyc10to20num = ();
        @cyc10to20inten = ();

        while (<IN>) {
            if ($_ =~ m/^\#/) {
                if ($_
                    =~ m/^\# Clusters: Filtered ([0-9]+) Original ([0-9]+)/) {

                    $PF_cluster_count = $1;
                    $raw_cluster_count = $2;
                }

                next;
            }

            @array = split('\t', $_);
            $combined_cycle_num = $array[1];

            my $curr_cycle = $array[1];
            $cycle_num++;

            for (my $sc=1; $sc<=@{$read_startcycles}; $sc++) {
               if ($combined_cycle_num>=@{$read_startcycles}[$sc-1]) {
                  $read_num=$sc;
                  $cycle_num = $combined_cycle_num - @{$read_startcycles}[$sc-1] + 1;
               } 
            }

            my $tile_read_ref = $results_ref->[$read_num][$tile_ind];

            $thisIntensity = (($array[2] + $array[3] + $array[4] + $array[5])
                              / $numChannels);

            if (($cycle_num >= 2) && ($cycle_num <= 4)) {
                $tile_read_ref->{signalAverage2to4} += $thisIntensity;
            }

            if (($cycle_num >= 2) && ($cycle_num <= 10)) {
                push @cyc2to10num, $cycle_num;
                push @cyc2to10inten, $thisIntensity;
            }

            if (($cycle_num >= 10) && ($cycle_num <= 20)) {
                push @cyc10to20num, $cycle_num;
                push @cyc10to20inten, $thisIntensity;
            }

            if ($cycle_num =~ /^[123]$/) {
                $tile_read_ref->{signalAverage1to3} += $thisIntensity;
                if ($cycle_num eq '1') {
                    $tile_read_ref->{oneSig} += $thisIntensity;
                }
            } elsif ($cycle_num eq '20') {
                $tile_read_ref->{signal20AsPctOf1} += $thisIntensity;
                # Skip remainder of file once all needed intensities read.
                # NB assumes cycle 20 data always comes after the other cycles
                #    of interest and (for multiple reads) that cycle 20 data
                #    for the last read always comes after that for the others. 
                last if ($read_num == $num_reads);
            }

        } # while loop over _all.txt lines corresponding to combined cycles

        # Finish calculating multi-cycle statistics.

        for ($read_num = 1; $read_num <= $num_reads; ++$read_num) {
            my $tile_read_ref = $results_ref->[$read_num][$tile_ind];

            if ($tile_read_ref->{oneSig} != 0) {
                $tile_read_ref->{signal20AsPctOf1}
                = ($tile_read_ref->{signal20AsPctOf1}
                   / $tile_read_ref->{oneSig}) * 100;
            } else {
                $tile_read_ref->{signal20AsPctOf1} = undef;
            }

            if (defined($tile_read_ref->{signalAverage1to3})) {
                $tile_read_ref->{signalAverage1to3} /= 3;
            }

            if (defined($tile_read_ref->{signalAverage2to4})) {
                $tile_read_ref->{signalAverage2to4} /= 3;
            }

            $tile_read_ref->{signalLoss2to10}
            = logLinearFit(\@cyc2to10num, \@cyc2to10inten);
            $tile_read_ref->{signalLoss10to20}
            = logLinearFit(\@cyc10to20num, \@cyc10to20inten);
        }
    } else {
        warn("`$file_name' not found.");
    }

    # Store cluster count statistics (same values for all reads).

    for ($read_num = 1; $read_num <= $num_reads; ++$read_num) {
        my $tile_read_ref = $results_ref->[$read_num][$tile_ind];

        $tile_read_ref->{clusterCountRaw} = $raw_cluster_count;
        $tile_read_ref->{clusterCountPF} = $PF_cluster_count;

        if ($PF_cluster_count != 0) {
            $lane_results_ref->[$read_num][$curr_tile_lane]{tileCountPF}++;
        }

        $tile_read_ref->{percentClustersPF}
        = (($raw_cluster_count != 0)
           ? ($PF_cluster_count / $raw_cluster_count * 100)
           : 0);
    }
}

#------------------------------------------------------------------------------

sub get_stat_names_by_type($$$)
{
    my $stats_raw_ref = shift;
    my $stats_PF_ref = shift;
    my $total_stats_ref = shift;

    my @prealign_stats_raw = (qw/clusterCountRaw/);

    my @prealign_stats_PF = qw/clusterCountPF percentClustersPF oneSig signal20AsPctOf1 
                               signalAverage1to3 signalAverage2to4 signalLoss2to10 signalLoss10to20/;

        @$stats_raw_ref = @prealign_stats_raw;
        @$stats_PF_ref = @prealign_stats_PF;

    @$total_stats_ref = qw/yield clusterCountRaw clusterCountPF/; #DAP-373: this should not be used
}

#------------------------------------------------------------------------------

sub calc_lane_results($$$$$$)
{
    my $laneResultsRef = shift;
    my $resultsRef = shift;
    my $maxNumLanes = shift;
    my $lane_num_reads_ref = shift;
    my $lane_type_ref = shift;
    my $numTiles = shift;

    my (@stats, @stats_raw, @stats_PF, @total_stats);

    my $thisVal;
    my $numSamples;

    for (my $lane_num = 1; $lane_num <= $maxNumLanes; $lane_num++) {
        my $num_reads = $lane_num_reads_ref->[$lane_num];

        for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
            my $lane_read_ref = $laneResultsRef->[$read_num][$lane_num]; 

            next unless defined($lane_read_ref);

            get_stat_names_by_type(\@stats_raw, \@stats_PF, \@total_stats);

            @stats = (@stats_raw, @stats_PF);

            if (exists $lane_type_ref->{$lane_num}) {
                for (my $tileInd = 0; $tileInd < $numTiles; $tileInd++) {
                    my $tile_read_ref = $resultsRef->[$read_num][$tileInd];

                    next if (!defined($tile_read_ref->{lane}));
                    next if ($tile_read_ref->{lane} != $lane_num);

                    if ($tile_read_ref->{percentClustersPF} != 0) {
                        foreach my $thisStat (@stats) {
                            $thisVal = 0;

                            if (defined($tile_read_ref->{$thisStat})) {
                                $thisVal = $tile_read_ref->{$thisStat};
                            }

                            $lane_read_ref->{$thisStat}{mean} += $thisVal;
                            $lane_read_ref->{$thisStat}{sumsq} += ($thisVal
                                                                   * $thisVal);
                        }
                    }
                } # for tileInd

                foreach my $thisStat (@stats) {
                    $numSamples = $lane_read_ref->{tileCountPF};

                    my $stat_ref = $lane_read_ref->{$thisStat};

                    if (defined($stat_ref)) {
                        if ($numSamples >= 2) {
                            $stat_ref->{stdev}
                            = ($stat_ref->{sumsq}
                               - (($stat_ref->{mean} * $stat_ref->{mean})
                                  / $numSamples));

                            $stat_ref->{stdev} /= ($numSamples - 1);
                           
                            if ($stat_ref->{stdev} < 0) {
                                $stat_ref->{stdev}=0;
                            }

                            $stat_ref->{stdev} = sqrt($stat_ref->{stdev});
                        } else {
                            $stat_ref->{stdev} = undef;
                        }

                        if ($numSamples > 0) {
                            $stat_ref->{mean} /= $numSamples;
                        } else {
                            $stat_ref->{mean} = undef;
                        }
                    }
                } # foreach (@stats)

                # Yield (bases) stats are obtained from clusterCountPF by 
                # simply multiplying by tileCountPF and by read length.
                if (defined($lane_read_ref->{clusterCountPF})
                    && defined($lane_read_ref->{length})) {
                    $lane_read_ref->{yield}{total}
                    = ($lane_read_ref->{clusterCountPF}{mean}
                       * $lane_read_ref->{tileCountPF}
                       * $lane_read_ref->{length}) if ($lane_read_ref->{length} ne $unknown);
                } else {
                    $lane_read_ref->{yield}{total} = 0;
                }

            } # if
        } # for read_num
    } #  for lane_num
}

#------------------------------------------------------------------------------

sub read_chip_file($$$)
{
    my $appName = shift;
    my $chipFile = shift;
    my $chipNameRef = shift;

    my $xmlData;

    if (defined($xmlData = readFileXML($chipFile))) {

        if (defined ($xmlData->{ChipID})) {
            $$chipNameRef = $xmlData->{ChipID};
        } else {
            warn "Unable to find field ChipID in file $chipFile";
        }
    }
}

#------------------------------------------------------------------------------

sub calc_chip_results($$$$$$)
{
    my $lane_results_ref = shift;
    my $num_lanes = shift;
    my $lane_num_reads_ref = shift;
    my $num_available_reads = shift;
    my $chip_stats_by_read_ref = shift;
    my $chip_totals_ref = shift;

    my (@stats_raw, @stats_PF, @total_stats);

    # Alignment stats are filtered out below if necessary (by lane),
    # so ask for all @stats_raw and @stats_PF by specifying $is_aligned as 1.
    get_stat_names_by_type(\@stats_raw, \@stats_PF, \@total_stats);

    my @stats = (@stats_raw, @stats_PF);

    foreach my $curr_stat (@stats) {
        # Does this stat require alignment to be meaningful?

        # Is this stat purity filtered?
        my $stat_is_PF = 0;

        foreach ('percentClustersPF', @stats_raw) {
            if ($curr_stat eq $_) {
                $stat_is_PF = 1;
                last;
            }
        }

        my @read_total_tiles_to_count = ();
        my @read_chip_sums = ();

        for (my $read_num = 1; $read_num <= $num_available_reads; ++$read_num) {
            $read_chip_sums[$read_num] = 0;
            $read_total_tiles_to_count[$read_num] = 0;
        }

        for (my $lane_ind = 1; $lane_ind <= $num_lanes; ++$lane_ind) {
            my $num_reads = $lane_num_reads_ref->[$lane_ind];

            for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
                my $lane_read_ref = $lane_results_ref->[$read_num][$lane_ind];

                # Some statistics are relevant only if alignment has been 
                # performed - assumed to correspond to the presence or absence
                # of _realign files for tiles in this lane and read.
                my $is_aligned = ($lane_read_ref->{alignedFlag} != 0);

                next if (!defined($lane_read_ref->{$curr_stat}));

                my $num_tiles_to_count = $lane_read_ref->{'tileCountPF'};

                next if ($num_tiles_to_count == 0);

                # The lane means are over relevant number of tiles,
                # so have to be multipled by those to get the original sums.
                $read_chip_sums[$read_num]
                    += ($lane_read_ref->{$curr_stat}{mean}
                        * $num_tiles_to_count);

                my $cluster_stat
                    = ($stat_is_PF ? 'clusterCountPF' : 'clusterCountRaw');
                my $num_clusters_to_count
                    = ($lane_read_ref->{$cluster_stat}{mean}
                       * $num_tiles_to_count);

                if (!defined($lane_read_ref->{$cluster_stat}{total})) {
                    $lane_read_ref->{$cluster_stat}{total}
                    = $num_clusters_to_count;
                }

                $read_total_tiles_to_count[$read_num] += $num_tiles_to_count;
            } # read_num
        } # lane_ind

        for (my $read_num = 1; $read_num <= $num_available_reads; ++$read_num) {
            $chip_stats_by_read_ref->[$read_num]{$curr_stat}{by_tile}
            = (($read_total_tiles_to_count[$read_num] != 0)
               ? ($read_chip_sums[$read_num]
                  / $read_total_tiles_to_count[$read_num])
               : undef);
        }
    } # stat

    # For some stats, totals across lanes and reads are meaningful.
    # foreach my $curr_stat (@total_stats) { # DAP-373
    foreach my $curr_stat ('clusterCountPF', 'clusterCountRaw') {
        my $chip_sum = 0;

        for (my $lane_ind = 1; $lane_ind <= $num_lanes; $lane_ind++) {
            my $num_reads = $lane_num_reads_ref->[$lane_ind];

            for (my $read_num = 1; $read_num <= 1; ++$read_num) {
                my $lane_read_ref = $lane_results_ref->[$read_num][$lane_ind];

                # Need this or the next check will cause it to be defined.
                next if (!defined($lane_read_ref->{$curr_stat}));

                if (defined($lane_read_ref->{$curr_stat}{total})) {
                    $chip_sum += $lane_read_ref->{$curr_stat}{total};
                }
            }
        }

        $chip_totals_ref->{$curr_stat} = $chip_sum;
    }
    foreach my $curr_stat ('yield') {
        my $chip_sum = 0;

        for (my $lane_ind = 1; $lane_ind <= $num_lanes; $lane_ind++) {
            my $num_reads = $lane_num_reads_ref->[$lane_ind];

            for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
                my $lane_read_ref = $lane_results_ref->[$read_num][$lane_ind];

                # Need this or the next check will cause it to be defined.
                next if (!defined($lane_read_ref->{$curr_stat}));

                if (defined($lane_read_ref->{$curr_stat}{total})) {
                    $chip_sum += $lane_read_ref->{$curr_stat}{total};
                }
            }
        }

        $chip_totals_ref->{$curr_stat} = $chip_sum;
    }
}

#------------------------------------------------------------------------------

# Read bustard path from xml file rather than from argument list
my $bustardPath = File::Spec->rel2abs($ARGV[0]);
my $configXMLFile = getConfigurationFilePath($bustardPath);

my $firecrestPath = File::Spec->catdir($bustardPath, File::Spec->updir());
my $signal_means_path = File::Spec->catdir($bustardPath, "SignalMeans");

my $runFolderPath = getRunFolderPath($bustardPath);

my $bustard_dir = basename($bustardPath);
my $param_is_legacy = param_file_is_legacy_style($bustard_dir);

my @read_start_cycles = ();

my ($runFolder, $machineName, $tileArea);
my $appName = $0;
extract_run_info($appName, $bustardPath, \$runFolder, \$machineName, \$bustard_dir, \$tileArea);

my %phasing_options;

if (!$param_is_legacy) {
    # FIXME : Multiple reads of $configXMLFile.
    if (!getPhasingOptions($configXMLFile, $bustard_dir,
                           %phasing_options)) {
        warn "Failed to get phasing options\n";
        exit(1);
    }

    if (!getReadStartCycles($configXMLFile, $bustard_dir,
                            @read_start_cycles)) {
        warn "Failed to extract read start cycles from $configXMLFile\n";
        exit(1);
    }
} else {
    # legacy
    @read_start_cycles = (1);
}

my $num_available_reads = @read_start_cycles;

my @laneResults;
my @lane_num_reads;
my @lane_num_tiles;

init_lane_stats(\@laneResults, \@lane_num_reads, \@lane_num_tiles,
                  $Gerald::Jerboa::maxNumLanes, $num_available_reads);

store_measured_phasings(\@laneResults, $bustardPath,
                        $Gerald::Jerboa::maxNumLanes, $num_available_reads,
                        \%phasing_options, $param_is_legacy);

store_applied_phasings($appName, \@laneResults,
                       $Gerald::Jerboa::maxNumLanes, $num_available_reads,
                       $bustardPath, $configXMLFile,
                       \%phasing_options, $param_is_legacy);


my $sampleFile = "$runFolderPath/samples.xml";
read_sample_file($appName, $sampleFile, \@laneResults, $num_available_reads);

#store the read lengths
my @read_lengths = ();
getReadLengths($configXMLFile, $bustard_dir, @read_lengths);

die "produceIntensityStats.pl: the numbers of reads and start cycles differ!" 
    unless (@read_lengths == @read_start_cycles);

for my $read ( 0 .. $#laneResults ) {
   if (defined($laneResults[$read])) {
      my $lanes = $laneResults[$read];
      for my $lane ( 0 .. $#{$lanes} ) {
          $laneResults[$read][$lane]{length}=$read_lengths[$read-1];
      }
   }
}

my @tiles = ();
my $tilesFile = File::Spec->catdir($bustardPath, "tiles.txt");
read_tiles_file($tilesFile, \@tiles);

if (@tiles == 0)
{
    warn "No tilenames found in `$tilesFile'\n";
    exit(1);
}

# Ensure sensible (or at least consistent) file ordering
@tiles = sort(@tiles);

my $tileInd=0;
my %undef_hash;

my @results;
my %lane_type;


foreach my $tile_name (@tiles) {

    ## Determine Lane and Tile
    my @array = split("_", $tile_name);
    my $curr_tile_lane = $array[1];
    my $tile_id = $array[2];
    my $num_reads = $lane_num_reads[$curr_tile_lane];

    my $lane_type_ref = \%lane_type;
$num_reads = $num_available_reads; #take all reads, no masking

    my $lo_read_num = (($num_reads > 1) ? 0 : 1);

    # Initialise tile values for the lane-specific number of reads.
    for (my $read_num = $lo_read_num; $read_num <= $num_reads; ++$read_num) {    

        # Force $results[$read_num][$tileInd] to exist before referencing it.
        %{$results[$read_num][$tileInd]} = %undef_hash;
        my $tile_read_ref = $results[$read_num][$tileInd];

        $tile_read_ref->{lane} = $curr_tile_lane;
        $tile_read_ref->{tile} = $tile_id;

        my $lane_read_ref = $laneResults[$read_num][$curr_tile_lane];

        $lane_read_ref->{tileCountRaw}++;

        # Set some defaults - mostly for Eland.

        $tile_read_ref->{template} = $lane_read_ref->{template};
        $tile_read_ref->{length} = $lane_read_ref->{length};
        $tile_read_ref->{type} = $lane_read_ref->{type};
        $tile_read_ref->{clusterCountRaw} = 0;
        $tile_read_ref->{errorPF} = 0;

        # Set lane type from single reads (i.e. read_num > 0) if available.
        if (($read_num > 0) && defined($tile_read_ref->{type})) {
            if ($tile_read_ref->{type} ne "NONE") {
                $lane_type_ref->{$curr_tile_lane} = $tile_read_ref->{type};
            }
        }
 
    }

    # Process _all.txt file for current tile.
    process_all_file("$signal_means_path/$tile_name$allSuffix",
                     $tileInd, $curr_tile_lane,
                     \@laneResults, \@results,
                     @read_start_cycles, $num_reads);

    $tileInd++;
}

#------------------------------------------------------------------------------

sub print_xml_chip($$$$$$) {
    my $machine_name = shift;
    my $run_folder = shift;
    my $chip_name = shift;
    my $tileArea = shift;
    my $chip_totals_ref = shift;
    my $xml_results_ref = shift;

    my $results = {};
    $results->{'Machine'} = $machine_name;
    $results->{'RunFolder'} = $run_folder;
    $results->{'ChipID'} = $chip_name;
    $results->{'TileArea'} = $tileArea if $tileArea;

    $xml_results_ref->{'ChipSummary'} = $results;

    $xml_results_ref->{'ChipResultsSummary'} = $chip_totals_ref;

}

#------------------------------------------------------------------------------

sub print_xml_expanded_lane($$$$$$) {
    my $laneResults_ref = shift;
    my $max_reads = shift;
    my $maxNumLanes = shift;
    my $lane_type_ref = shift;
    my $lane_num_reads_ref = shift;
    my $xml_results_ref = shift;

    my $results = {};
    for (my $read_num=1; $read_num <= $max_reads; $read_num++) {
       my $arrayhash = {};
       for (my $lane = 1; $lane <= $maxNumLanes; $lane++) {

          next if ($lane_num_reads_ref->[$lane] < $read_num);

          my $lanehash = {};
          if (exists $lane_type_ref->{$lane}) {
             my $lane_read_ref = $laneResults_ref->[$read_num][$lane];
             die "Undefined lane read at: read_num = $read_num, lane = $lane"
                 unless defined $lane_read_ref;
             unless (exists $lane_read_ref->{clusterCountRaw})
             {
                 warn "Fixing corrupted lane read at: read_num = $read_num, lane = $lane";
                 $lane_read_ref->{clusterCountRaw}={};
             }

             #stdev defined only when clusterCountPF>=2
             $lane_read_ref->{clusterCountRaw}{stdev}=0  unless defined $lane_read_ref->{clusterCountRaw}{stdev};
             $lane_read_ref->{clusterCountRaw}{mean}=0   unless defined $lane_read_ref->{clusterCountRaw}{mean};
             $lane_read_ref->{clusterCountRaw}{sumsq}=0  unless defined $lane_read_ref->{clusterCountRaw}{sumsq};

             $lanehash->{clusterCountRaw}{stdev}=sprintf("%d", $lane_read_ref->{clusterCountRaw}{stdev});
             $lanehash->{clusterCountRaw}{mean}=sprintf("%d", $lane_read_ref->{clusterCountRaw}{mean});
             $lanehash->{clusterCountRaw}{sumsq}=sprintf("%d", $lane_read_ref->{clusterCountRaw}{sumsq});

             foreach my $stat (qw/clusterCountRaw percentClustersPF signalAverage2to4 signalLoss2to10 signalLoss10to20 signalAverage1to3/) {
                #stdev defined only if clusterCountPF>=2
                $lane_read_ref->{$stat}{stdev}=0  unless defined $lane_read_ref->{$stat}{stdev};
                $lane_read_ref->{$stat}{mean}=0   unless defined $lane_read_ref->{$stat}{mean};
                $lane_read_ref->{$stat}{sumsq}=0  unless defined $lane_read_ref->{$stat}{sumsq};

                $lanehash->{$stat}{stdev}=sprintf("%3.2f", $lane_read_ref->{$stat}{stdev});
                $lanehash->{$stat}{mean}=sprintf("%3.2f", $lane_read_ref->{$stat}{mean});
                $lanehash->{$stat}{sumsq}=sprintf("%3.2f", $lane_read_ref->{$stat}{sumsq});
             }

             foreach my $stat (qw/phasingApplied prephasingApplied/) {
                   $lanehash->{$stat} = sprintf("%6.4f", (defined($lane_read_ref->{$stat}) ? $lane_read_ref->{$stat} : 0) * 100);
             }

          }
          push @{$arrayhash->{'Lane'}}, { %$lanehash,
                                        'laneNumber' => $lane };

       }
       push @{$arrayhash->{'readNumber'}}, $read_num;
       push @{$results->{'Read'}}, $arrayhash;
    }

    $xml_results_ref->{'ExpandedLaneSummary'} = $results;
}

#------------------------------------------------------------------------------

sub print_xml_lane_results($$$$$$$) {
    my $laneResults_ref = shift;
    my $max_reads = shift;
    my $maxNumLanes = shift;
    my $lane_type_ref = shift;
    my $chip_stats_by_read_ref = shift;
    my $lane_num_reads_ref = shift;
    my $xml_results_ref = shift;


    my @irslt_stats = qw/clusterCountRaw clusterCountPF oneSig/;
    my @frslt_stats = qw/signal20AsPctOf1 percentClustersPF errorPF/;
    my @rslt_stats = (@irslt_stats, @frslt_stats);

    my $stat_formats = {};
    foreach my $stat (@irslt_stats) {
       $stat_formats->{$stat}='%d';
    }
    foreach my $stat (@frslt_stats) {
       $stat_formats->{$stat}='%3.2f';
    }

    my $results;
    for (my $read_num=1; $read_num <= $max_reads; $read_num++) {
       my $arrayhash;
       for (my $lane = 1; $lane <= $maxNumLanes; $lane++) {

          next if ($lane_num_reads_ref->[$lane] < $read_num);
  
          my $lanehash = {};
          if (exists $lane_type_ref->{$lane}) {
             my $lane_read_ref = $laneResults_ref->[$read_num][$lane];

             foreach my $stat (qw/clusterCountRaw clusterCountPF 
                               oneSig signal20AsPctOf1 percentClustersPF/) {
                if (defined($lane_read_ref->{$stat})) {
                   my $form = $stat_formats->{$stat};
                   $lanehash->{$stat}{mean} = sprintf($form, $lane_read_ref->{$stat}{mean});
                   $lanehash->{$stat}{sumsq} = sprintf($form, $lane_read_ref->{$stat}{sumsq});

                   #when tileCountPF<2, stdev is not defined
                   if (defined($lane_read_ref->{$stat}{stdev})) {
                      $lanehash->{$stat}{stdev} = sprintf($form, $lane_read_ref->{$stat}{stdev});
                   } else {
                      $lanehash->{$stat}{stdev} = 0;
                   }
                } 
             }

             my $yieldval = int($lane_read_ref->{yield}{total} / 1000);

             push @{$arrayhash->{'Lane'}}, { %$lanehash,
                                           'laneNumber' => $lane,
                                           'laneYield' => $yieldval };
          } else {
             push @{$arrayhash->{'Lane'}}, { %$lanehash,
                                           'laneNumber' => $lane };
          }
       }

       push @{$arrayhash->{'readNumber'}}, $read_num;
       push @{$results->{'Read'}}, $arrayhash;

   }

   $xml_results_ref->{'LaneResultsSummary'} = $results;
}

#------------------------------------------------------------------------------

sub print_xml_tile_results($$$$$$$) {
    my $laneResults_ref = shift;
    my $results_ref = shift;
    my $maxNumLanes = shift;
    my $lane_type_ref = shift;
    my $numTiles = shift;
    my $lane_num_reads_ref = shift;
    my $xml_results_ref = shift;

    my $results = {};
    my $aligned = 0;
    for (my $lane_num = 1; $lane_num <= $maxNumLanes; ++$lane_num) {
       my $lanehash = {};
       if (exists $lane_type_ref->{$lane_num}) {
          my $analysis_type = $lane_type_ref->{$lane_num};
          my $num_reads = $lane_num_reads_ref->[$lane_num];
          for (my $read_num = 1; $read_num <= $num_reads; ++$read_num) {
             my $readhash = {};

             my $lane_read_ref = $laneResults_ref->[$read_num][$lane_num];
             $aligned = $lane_read_ref->{alignedFlag};

             for (my $tileInd = 0; $tileInd < $numTiles; $tileInd++) {
                my $tilehash = {};
                my $tile_read_ref = $results_ref->[$read_num][$tileInd];

                next if (!defined($tile_read_ref->{lane}));
                next if ($tile_read_ref->{lane} != $lane_num);

                $tilehash->{'tileNumber'} = $tile_read_ref->{'tile'};

                my $percent_PF = $tile_read_ref->{percentClustersPF};

                $tilehash->{'clusterCountRaw'} = $tile_read_ref->{clusterCountRaw};
                $tilehash->{'clusterCountPF'} = $tile_read_ref->{clusterCountPF};
                $tilehash->{'oneSig'} = (sprintf "%3.2f", $tile_read_ref->{oneSig});

                if (defined($tile_read_ref->{signal20AsPctOf1})) {
                   $tilehash->{'signal20AsPctOf1'} = (sprintf "%3.2f", $tile_read_ref->{signal20AsPctOf1});
                } else {
                   $tilehash->{'signal20AsPctOf1'} = $not_avail_str;
                }

                $tilehash->{'percentClustersPF'} = (sprintf "%3.2f", $percent_PF);

                push @{$readhash->{'Tile'}}, $tilehash;
             }

             push @{$readhash->{'readNumber'}}, $read_num;
             push @{$lanehash->{'Read'}}, $readhash;
          }
                                       
      }

      $lanehash->{'laneNumber'} = $lane_num;
      push @{$results->{'Lane'}}, $lanehash;
   }

   $xml_results_ref->{'TileResultsByLane'} = $results;
}

#------------------------------------------------------------------------------

sub print_bustard_xml($$) {
   my $bustard_path = shift;
   my $xml_results_ref = shift;

   $xml_results_ref->{'Date'} = ( strftime $timeStampFormat, localtime );
   $xml_results_ref->{'Software'} = 'bcl2fastq-1.8.4';

   my $decl = "<?xml version=\"1.0\" ?>\n<?xml-stylesheet type=\"text\/xsl\"\n\thref=\"BustardSummary.xsl\" ?>\n";

   open (XML, ">$bustard_path/BustardSummary.xml") or die "Failed to open file $bustard_path/BustardSummary.xml";

   print XML XMLout($xml_results_ref, 
                 NoAttr=>1,
                 KeyAttr=>[],
                 XMLDecl=>$decl,
                 RootName=>'BustardSummary');

   close (XML);
}

#------------------------------------------------------------------------------

my $chipFile = "$runFolderPath/chip.xml";
my $chipName = $unknown;
read_chip_file($appName, $chipFile, \$chipName);

# Ensure sensible (or at least consistent) file ordering
@tiles = sort(@tiles);

my $numTiles = $tileInd;

# Calculate the Lane Results.
calc_lane_results(\@laneResults, \@results, $Gerald::Jerboa::maxNumLanes, \@lane_num_reads,
                  \%lane_type, $numTiles);

# Calculate weighted averages of statistics across the chip
# (by relevant tiles and by relevant clusters) and totals where relevant.
my @chip_stats_by_read = ();
my %chip_totals = ();
calc_chip_results(\@laneResults, $Gerald::Jerboa::maxNumLanes,
                  \@lane_num_reads, $num_available_reads,
                  \@chip_stats_by_read, \%chip_totals);

#now writes the hash and then the xml file
my %xml_results = ();

# Add the information from the sample sheet if available

my $defaultSampleSheet;
# default is csv file if it exists
foreach my $extension ('xml', 'csv')
{
    my $path = File::Spec->catfile($bustardPath, "SampleSheet.$extension");
    $defaultSampleSheet = $path if -e $path;
}
if ($defaultSampleSheet)
{
    my $sampleSheet = Casava::Demultiplex::SampleSheet::create($defaultSampleSheet);
    open my $handle, "<", $defaultSampleSheet
      or croak "ERROR: Could not open '$defaultSampleSheet' for reading: $!\n    ";
    $sampleSheet->load($handle);
    $xml_results{Samples}->{Lane} = [];
    foreach my $lane (sort $sampleSheet->laneNumbers)
    {
        next unless 1 == $sampleSheet->barcodes($lane);
        my ($barcode) = $sampleSheet->barcodes($lane);
        my $sampleId = $sampleSheet->sampleId($lane, $barcode);
        my $species = $sampleSheet->species($lane, $barcode);
        push @{$xml_results{Samples}->{Lane}}, {laneNumber => $lane, sampleId => $sampleId, barcode => $barcode, species => $species};
    }
}

print_xml_chip($machineName, $runFolder, $chipName, $tileArea,
               \%chip_totals, \%xml_results);

print_xml_lane_results(\@laneResults, $num_available_reads, $Gerald::Jerboa::maxNumLanes, \%lane_type,
                       \@chip_stats_by_read, \@lane_num_reads, \%xml_results);

print_xml_expanded_lane(\@laneResults, $num_available_reads, $Gerald::Jerboa::maxNumLanes, \%lane_type,
                        \@lane_num_reads, \%xml_results);

print_xml_tile_results(\@laneResults, \@results, $Gerald::Jerboa::maxNumLanes, \%lane_type, 
                        $numTiles, \@lane_num_reads, \%xml_results);

print_bustard_xml($bustardPath, \%xml_results);

