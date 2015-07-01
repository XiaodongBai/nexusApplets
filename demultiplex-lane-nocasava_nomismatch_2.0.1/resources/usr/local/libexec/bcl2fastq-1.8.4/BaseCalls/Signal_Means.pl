#!/usr/bin/env perl

# Copyright (c) 2006 Solexa
# Authors: 
# M. Parkinson, L. J. Davies
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).

use warnings;
use strict;
use Cwd;

# Plot average intensity versus cycle

my $usage=<<"END";
Usage: $0 [tileStatsDir] tileStatsFile[s]
tileStatsDir - directory containing tileStatsFiles (s_<lane>_<tile>_all.txt) 
END

die $usage unless (@ARGV >= 1);

my $numTiles;# later=scalar(@ARGV);



my ($laneNum,$cycleNum,$statsLine,$tileInd,$baseInd);
my $offset;

my $tileStatsDir=undef;
my $tileStatsFilePath;

my (@array,@stats, @tileAllInts, @tileCallInts, @tileCallNums);
my ($tileCurrBaseCallNum, $tileNumNotCalled, $tileNumClusters);
my ($numClusters, $currBaseCount, $currBasePercent);
my (@allInts,@callInts,@baseCounts);
my (%lanes, %cycleNums);
my @bases = qw(A C G T);
my $numBases = scalar(@bases);

my %baseIndMap =
    ( 'A'=>0,
      'C'=>1,
      'G'=>2,
      'T'=>3 );
my $uncalledBase = 'X';

my $thisCallKey;


for (my $argInd = 0; $argInd < @ARGV; ++$argInd){
    print STDERR "$ARGV[$argInd]\n";
}

if (-d $ARGV[0])
{
    $tileStatsDir=shift (@ARGV);
    $tileStatsDir.="/" unless ($tileStatsDir =~ /\/$/);
} # if

$numTiles = scalar(@ARGV);

if ($numTiles < 1)
{
    warn "No tile statistics files (s_*_all.txt) specified";
    exit (0);
}


for ($tileInd=0; $tileInd < $numTiles; $tileInd++)
{
    $tileStatsFilePath = $ARGV[$tileInd];

    if (defined($tileStatsDir))
    {
	$tileStatsFilePath = $tileStatsDir.$tileStatsFilePath;
    }

    print STDERR "Opening $tileStatsFilePath\n";
    open (TILE_STATS, $tileStatsFilePath)
	|| die "Error: failed to open tile stats file $tileStatsFilePath: $!";

    while (defined($statsLine=<TILE_STATS>)) {
	next if($statsLine =~ m/^\#/);

	@stats = split('\t',$statsLine);

	$laneNum = shift @stats;
	$cycleNum = shift @stats;
	@tileAllInts = splice(@stats, 0, $numBases);
	@tileCallInts = splice(@stats, 0, $numBases);
	@tileCallNums = splice(@stats, 0, $numBases);
	$tileNumNotCalled = shift @stats;

	$tileNumClusters = 0;
	$lanes{$laneNum}++;
	$cycleNums{$cycleNum}++;

	foreach (@bases) {
	    $baseInd = $baseIndMap{$_};
	    $tileCurrBaseCallNum = $tileCallNums[$baseInd];

	    $baseCounts[$laneNum][$cycleNum]{$_} += $tileCurrBaseCallNum;
	    $tileNumClusters += $tileCurrBaseCallNum;

	    $callInts[$laneNum][$cycleNum]{$_} 
	    += ($tileCallInts[$baseInd] * $tileCurrBaseCallNum);
	}

	$baseCounts[$laneNum][$cycleNum]{$uncalledBase} += $tileNumNotCalled;
	$tileNumClusters += $tileNumNotCalled;

	foreach (@bases) {
	    $baseInd = $baseIndMap{$_};

	    $allInts[$laneNum][$cycleNum]{$_} 
	    += ($tileAllInts[$baseInd] * $tileNumClusters);
	}
    }

    close (TILE_STATS)
        || die "Error: failed to close tile stats file $tileStatsFilePath: $!";
}

### Now need to output values 
## print a header line

print "#Lane\tCycle\tAll A\tAll C\tAll G\tAll T\tCall A\tCall C\tCall G\tCall T\tBase A\tBase C\tBase G\tBase T\n";

foreach $laneNum (sort {$a <=> $b} keys %lanes){
    ## only going to print out relevant lines

    foreach $cycleNum (sort {$a <=> $b} keys %cycleNums) {
	print "$laneNum\t$cycleNum\t";


	## Calculate divisors

	$numClusters = 0;

	foreach (@bases){
	    $baseCounts[$laneNum][$cycleNum]{$_} = 0
		unless defined ($baseCounts[$laneNum][$cycleNum]{$_});

	    $numClusters += $baseCounts[$laneNum][$cycleNum]{$_};
	}

	$baseCounts[$laneNum][$cycleNum]{$uncalledBase} = 0
	    unless defined ($baseCounts[$laneNum][$cycleNum]{$uncalledBase});
	
	$numClusters += $baseCounts[$laneNum][$cycleNum]{$uncalledBase};


	## all intensities - whether called base or not

	foreach (@bases){
	    $allInts[$laneNum][$cycleNum]{$_} = 0 
		unless defined ($allInts[$laneNum][$cycleNum]{$_});

	    if ($numClusters != 0){
		$allInts[$laneNum][$cycleNum]{$_} /= $numClusters;
	    }
	    printf ("%0.1f\t",$allInts[$laneNum][$cycleNum]{$_});
	}


	## intensities filtered by called base

	foreach (@bases) {
	    $callInts[$laneNum][$cycleNum]{$_} = 0 
		unless defined ($callInts[$laneNum][$cycleNum]{$_});

	    $currBaseCount = $baseCounts[$laneNum][$cycleNum]{$_};

	    if ($currBaseCount != 0) {
		$callInts[$laneNum][$cycleNum]{$_} /= $currBaseCount;
	    }
	    printf ("%0.1f\t",$callInts[$laneNum][$cycleNum]{$_});
	}


	## called bases by base as percentages of number of clusters
	## i.e. sum will be less than 100 if bases could not be called
	## for all clusters

	foreach (@bases){
	    $currBasePercent = 0;

	    if ($numClusters != 0){
		$currBasePercent 
		    = ($baseCounts[$laneNum][$cycleNum]{$_} / $numClusters)
		    * 100;
	    }

	    printf ("%0.1f\t", $currBasePercent);
	}

	print "\n";
    } # cycleNum loop
} # lane loop

