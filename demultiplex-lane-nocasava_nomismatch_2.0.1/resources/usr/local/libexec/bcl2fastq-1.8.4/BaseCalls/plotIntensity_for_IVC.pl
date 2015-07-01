#!/usr/bin/env perl

# File: plotIntensity_for_IVC.pl
# Copyright (c) 2006 Solexa
# Author: A. J. Cox, L. J. Davies
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).

# Description:
# Produces intensity cycle by cycle plot for lanes 

use warnings;
use strict;
use File::Spec;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Common::Utils;


my $maxNumLanes=8;
my $plotdir = "Plots";
my $tmpdir  = "Temp";

sub getTempName() {
    return File::Spec->catfile($tmpdir, "tmp-plotIntensity-$>-$$");
}

my $linewidth=3;

sub draw_graphs {
    my ($lane,$inputName,$outName,$title,$max) = @_;
    
    my $tmpName = getTempName().".txt";

    open (OUT, ">$tmpName") 
	|| die "Problem opening $tmpName for writing: $!";

    print OUT << "END";
set terminal postscript portrait color solid;
set yrange [0:$max];
set ytics 500;
set size 1.0, 1.0;
set origin 0.0, 0.0;
set title "$title";
set xlabel "Cycle Number";
set ylabel "Intensity";
plot '$inputName' using 1:2 title 'A' with linespoints lt 9 lw $linewidth pt 7 ps 2, '$inputName' using 1:3 title 'C' with linespoints lt 1 lw $linewidth pt 7 ps 2,    '$inputName' using 1:4 title 'G' with linespoints lt 3 lw $linewidth pt 7 ps 2, '$inputName' using 1:5 title 'T' with linespoints lt 2 lw $linewidth pt 7 ps 2;

END
    close (OUT);
    gnuplotImage($tmpName, $outName);
    unlink $tmpName;
}

sub draw_graphs_percent {
    my ($lane,$inputName,$outName,$title) = @_;

    my $tmpName = getTempName().".txt";

    open (OUT, ">$tmpName") 
	|| die "Problem opening $tmpName for writing: $!";

    print OUT << "END";
set terminal postscript portrait color solid;
set yrange [0:100];
set ytics 10;
set size 1.0, 1.0;
set origin 0.0, 0.0;
set title "$title";
set xlabel "Cycle Number";
set ylabel "Percent";
  plot '$inputName' using 1:2 title 'A' with linespoints lt 9 lw $linewidth pt 7 ps 2, '$inputName' using 1:3 title 'C' with linespoints lt 1 lw $linewidth pt 7 ps 2,    '$inputName' using 1:4 title 'G' with linespoints lt 3 lw $linewidth pt 7 ps 2, '$inputName' using 1:5 title 'T' with linespoints lt 2 lw $linewidth pt 7 ps 2;

END
    close (OUT);
    gnuplotImage($tmpName, $outName);
    unlink $tmpName;
}



### MAIN
die "Usage: $0 Signal_Means.txt Directory_Name\n"
    unless (@ARGV==2);

my $file      = shift(@ARGV);
my $directory = shift(@ARGV);
$directory = File::Spec->catdir($directory, $plotdir);

my @intensities;

open (IN, $file) 
    || die "Problem opening $file: $!";
while (<IN>)
{
  
    next if($_ =~ m/\#/);
    push (@intensities,[split('\t',$_)]);

}
close(IN);


my %lanes;
for(my $i=0;$i<scalar(@intensities);$i++){
    $lanes{$intensities[$i][0]}++;
}


### Going to plot graphs now


### Going to do CALL 
### Going to do call stuff first as will be highest intensities

my $max    = 0;
my $setmax = 0;

my $lane;
for($lane=1; $lane<=$maxNumLanes; $lane++){
    
    if (exists $lanes{$lane}){

	my $fileName = "s_".$lane."_call.png";
	my $outName  = File::Spec->catfile($directory, $fileName);
	print STDERR "$outName\n";

	my $tmpName=getTempName().".dat";
	open (OUT, ">$tmpName")
	    || die "Problem opening $tmpName for writing: $!";

	for(my $j=0;$j<scalar(@intensities);$j++){
	    next if ($intensities[$j][0] != $lane);

	    print  OUT "$intensities[$j][1]\t";
	    for(my $k=0;$k<4;$k++){
	 
		print  OUT "$intensities[$j][$k+6]\t";
		if($setmax == 0){
		
		    ### going to defined a max value based on first tile
		    if($intensities[$j][$k+6] > $max){
			$max = $intensities[$j][$k+6];
		    }	
		} 
	    }
	    print OUT "\n";
	}
	close (OUT);

	if($setmax == 0){
	    $max = $max * 1.5;
	    ## wanted to round to nearest 1000 so will cheat until i have more time
	    $max = $max / 1000;
	    $max = int($max);
	    $max = ($max+1) * 1000;
	    $setmax++;
	}
	draw_graphs($lane,$tmpName,$outName,
		    "Called intensities lane ".$lane, $max);
	unlink $tmpName;
    }
}


### Going to do ALL 
for($lane=1; $lane<=$maxNumLanes; $lane++){
    
    if (exists $lanes{$lane}){
	my $fileName = "s_".$lane."_all.png";
	my $outName  = File::Spec->catfile($directory, $fileName);
	print STDERR "$outName\n";

	my $tmpName=getTempName().".dat";
	open (OUT, ">$tmpName")
	    || die "Problem opening $tmpName for writing: $!";

	for(my $j=0;$j<scalar(@intensities);$j++){
	    next if ($intensities[$j][0] != $lane);

	    for(my $k=0;$k<5;$k++){
		print  OUT "$intensities[$j][$k+1]\t";
	    }
	    print OUT "\n";
	}
	close (OUT);
	draw_graphs($lane,$tmpName,$outName,
		    "All intensities lane ".$lane, $max);
	unlink $tmpName;
    }
}


### Going to do %ALL 
for($lane=1; $lane<=$maxNumLanes; $lane++){
   
    if (exists $lanes{$lane}){

	my $fileName = "s_".$lane."_percent_all.png";
	my $outName  = File::Spec->catfile($directory, $fileName);
	print STDERR "$outName\n";

	my $tmpName=getTempName().".dat";
	open (OUT, ">$tmpName")
	    || die "Problem opening $tmpName for writing: $!";

	for(my $j=0;$j<scalar(@intensities);$j++){
	    next if ($intensities[$j][0] != $lane);
	    my $total = 0;
	    for(my $k=0;$k<4;$k++){
		$total+= $intensities[$j][$k+2];
	    }
	    print  OUT "$intensities[$j][1]\t";
	    for(my $k=0;$k<4;$k++){
		if($total != 0){
		    $intensities[$j][$k+2] =  $intensities[$j][$k+2] / $total * 100;
		}
		else {
		    $intensities[$j][$k+2] = 0;
		}   
		print  OUT "$intensities[$j][$k+2]\t";
	    }
	    print OUT "\n";
	}
	close (OUT);

	draw_graphs_percent($lane,$tmpName,$outName,
			    "Percent All intensities lane ".$lane);
	unlink $tmpName;
    }
}


### Going to do %Call 
for($lane=1; $lane<=$maxNumLanes; $lane++){
    
    if (exists $lanes{$lane}){

	my $fileName = "s_".$lane."_percent_call.png";
	my $outName  = File::Spec->catfile($directory, $fileName);
	print STDERR "$outName\n";

	my $tmpName=getTempName().".dat";
	open (OUT, ">$tmpName")
	    || die "Problem opening $tmpName for writing: $!";

	for(my $j=0;$j<scalar(@intensities);$j++){
	    next if ($intensities[$j][0] != $lane);
	    my $total = 0;
	    for(my $k=0;$k<4;$k++){
		$total+= $intensities[$j][$k+6];
	    }
	    print  OUT "$intensities[$j][1]\t";
	    for(my $k=0;$k<4;$k++){
		if($total != 0){
		    $intensities[$j][$k+6] =  $intensities[$j][$k+6] / $total * 100;
		}
		else {
		    $intensities[$j][$k+6] =  0;
		}	
		print  OUT "$intensities[$j][$k+6]\t";
	    }
	    print OUT "\n";
	}
	close (OUT);
	draw_graphs_percent($lane,$tmpName,$outName,
			    "Percent Called intensities lane ".$lane);
	unlink $tmpName;
    }
}


### Going to do %Base Calls
for($lane=1; $lane<=$maxNumLanes; $lane++){
    
    if (exists $lanes{$lane}){

	my $fileName = "s_".$lane."_percent_base.png";
	my $outName  = File::Spec->catfile($directory, $fileName);
	print STDERR "$outName\n";

	my $tmpName=getTempName().".dat";
	open (OUT, ">$tmpName")
	    || die "Problem opening $tmpName for writing: $!";
    
	for(my $j=0;$j<scalar(@intensities);$j++){
	    next if ($intensities[$j][0] != $lane);

	    print  OUT "$intensities[$j][1]\t";
	    for(my $k=0;$k<4;$k++){
		print  OUT "$intensities[$j][$k+10]\t";
	    }
	    print OUT "\n";
	}
	close (OUT);

	draw_graphs_percent($lane,$tmpName,$outName,
			    "Percent Basecalls lane ".$lane);
	unlink $tmpName;
    }
}
