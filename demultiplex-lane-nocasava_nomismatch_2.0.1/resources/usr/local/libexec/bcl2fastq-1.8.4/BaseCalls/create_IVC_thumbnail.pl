#!/usr/bin/env perl


## This program takes the individual IVC plots and combines them
## into a single htm page which links to orginal files
## The page format has been designed by Tony Cox

# Copyright (c) 2006 Solexa
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).

### The program is called using directory name

use warnings;
use strict;
use POSIX qw(strftime);

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Common::Utils;

my $id = '';

my $timeStampFormat = "%Y-%m-%d %H:%M:%S";

my ($file,$i,$j,$output_name,$expt_name,$deblock);
my (%lane,%tile,%results);


if(@ARGV != 1){
    print STDERR "create_IVC_thumbnail.pl Expt_dir\n";
    die;
}

my $directory = $ARGV[0];

## Look for experimental number
$expt_name=getRunFolderPath($directory);
if (! defined($expt_name)) {
    warn "Could not parse experiment name from ${directory}";
    $expt_name = "";
}

print  "<html><head>\n";

print q(<!--RUN_TIME ) . ( strftime $timeStampFormat, localtime ) . q(-->) . "\n";
print q(<!--SOFTWARE_VERSION ) . $id . q(-->) . "\n";

print  "</head><body>\n";

print  qq(<title>$expt_name IVC</title>);

print  qq(<h1 align="center">$expt_name Intensity Plots<br><br> </h1>);


## Starting the html output with the body name


opendir (DIR,$ARGV[0]."/Plots/") or die "can't open $ARGV[0]\n";

$deblock = 0;

while (defined($file = readdir (DIR))){

    if($file =~ m/_deblock.png/){
	$deblock++;
    }

    next if ($file !~ m/s_[0-9]_all\.png$/);
    next if ($file =~ m/[%]/);
    #print "$file\n";
 
 #   open(IN,$file) || die "Can't open $file\n";
 #   close(IN);

    my @array = split("_",$file);
    
    $results{$array[1]} = $file;
 #   print "$array[1]\n";
    $lane{$array[1]}++;
 
}


print  qq(<table border="1" cellpadding="5"><tr><td><b>Lane</b></td>);
print  qq(<td><b>All</b></td>);
print  qq(<td><b>Called</b></td>);
if($deblock > 0){
    print  qq(<td><b>Deblock</b></td>);
}
print  qq(<td><b>% Base Calls</b></td>);
print  qq(<td><b>% All</b></td>);
print  qq(<td><b>% Called</b></td>);

print  qq(</tr>)."\n";

    foreach $j (sort {$a<=>$b} keys %lane){


	print  qq(<tr><td><b>$j</b></td>);
    
	if(defined $results{$j}){
	    $file = $results{$j};
	    print  qq(<td><a href="Plots/$file"> <img height=84 width=126 src="Plots/$file"></a></td>);
	}
	else {
	    print  qq(<td>&nbsp;</td>);
	}
	$file = $results{$j};
	$file =~ s/_all.png/_call.png/;
	if(defined $file){
	    print  qq(<td><a href="Plots/$file"> <img height=84 width=126 src="Plots/$file"></a></td>);
	}
	else {
	    print  qq(<td>&nbsp;</td>);
	}


	if($deblock > 0){

	    $file = $results{$j};
	    $file =~ s/_all.png/_deblock.png/;
	    if(defined $file){
		print  qq(<td><a href="Plots/$file"> <img height=84 width=126 src="Plots/$file"></a></td>);
	    }
	    else {
		print  qq(<td>&nbsp;</td>);
	    }

	}
	
	
	$file = $results{$j};
	$file =~ s/_all.png/_percent_base.png/;
	if(defined $file){
	    print  qq(<td><a href="Plots/$file"> <img height=84 width=126 src="Plots/$file"></a></td>);
	}
	else {
	    print  qq(<td>&nbsp;</td>);
	}



	$file = $results{$j};
	$file =~ s/_all.png/_percent_all.png/;
	if(defined $file){
	    print  qq(<td><a href="Plots/$file"> <img height=84 width=126 src="Plots/$file"></a></td>);
	}
	else {
	    print  qq(<td>&nbsp;</td>);
	}
    	$file = $results{$j};
	$file =~ s/_all.png/_percent_call.png/;

	if(defined $file){
	    print  qq(<td><a href="Plots/$file"> <img height=84 width=126 src="Plots/$file"></a></td>);
	}
	else {
	    print  qq(<td>&nbsp;</td>);
	}

    
    print  qq(</tr>)."\n";
}


print "</table></body></html>\n";
