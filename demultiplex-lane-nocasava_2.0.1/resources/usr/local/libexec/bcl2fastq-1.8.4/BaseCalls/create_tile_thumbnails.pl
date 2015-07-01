#!/usr/bin/env perl

## This program takes the individual error plots and combines them
## into a single htm page which links to original files
## The page format has been designed by Tony Cox
# Copyright (c) Solexa 2006
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
use Cwd;
use Getopt::Long qw(:config no_ignore_case);
use POSIX qw(strftime);

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Common::Utils;

my $timeStampFormat = "%Y-%m-%d %H:%M:%S";

# To add support for new thumbnail types, add an appropriate entry 
# to %suffix (filename suffix of files to thumbnail)
# and %title (title for generated thumbnail page)

# hash of suffixes
my %suffix=
( "error"=>"_rescore.png",
  "base"=>"_base.png",
  "call"=>"_call.png",
  "all"=>"_all.png",
  "perfect"=>"_errors.png",
  "hist"=>"_hist.png",
  "info"=>"_info.png" );

# hash of page titles
my %title=
( "error"=>"Errors",
  "base"=>"Base",
  "call"=>"Calls",
  "all"=>"All",
  "perfect"=>"Error Curves",
  "hist"=>"Quality Score Histograms",
  "info"=>"Information Content" );
            
my $maxTiles = 0;
my @links;
my @optionKeys = ("maxTiles=i", "link=s@");
my %options;
GetOptions(\%options, @optionKeys);
if (defined($options{"maxTiles"})) {
    $maxTiles = $options{"maxTiles"};
}
if (defined($options{"link"})) {
    @links = @{$options{"link"}};
}

if ((@ARGV!=1)&&(@ARGV!=2))
{
    $"='/';
    my @array=sort keys %suffix;
    print STDERR 
        "Usage: $0 create_error_thumbnails.pl type [dir], where\n";
    print STDERR 
        "type = @array\n"; 
    print STDERR 
        "dir  = path of GERALD dir or omit to use current working dir\n";
    die;
} # if

my ($suffixToFind, $titleToUse);

if (defined($suffix{$ARGV[0]}))
{
    $suffixToFind=$suffix{$ARGV[0]};
    $titleToUse=$title{$ARGV[0]};
} # if
else
{
    $"='/';
    my @array=sort keys %suffix;
    print STDERR "$ARGV[0] not one of @array, assuming it's a suffix\n";
    $suffixToFind=$ARGV[0];
    unless ($suffixToFind=~/\.png$/)
    {
        $suffixToFind.='.png';
    } # unless
    unless ($suffixToFind=~/^_/)
    {
        $suffixToFind='_'.$suffixToFind;
    } # unless

    $titleToUse="Files of type *$ARGV[0]";
} # else


$ARGV[1]=getcwd() unless (@ARGV==2);


## Look for experimental number
my $expt_name=getRunFolderFromGerald($ARGV[1]);
$expt_name=getRunFolderFromBustard($ARGV[1]) unless defined $expt_name;
if (! defined($expt_name)) {
    warn "Could not parse experiment name from $ARGV[1]";
    $expt_name = "";
}

#= Look for list of tiles
my $numCols = 0;  # going to use this to work out whether manteia, 1 column or 2 column
my (%lane,%tile,%fileList);
%fileList = getTileList($ARGV[1]."/Plots/", $suffixToFind);
foreach my $key (keys %fileList) {
    $lane{$key} ++;
    foreach my $key2 (keys %{$fileList{$key}}) {
	$tile{$key2} ++;
    }
#   my $val;
#   my @array = split("_",$fileList{$key}{$key2});

# Next bit attempts to make the web page layout mimic the configuration
# of the tiles on the chip. Difficult to do this while a) chip layout is
# subject to change and b) there is no reliable way of determining what
# the layout is from the instrument output. Until these issues are resolved,
# always print tiles of each lane as one long column.  

#    if($array[2] < 8){
#	$val = 0;
#    }
#    elsif ($array[2] < 76){
#	$val = 1;
#    }
#    elsif ($array[2] < 151){
#	$val = 2;
#    }
#    else {
#	print "Unknown Tile Configuration tile = $val\n";
#	die;
#    }
#
#    if($val > $numCols){
#	$numCols = $val;
#    }
} # foreach

$numCols = 1;
my $step = 1;
if ($maxTiles > 0) {
    $step = 1 + ((scalar(keys %tile)-1) / ($maxTiles)); 
}

## Define output file name
print  "<html>\n<head>\n";

print  "<!--RUN_TIME ".( strftime $timeStampFormat, localtime )."-->\n";
print  '<!--SOFTWARE_VERSION bcl2fastq-1.8.4-->'."\n";

print  qq(</head>)."\n";

print  qq(<title>$expt_name $titleToUse</title>);
print  qq(<h1 align="center">$expt_name $titleToUse</h1>);

print "\n<body>\n";

## Starting the html output with the body name

foreach (@links) {
    print qq($_)."<br><br>\n";
}

print  qq(<table border="1" cellpadding="5"><tr><td><b>Tile</b></td>)."\n";

foreach (sort {$a<=>$b} keys %lane)
{
    my $l = $_;
    for (my $col = 0; $col < $numCols; $col++) {
	my $colId = "";
	if ($numCols > 1) {
	    $colId = chr(65+$col);
	}
        print  qq(<td><b>Lane ${_}${colId}</b></td>)."\n";
    } # for
} # foreach

print  qq(</tr>)."\n";

my $count = 0;

if($numCols < 2)
{
    my ($i, $j);
    foreach $i (sort {$a<=> $b} keys %tile)
    {
	$count ++;
	print  qq(<tr><td><b>$i</b></td>)."\n";
	foreach $j (sort {$a<=>$b} keys %lane)
	{
	    if(defined $fileList{$j}{$i})
	    {
		if ($count % $step == 0) {
		    print  qq(<td><a href="Plots/$fileList{$j}{$i}"> <img height=84 width=84 src="Plots/$fileList{$j}{$i}"></a></td>)."\n";
		}
		else {
		    print  qq(<td><a href="Plots/$fileList{$j}{$i}">Plots/$fileList{$j}{$i}</a></td>)."\n";
		}
	    } # if
	    else 
	    {
		print  qq(<td>&nbsp;</td>)."\n";
	    } # else
	} # foreach
	print  qq(</tr>)."\n";
    } # foreach

} # if
else 
{
    # Attempts to do two sublanes per lane. Needs rewriting
    my $numRows    = 75;
    my ($i, $j, $k);
    for ($i=1;$i<=$numRows;$i=$i+$numCols)
    {
	my @t = ();
	for ($k= 0; $k < $numCols; $k++) {
	    push (@t, sprintf("%04i", $i+$k));
	}
	my $flag = 0;
	for ($k= 0; $k < $numCols; $k++) {
	    if (defined $tile{$t[$k]}) {
		$flag = 1;
	    }
	}
	if($flag) 
	{
	    my $label = join(' / ', @t);
	    print  qq(<tr><td><b>$label </b></td>)."\n";

	    foreach $j (sort {$a<=>$b} keys %lane)
	    {
		for ($k = 0; $k < $numCols; $k++) {
		    if(defined $fileList{$j}{$t[$k]}) {
			print  qq(<td><a href="Plots/$fileList{$j}{$t[$k]}"> <img height=84 width=84 src="Plots/$fileList{$j}{$t[$k]}"></a></td>)."\n";
		    } # if
		    else 
		    {
			print  qq(<td>&nbsp;</td>)."\n";
		    } # else
		}
	    } # foreach $j
	    print  qq(</tr>)."\n";
	} # if
    } # for i
} # else


print qq(</table>)."\n";
print qq(</body>)."\n";
print qq(</html>)."\n";


