#!/usr/bin/env perl

# Written by Xiaodong Bai
# To combine all stats of multiple samples and convert to JSON format for uploading to LIMS

use strict;
use warnings;
use File::Basename;
use File::Find;
use lib '/usr/local/lib';
use JSON;
use vars qw(@captureFiles @vcfStatFiles @insertStatFiles @alnStatFiles %prejson %json %nexusids);
use vars qw($totalUnassignedReads $laneTotalReadNum %sampleReadNum %sampleReadPerc);
use vars qw(@runJsonArray @flowcellJsonArray @laneJsonArray @sampleJsonArray);

my ($PROJECT,$demulFile,$runId,$flowcell,$laneid,$startDate,$endDate,$idfile) = @ARGV;
my $folder = "/data/runs/$runId";
$runId =~ /.+\_(\w{10})$/;
$flowcell = $1;
if ($laneid !~ /lane/) {
    $laneid = "lane.$laneid";
}

open(I,"$idfile");
while(<I>) {
    chomp;
    my @t = split(" ",$_);
    $nexusids{$t[0]} = $t[1];
}
close I;

# Hash with the suffixes for all lane-specific stats files
my %laneSpecificStatsSuffix = (
    "BustardSummary.xml" => ".BustardSummary.xml",
    "Demultiplex_Stats.htm" => ".Demultiplex_Stats.htm",
    "Flowcell_demux_summary.xml" => ".Flowcell_demux_summary.xml",
);

my $outjson = "outputJSON";

# Process the lane-specific information from the demultiplexing stats file
my @sampleReadNumbers = ();
my @demulHeaders = ();
my @sampleIds = ();
my $sampleId = "";
open(H,"$demulFile");
my $colCount = 0;
while(<H>) {
    chomp;
    if (/\<th\>(.+)\<\/th\>/) {
	push @demulHeaders, $1;
    }
    if (/\<td\>(.+)\<\/td\>/) {
	my $info = $1;
	if ($demulHeaders[$colCount] eq "Sample ID") {
	    $sampleId = $info;
	    push @sampleIds, $sampleId, unless ($sampleId =~ /lane\d/);
	    $colCount++;
	    next;
	}
	if ($demulHeaders[$colCount] eq "# Reads") {
	    $info =~ s/\,//g;
	    $info = int($info);
	    $laneTotalReadNum += $info;
	    if ($sampleId =~ /lane\d/) {
		$totalUnassignedReads = $info;
	    } else {
		push @sampleReadNumbers, $info;
	    }
	    $colCount++;
	    next;
	}
	$colCount++;
    }
    if (/\<\/tr\>/) {
	$colCount = 0;
	$sampleId = "";
    }
    last if (/\<h2\>Sample\sinformation\<\/h2\>/);
}
close H;

for my $i (0..$#sampleIds) {
    my $sampleid = $sampleIds[$i];
    my $samplereadcnt = $sampleReadNumbers[$i];
    my $samplereadperc = sprintf("%.6f",$samplereadcnt/$laneTotalReadNum);
    $prejson{$sampleid} = $samplereadperc;
}
my $percUnassignedReads = sprintf("%.6f",$totalUnassignedReads/$laneTotalReadNum);
#$prejson{$sample}{"TotalReads"} = $sampleReadNumber;
#my $percReadsPerSample = sprintf("%.6f",$sampleReadNumber/$laneTotalReadNum);
#$prejson{$sample}{"PercentReadsPerSample"} = $percReadsPerSample;

# Run-specific information
$json{"WorkflowName"} = "Illumina Sequence Analysis";
$json{"JSONversion"} = "1.0";
$json{"JSONtype"} = "Lane";
$json{"RunDetails"} = {
    "Version" => "American Bobtail",
    "Release" => "2.0.0",
    "Status" => "Success",
    "Description" => "",
};
$json{"FlowcellId"} = $flowcell;
$json{"RunId"} = $runId;
$laneid =~ /lane\.(\d)/;
$json{"LaneId"} = $1;

# prepare perl data structure for JSON output
my %laneSpecific = ();
$laneSpecific{"DemuxStartDate"} = $startDate;
$laneSpecific{"DemuxEndDate"} = $endDate;
$laneSpecific{"Stats"} = {
    "TotalReads" => "$laneTotalReadNum",
    "PercentSampleReadsPerLane" => \%prejson,
    "TotalUnassignedReads" => "$totalUnassignedReads",
    "PercentUnassignedTotalReads" => "$percUnassignedReads",
};
my @laneSpecificFiles = ();
#foreach my $laneKey (sort keys %laneSpecificStatsSuffix) {
foreach my $laneKey (sort keys %laneSpecificStatsSuffix) {
    my %laneFileSpecific = ();
    $laneFileSpecific{"FileType"} = $laneKey;
    my $localLaneFileName = "run.".$runId.".".$laneid.$laneSpecificStatsSuffix{$laneKey};
    print STDERR $localLaneFileName,"\t";
#    chomp(my $fileid = `dx find data --project=${PROJECT} --name=${localLaneFileName} --brief | cut -d: -f2`);
    my $fileid = $nexusids{$localLaneFileName};
    print STDERR $fileid,"\n";
    $laneFileSpecific{"Location"} = {
	"\$dnanexus_link" => "$fileid",
	"local_name" => "$localLaneFileName",
    };
    push @laneSpecificFiles, \%laneFileSpecific;
}
$laneSpecific{"FilesToUpload"} = \@laneSpecificFiles;
$json{"Lane"} = \%laneSpecific;
my $jsontext = to_json(\%json, { ascii => 1, pretty => 1, canonical => 1} );
open(J,">$outjson");
print J $jsontext;
close J;
