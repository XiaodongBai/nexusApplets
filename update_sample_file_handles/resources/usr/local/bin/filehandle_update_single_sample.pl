#!/usr/bin/env perl

# Written by Xiaodong Bai
# To combine all stats of multiple samples and convert to JSON format for uploading to LIMS

use strict;
use warnings;
use File::Basename;
use File::Find;
use lib '/usr/local/lib';
use JSON;
use vars qw(%json);

my ($PROJECT,$runId,$laneid,$sample,$tag) = @ARGV;
my $folder = "/data/runs/$runId";
$runId =~ /.+\_(\w{10})$/;
my $flowcell = $1;

# Hash with the suffixes for all sample-specific stats files
my %sampleSpecificStatsSuffix = (
    "realign.bam" => ".$tag.bam",
    "realign.bai" => ".$tag.bai",
    "gvcf.gz" => ".gvcf.gz",
    "gvcf.gz.tbi" => ".gvcf.gz.tbi",
    "vcf.gz" => ".vcf.gz",
    "vcf.gz.tbi" => ".vcf.gz.tbi",
    "stats-vcf.txt" => ".stats-vcf.txt",
    "stats-vcf.pdf" => ".stats-vcf.pdf",
    "stats-alignment-summary.txt" => ".$tag.stats-alignment-summary.txt",
    "stats-capture-metrics.txt" => ".$tag.stats-capture-metrics.txt",
    "stats-fastqc.html" => ".$tag.stats-fastqc.html",
    "stats-fastqc.txt" => ".$tag.stats-fastqc.txt",
    "stats-insertsize.pdf" => ".$tag.stats-insertsize.pdf",
    "stats-insertsize.txt" => ".$tag.stats-insertsize.txt",
    "stats-pertarget-coverage.txt" => ".$tag.stats-pertarget-coverage.txt",
    "gatk_readDepth_1x_q30.out.sample_cumulative_coverage_counts" => ".gatk_readDepth_1x_q30.out.sample_cumulative_coverage_counts",
    "gatk_readDepth_1x_q30.out.sample_cumulative_coverage_proportions" => ".gatk_readDepth_1x_q30.out.sample_cumulative_coverage_proportions",
    "gatk_readDepth_1x_q30.out.sample_interval_statistics" => ".gatk_readDepth_1x_q30.out.sample_interval_statistics",
    "gatk_readDepth_1x_q30.out.sample_interval_summary" => ".gatk_readDepth_1x_q30.out.sample_interval_summary",
    "gatk_readDepth_1x_q30.out.sample_statistics" => ".gatk_readDepth_1x_q30.out.sample_statistics",
    "gatk_readDepth_1x_q30.out.sample_summary" => ".gatk_readDepth_1x_q30.out.sample_summary",
    "gatk_readDepth_1x_q30.out" => ".gatk_readDepth_1x_q30.out",
);

my $outjson = "outputJSON";

# Run-specific information
$json{"JSONtype"} = "Sample";
$json{"FlowcellId"} = $flowcell;
$json{"RunId"} = $runId;
$laneid =~ /lane\.(\d)/;
$json{"LaneId"} = $1;

my $searchfolder = $folder."/$laneid";

my %sampleSpecific = ();
$sampleSpecific{"SampleId"} = $sample;
my @sampleSpecificFiles = ();
foreach my $key (sort keys %sampleSpecificStatsSuffix) {
    my %sampleFileSpecific = ();
    $sampleFileSpecific{"FileType"} = $key;
    my $localSampleFileName = $sample.$sampleSpecificStatsSuffix{$key};
    print STDERR $localSampleFileName,"\t";
    chomp(my $fileid = `dx find data --project=${PROJECT} --folder=${searchfolder} --name=${localSampleFileName} --brief | cut -d: -f2`);
    print STDERR $fileid,"\n";
    $sampleFileSpecific{"Location"} = {
	"\$dnanexus_link" => "$fileid",
	"local_name" => "$localSampleFileName",
    };
    push @sampleSpecificFiles, \%sampleFileSpecific;
}
$sampleSpecific{"FilesToUpload"} = \@sampleSpecificFiles;
$json{"Sample"} = \%sampleSpecific;

my $jsontext = to_json(\%json, { ascii => 1, pretty => 1, canonical => 1} );
open(J,">$outjson");
print J $jsontext;
close J;
