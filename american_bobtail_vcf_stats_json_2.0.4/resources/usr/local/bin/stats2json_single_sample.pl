#!/usr/bin/env perl

# Written by Xiaodong Bai
# To combine all stats of multiple samples and convert to JSON format for uploading to LIMS

use strict;
use warnings;
use File::Basename;
use File::Find;
use lib '/usr/local/lib';
use JSON;
use vars qw(@captureFiles @vcfStatFiles @insertStatFiles @alnStatFiles %capprejson %insprejson %vcfprejson %json %nexusids);
use vars qw($totalUnassignedReads $laneTotalReadNum %sampleReadNum %sampleReadPerc);
use vars qw(@runJsonArray @flowcellJsonArray @laneJsonArray @sampleJsonArray);

my ($PROJECT,$captureFile,$insertFile,$vcfStatsFile,$runId,$flowcell,$laneid,$sample,$startDate,$endDate,$idfile,$novelvarnum,$allvcffile,$obsgender) = @ARGV;
my $folder = "/data/runs/$runId";
$runId =~ /.+\_(\w{10})$/;
$flowcell = $1;

$captureFile =~ /\.(\w+)\.stats-capture-metrics/;
my $tag = $1;

# Get the DNAnexus IDs for the files being generated in the current applet
open(I,"$idfile");
while(<I>) {
    chomp;
    my @t = split /\s/;
    $nexusids{$t[0]} = $t[1];
}
close I;

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

# Process information in capture metrices files
my @captureStats = ();
my @captureHeaders = ();
open(F,"$captureFile");
while(<F>) {
    next if /^\#/;
    next unless /\t/;
    chomp;
    if (/^BAIT_SET/) {
	@captureHeaders = split /\t/;
    } else {
	@captureStats = split /\t/; # Extract the stats in an array
    }
}
close F;

my $captureHeadersRef = &wordCaseChange(\@captureHeaders);

my $lastCol = $#captureHeaders-3;
for my $col (0..$lastCol) {
    my $header = $$captureHeadersRef[$col];
    $captureStats[$col] = "na", unless ($captureStats[$col]); # Assign "na" to the missing stats
    $capprejson{$sample}{$header} = $captureStats[$col]; # Populate the stats hash
}

# Process information in insert size stat files
my @insertStats = ();
my @insertHeaders = ();
open(I,"$insertFile");
while(<I>) {
    next if /^\#/;
    next unless /\t/;
    last if /^insert_size/;
    chomp;
    if (/^MEDIAN_INSERT_SIZE/) {
	@insertHeaders = split /\t/;
    } elsif (/^\d/) {
	@insertStats = split /\t/;
    }
}
close I;

my $insertHeaderRef = &wordCaseChange(\@insertHeaders);

my $lastIsCol = $#insertHeaders-3;
for my $iscol (0..$lastIsCol) {
    my $header = $$insertHeaderRef[$iscol];
    $insertStats[$iscol] = "na", unless ($insertStats[$iscol]);
    $insprejson{$sample}{$header} = $insertStats[$iscol];
}

# Process information in VCF stats file
open(V,"$vcfStatsFile");
while(<V>) {
    chomp;
    if (/SNPs\:\t(\d+)/) {
        $vcfprejson{$sample}{"NumberSNPs"} = $1;
    }
    if (/indels\:\t(\d+)/) {
        $vcfprejson{$sample}{"NumberIndels"} = $1;
    }
    if (/multiallelic sites\:\t(\d+)/) {
        $vcfprejson{$sample}{"MultiallelicSites"} = $1;
    }
    if (/multiallelic SNP sites\:\t(\d+)/) {
        $vcfprejson{$sample}{"MultiallelicSNPSites"} = $1;
    }
    if (/^TSTV/) {
        my @t = split /\t/;
        my @tstvheaders = ("","","Ts","Tv","TsTvRatio","Ts1stALT","Tv1stALT","TsTvRatio1stALT");
        for my $col (2..7) {
            $vcfprejson{$sample}{$tstvheaders[$col]} = $t[$col];
        }
    }
    if (/^AF\t0\t0.000000\t(\d+)\t(\d+)\t(\d+)\t/) {
        $vcfprejson{$sample}{"HetSNPs"} = $1;
        $vcfprejson{$sample}{"HetTs"} = $2;
        $vcfprejson{$sample}{"HetTv"} = $3;
    }
    if (/^AF\t0\t99.000000\t(\d+)\t(\d+)\t(\d+)\t/) {
        $vcfprejson{$sample}{"HomSNPs"} = $1;
        $vcfprejson{$sample}{"HomTs"} = $2;
        $vcfprejson{$sample}{"HomTv"} = $3;
	$vcfprejson{$sample}{"HetHomRatio"} = sprintf("%.2f",$vcfprejson{$sample}{"HetSNPs"}/$1)
    }
    if (/^ST/) {
        my @t = split /\t/;
	$t[2] =~ s/>/to/;
        $vcfprejson{$sample}{$t[2]} = $t[3];
    }
    $vcfprejson{$sample}{"NumberNovelVariants"} = $novelvarnum;
    
}
close V;

for(my $i=1; $i<=60; $i++) {
    my $possize = "ISD".$i;
    my $negsize = "ISD_".$i;
    $vcfprejson{$sample}{$possize} = "0";
    $vcfprejson{$sample}{$negsize} = "0";
}

open(AV,"$allvcffile");
while(<AV>) {
    chomp;
    next unless (/^IDD/);
    my @t = split /\t/;
    my $size = $t[2];
    if ($size =~ /-/) {
	$size =~ s/-/_/;
    }
    $size = "ISD".$size;
    $vcfprejson{$sample}{$size} = $t[3];
}
close AV;

# subroutine to convert uppercase words in the original stats file to captilized first letter only
sub wordCaseChange {
    my $arrayref = shift;
    my @newarray = ();
    foreach my $elem (@$arrayref) {
	$elem = lc($elem);
	my @words = split("_",$elem);
	for my $i (0..$#words) {
	    $words[$i] =~ s/^(\w)/\u$1/;
	}
	push @newarray, join("",@words);
    }
    return \@newarray;
}

# Run-specific information
$json{"WorkflowName"} = "Illumina Sequence Analysis";
$json{"JSONversion"} = "2.0";
$json{"JSONtype"} = "Sample";
$json{"RunDetails"} = {
    "Version" => "American Bobtail",
    "Release" => "2.1.0",
    "Status" => "Success",
    "Description" => "",
};
$json{"FlowcellId"} = $flowcell;
$json{"RunId"} = $runId;
$laneid =~ /lane\.(\d)/;
$json{"LaneId"} = $1;

my $searchfolder = $folder."/$laneid";
print STDERR $PROJECT,"\n";
print STDERR $searchfolder,"\n";

my %sampleSpecific = ();
$sampleSpecific{"SampleStartDate"} = $startDate;
$sampleSpecific{"SampleEndDate"} = $endDate;
$sampleSpecific{"SampleId"} = $sample;
$sampleSpecific{"SequencedGender"} = $obsgender;
$sampleSpecific{"CaptureMetricsStats"} = $capprejson{$sample};
$sampleSpecific{"InsertSizeStats"} = $insprejson{$sample};
$sampleSpecific{"VCFStats"} = $vcfprejson{$sample};
my @sampleSpecificFiles = ();
foreach my $key (sort keys %sampleSpecificStatsSuffix) {
    my %sampleFileSpecific = ();
    $sampleFileSpecific{"FileType"} = $key;
    my $localSampleFileName = $sample.$sampleSpecificStatsSuffix{$key};
    print STDERR $localSampleFileName,"\t";
    my $fileid = "";
    if ($nexusids{$localSampleFileName}) {
	$fileid = $nexusids{$localSampleFileName};
    } else {
	chomp($fileid = `dx find data --project=${PROJECT} --folder=${searchfolder} --name=${localSampleFileName} --brief | cut -d: -f2`);
    }
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
