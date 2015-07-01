#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2007-2011 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And
VAriation (CASAVA) software package.

=cut

=head1 SYNOPSIS

configureQseqToFastq.pl --input-dir DIR [options] | --help

=head2 SUMMARY

Create a makefile to convert a directory of qseq files to a directory
tree of compressed fastq files following CASAVA 1.8 filename and
directory structure conventions. If detected, configuration data
used by GERALD are also transfered to the output directory.

This script will not configure demultiplexing. The input directory
must contain qseq files which are either non-demultiplexed or already
demultiplexed by another utility.

=head2 OPTIONS

=over 4

=item --input-dir DIRECTORY

Path to qseq directory (no default)

=item --output-dir DIRECTORY

Path to root of CASAVA 1.8 unaligned directory structure. Directory
will be created if it does not exist (default:
'<input-dir>/QseqToFastq/Unaligned')

=item --fastq-cluster-count INTEGER

Maximum number of fastq records per fastq file (default: 4000000)

=item --config-file FILENAME

Specify the Bustard config file to be copied to the fastq
directory (default: '<input-dir>/config.xml')

=item --summary-file FILENAME

Specify the Bustard summary file to be copied to the fastq
directory (default: '<input-dir>/BustardSummary.xml')

=item --flowcell-id STRING

Use the specified string as the flowcell id. (default value is parsed
from the config-file)

=back

=cut

use warnings "all";
use strict;

use Carp;
use File::Copy qw(copy);
use File::Spec;
use File::Spec::Win32;
use Getopt::Long;
use Pod::Usage;
use XML::Simple;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Common::Utils qw(parseRunFolderName reallyRealPath);
use Casava::Demultiplex::SampleSheet;
use Casava::Demultiplex::DemuxUtil;

my $defaultOutDirSuffix=File::Spec->catdir("QseqToFastq","Unaligned");
my $projectName="default";
my $projectDirPrefix="Project_";

my $scriptName = (File::Spec->splitpath($0))[2];
my $argCount=scalar(@ARGV);
my $help;

my $inputDir;
my $outputDir;
my $fastqClusterCount=4000000;
my $configFile;
my $summaryFile;
my $userFlowCellID;


my $result = GetOptions('help|h' => \$help,
                        'input-dir=s' => \$inputDir,
                        'output-dir=s' => \$outputDir,
                        'fastq-cluster-count=i' => \$fastqClusterCount,
                        'config-file=s' => \$configFile,
                        'summary-file=s' => \$summaryFile,
                        'flowcell-id=s' => \$userFlowCellID) or pod2usage(2);

pod2usage(2) if($help or (not $argCount));

croak "ERROR: must specify --input-dir argument\n" unless(defined $inputDir);
croak "ERROR: input directory argument: '$inputDir' does not exist or is not a directory\n"  unless -d $inputDir;
croak "ERROR: --fastq-cluster-count must be greater than 0\n"  unless ($fastqClusterCount>0)  ;

if(reallyRealPath($inputDir)) { croak("ERROR: Can't interpret input directory argument: '$inputDir'.\n"); }

unless (defined $outputDir) {
    $outputDir=File::Spec->catdir($inputDir,$defaultOutDirSuffix);
}

if((defined $configFile) and (not -f $configFile)) {
    croak "ERROR: cannot find specified config-file: '$configFile'\n";
}

if((defined $summaryFile) and (not -f $summaryFile)) {
    croak "ERROR: cannot find specified summary-file: '$summaryFile'\n";
}

if(defined $userFlowCellID) {
    croak "ERROR: invalid flowcell id name: '$userFlowCellID'" if(($userFlowCellID eq '') or 
                                                                  ($userFlowCellID =~ /[\s\*\/]/));
}



sub sort_qseq_names($$) {
    my $qseq_pattern = '^s_[1-8]_[1-3]_([0-9]{4})_qseq.txt$';
    my ($v0) = ($_[0] =~ /$qseq_pattern/);
    my ($v1) = ($_[1] =~ /$qseq_pattern/);
    return $v0 <=> $v1;
}



# 0: enumerate qseqs and lane/read combinations:
#
my @qseqs;
my %lanes;
{
    opendir(my $dh, $inputDir);
    @qseqs = sort sort_qseq_names grep { /^s_[1-8]_[1-3]_[0-9]{4}_qseq.txt$/ && -f "$inputDir/$_" } readdir($dh);
    closedir($dh);

    unless(@qseqs) {
        croak("ERROR: no qseq files found in directory: $inputDir\n");
    }

    # get total list of lanes:
    #
    my %tileCheck;
    for my $qseqFile (@qseqs) {
        my ($lane,$read,$tile) = ($qseqFile =~ /^s_([1-8])_([1-3])_([0-9]{4})_qseq.txt$/);
        $lanes{int($lane)}{int($read)} = undef;
        $tileCheck{int($lane)}{int($tile)}{int($read)} = undef;
    }

    # sanity check for obvious missing tiles:
    #
    for my $lane (keys %lanes) {
        my @laneReads = keys %{$lanes{$lane}};
        for my $tile (keys %{$tileCheck{$lane}}) {
            for my $read (@laneReads) {
                next if(exists $tileCheck{$lane}{$tile}{$read});
                my $missingFile = "$inputDir/s_$lane\_$read\_$tile\_qseq.txt";
                croak("ERROR: The following qseq file appears to be missing: $missingFile\n");
            }
        }
    }
}


# 1: create output directory:
#
checkMakeDir($outputDir);
if(reallyRealPath($outputDir)) { croak("ERROR: Can't interpret output directory argument: '$outputDir'.\n"); }


# 2: transfer configuration data:
#
my $configxmlPath=(defined $configFile ? $configFile : File::Spec->catfile($inputDir,"config.xml"));
my $flowCellID="";

if(-f $configxmlPath) {
    my $newConfigPath=File::Spec->catfile($outputDir,"DemultiplexedBustardConfig.xml");
    copy($configxmlPath,$newConfigPath)
      or croak("ERROR: failed to copy file '$configxmlPath' to '$newConfigPath' $!\n");

    # attempt to recover flowCellID:
    unless(defined $userFlowCellID)
    {
        my $configRef = XMLin($configxmlPath, SuppressEmpty => 1, ForceArray => ['Lane', 'SelectedTiles', 'Tile', 'TileRange'])
          or croak "ERROR: couldn't load config file: '$configxmlPath' $!\n";
        croak "ERROR: 'Run' element missing from config file: '$configxmlPath'\n" unless exists $configRef->{Run};
        my $run = $configRef->{Run};
        if((defined  $run->{RunParameters}) and (defined $run->{RunParameters}->{RunFolder})){
            my $runFolder = $run->{RunParameters}->{RunFolder};
            my $runFolderDirName = (File::Spec->splitpath($runFolder))[2];
            my $parsedName = parseRunFolderName($runFolderDirName);
            if($runFolder and (not $parsedName)) {
                # assume flowcell-id parse failed because
                # runFolder is a windows path from RTA:
                $runFolderDirName = (File::Spec::Win32->splitpath($runFolder))[2];
                $parsedName = parseRunFolderName($runFolderDirName);
            }
            $flowCellID=$parsedName->[3] if($parsedName);
        }
    } else {
        $flowCellID = $userFlowCellID;
    }
    $projectName = $flowCellID if($flowCellID);
}


my $summaryxmlPath=(defined $summaryFile ? $summaryFile : File::Spec->catfile($inputDir,"BustardSummary.xml"));
if(-f $summaryxmlPath) {
    my $newSummaryPath=File::Spec->catfile($outputDir,"DemultiplexedBustardSummary.xml");
    copy($summaryxmlPath,$newSummaryPath)
      or croak("ERROR: failed to copy file '$summaryxmlPath' to '$newSummaryPath' $!\n");
}


# create default SampleSheet:
{
    my $sampleSheetData = "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n";
    for my $lane (sort keys %lanes) {
        my $sampleName = defaultSampleName($lane);
        $sampleSheetData .="$flowCellID,$lane,$sampleName,,,,N,,,$projectName\n";
    }
    my $sampleSheetPath = File::Spec->catfile($outputDir,"SampleSheet.csv");
    open(my $sampleSheetFH,">$sampleSheetPath");
    print $sampleSheetFH $sampleSheetData;
    close($sampleSheetFH);

}

# 3: create directory structure for new fastq files and temporaries:
#
my $projectDirName=$projectDirPrefix . $projectName;
my $projectPath=File::Spec->catdir($outputDir,$projectDirName);
my $projectPathMake=File::Spec->catdir("\$(ROOT_DIR)","\$(PROJECT_DIRNAME)");
checkMakeDir($projectPath);
my $tempProjectPath=File::Spec->catdir($outputDir,"Temp",$projectDirName);
my $tempProjectPathMake=File::Spec->catdir("\$(ROOT_DIR)","Temp","\$(PROJECT_DIRNAME)");
checkMakeDir($tempProjectPath);
{
    for (keys %lanes) {
        my $sampleDirName = defaultSampleDirName($_);
        my $samplePath=File::Spec->catdir($projectPath,$sampleDirName);
        checkMakeDir($samplePath);
    }
}



# 4. create makefile to convert qseqs to fastqs:
#
my $DATA_DIR = '/usr/local/share/bcl2fastq-1.8.4';

my $makeFile = File::Spec->catfile($outputDir,"Makefile");
open(my $MAKEFH, "> $makeFile") || croak("ERROR: failed to open file: $makeFile\n");


my $fccmd = ($flowCellID ? "--fc=$flowCellID " : "" );


print $MAKEFH <<ENDE;
# This makefile was automatically generated by $scriptName
# Please do not edit.

MAKEFILES_DIR:=$DATA_DIR/makefiles

# Import the global configuration
include \$(MAKEFILES_DIR)/Config.mk

# Import the debug functionalities
include \$(MAKEFILES_DIR)/Debug.mk

FQCONV := \$(BIN_DIR)/FastqConverter
FQFINISH := \$(LIBEXEC_DIR)/BaseCalls/finishFastq.pl

QSEQ_DIR := $inputDir
ROOT_DIR := $outputDir
PROJECT_DIRNAME := $projectDirName
PROJ_DIR := $projectPathMake
TEMP_PROJ_DIR := $tempProjectPathMake

FASTQ_SUFFIX := $fastqSuffix
COMPRESSED_FASTQ_SUFFIX := $compressedFastqSuffix

all: all_lanes

.PHONY : all all_lanes

ENDE

my @laneTargets;

for my $lane (sort keys %lanes) {
    my $laneTarget = "l$lane";
    push @laneTargets, $laneTarget;

    my $sampleDirName = defaultSampleDirName($lane);
    my $tempSamplePathMake = "\$(TEMP_PROJ_DIR)/$sampleDirName";
    my $samplePathMake="\$(PROJ_DIR)/$sampleDirName";

    my @laneReadTargets;

    for my $read (sort keys %{$lanes{$lane}}) {
        my $laneReadTarget = "l$lane\_r$read";
        push @laneReadTargets, $laneReadTarget;

        my $fastqPrefix=defaultFastqPrefix($lane,$read);
        my $fastqPathPrefix=File::Spec->catfile($tempSamplePathMake,$fastqPrefix);
        my $finalPathPrefix=File::Spec->catfile($samplePathMake,$fastqPrefix);

        my @lrqseqs = grep { /^s_$lane\_$read\_[0-9]{4}_qseq.txt$/ } @qseqs;
        my @lrfastqs = map { s/_qseq.txt/$fastqSuffix/; $_ } @lrqseqs;

        print $MAKEFH <<ENDE;

$tempSamplePathMake/s_$lane\_$read\_\%\$(FASTQ_SUFFIX): \$(QSEQ_DIR)/s_$lane\_$read\_\%_qseq.txt
\tmkdir -p $tempSamplePathMake; \\
\t\$(FQCONV) --read=$read --it qseq $fccmd--no-compression --in=\$^ --out=\$@

l$lane\_r$read\_fastq_filenames := @lrfastqs
l$lane\_r$read\_fastqs := \$(patsubst \%, $tempSamplePathMake/\%, \$(l$lane\_r$read\_fastq_filenames) )

.SECONDARY : \$(l$lane\_r$read\_fastqs)

${finalPathPrefix}001\$(COMPRESSED_FASTQ_SUFFIX): \$(l$lane\_r$read\_fastqs)
\t\$(FQFINISH) --sample-tmp-dir=$tempSamplePathMake --sample-output-dir=$samplePathMake --lane=$lane --read=$read --fastq-cluster-count=$fastqClusterCount \$(l$lane\_r$read\_fastq_filenames)

$laneReadTarget: ${finalPathPrefix}001\$(COMPRESSED_FASTQ_SUFFIX)

ENDE

    }
    print $MAKEFH "$laneTarget: @laneReadTargets \n\n";
    print $MAKEFH ".PHONY : @laneReadTargets\n\n";
}

print $MAKEFH "all_lanes: @laneTargets \n\n";
print $MAKEFH ".PHONY : @laneTargets\n\n";
close($MAKEFH);

print STDERR "\n$scriptName: Successfully created makefile: '$makeFile'\n\n";

1;

__END__

