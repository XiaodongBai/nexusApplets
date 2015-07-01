#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2009-2011 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And VAriation
(CASAVA) software package.

=cut

=head1 NAME

configureBclToFastq.pl

Demultiplex a Base Call directory based upon a given set of barcodes.

=cut

=head1 SYNOPSIS

=over 12

=item B<configureBclToFastq.pl>

[S<B<--adapter-sequence> I<adapter_fasta_file_path>>]
[S<B<--use-bases-mask> I<mask>>]
[S<B<--no-eamss>>]
[S<B<--with-failed-reads>>]
[S<B<--input-dir> I<base_calls_dir>>]
[S<B<--intensities-dir> I<intensities_dir>>]
[S<B<--positions-dir> I<positions_dir>>]
[S<B<--positions-format> I<.locs .clocs or _pos.txt>>]
[S<B<--filter-dir> I<filter_dir>>]
[S<B<--output-dir> I<output_dir>>]
[S<B<--sample-sheet> I<sample_sheet>>]
[S<B<--mismatches> I<num_of_mismatches per barcode component>>]
[S<B<--fastq-cluster-count> I<cluster_count>>]
[S<B<--ignore-missing-stats>>]
[S<B<--ignore-missing-bcl>>]
[S<B<--ignore-missing-control>>]
[S<B<--flowcell-id> I<flow_cell_id>>]
[S<B<--tiles> I<tile_selection>>]

=item B<configureBclToFastq.pl>

B<--help> or B<--man>

=back

=head1 OPTIONS AND ARGUMENTS

=over 8

=item B<--adapter-sequence> I<adapter_fasta_file_path>

Path to a multi-contig FASTA file (multiple paths allowed), where each contig specifies an adapter. 
If only one FASTA file is provided, the same (set of) adapter(s) will be used across all reads.
Otherwise, Read 1 take its adapters from file #1, Read 2 takes its adapters from file #2, and so on.

=item B<--adapter-stringency> I<matching_fraction>

When applying adapter trimming/masking, this is the fraction of bases which must match the adapter 
in order to be applied. Set to 1.0 to allow no mismatches at all (more selective), or to 0.66...7 
(2/3) to replicate CASAVA's behaviour. 

(defaults to 0.9)

=item B<--use-bases-mask> I<mask>[[I<,mask>]...]

Conversion mask characters:

  - Y or y: use
  - N or n: discard
  - I or i: use for indexing

If not given, the mask will be guessed from the B<RunInfo.xml> file in the run folder.

For instance, in a 2x76 indexed paired end run, the mask I<Y76,I6n,y75n> means: 
"use all 76 bases from the first end, discard the last base of the indexing read, and 
use only the first 75 bases of the second end". 

=item B<--no-eamss> Disable the masking of the quality values with EAMSS.

Z<>

=item B<--with-failed-reads> Include failed reads into the FASTQ files.

(by default, only reads passing filter are included)

Z<>

=item B<--input-dir> I<base_calls_dir>

Path to a valid BaseCalls directory (defaults to current dir)

=item B<--intensities-dir> I<intensities_dir>

Path to a valid Intensities directory (defaults to parent of I<base_calls_dir>)

=item B<--positions-dir> I<positions_dir>

Path to a directory containing positions files.
(defaults depending on RTA version)

=item B<--positions-format> I<.locs .clocs or _pos.txt>

Format of the input cluster positions information. (defaults to I<.clocs>)

=item B<--filter-dir> I<filter_dir>

Path to a directory containing filter files. (defaults depending on RTA version)

=item B<--output-dir> I<output_dir>

Path to the demultiplexed output (defaults to I<base_calls_dir>/../../../Unaligned)

=item B<--sample-sheet> I<sample_sheet>

Path to SampleSheet.csv (defaults to I<base_call_dir>/SampleSheet.csv)

=item B<--mismatches> I<num_of_mismatches per barcode component>

Comma-separated list of mismatches allowed for each barcode component. 
Either 0 or 1 (defaults to 0)

=item B<--fastq-cluster-count> I<cluster_count>

Maximum number of fastq records per fastq file.

(default: 4000000)

=item B<--ignore-missing-stats>

Fill in with zeros when *.stats files are missing

=item B<--ignore-missing-bcl>

Interpret missing *.bcl files as no call

=item B<--ignore-missing-control>

Interpret missing control files as not-set control bits

=item B<--tiles> I<regex>[[I<,regex>]...]

Comma-separated list of regular expressions to select only a subset of the tiles available in the flow-cell.

 - to select all the tiles ending with "5" in all lanes: --tiles [0-9][0-9][0-9]5
 - to select tile 2 in lane 1 and all the tiles in the other lanes: --tiles s_1_0002,s_[2-8]

=item B<--flowcell-id> I<flow_cell_id>

Use the specified string as the flowcell id. (default value is parsed
from the config-file)

=item B<--help> Print a brief help message and exit.

Z<>

=item B<--man> View this help formatted in "man" style.

Z<>

=item B<Running the small test data set included with CASAVA>

Note that with the example below, if the SampleSheet.csv file is not found in the 
BaseCalls folder, the data is assumed to be non-multiplexed. Use L</--sample-sheet>
option to override default SampleSheet.csv location.

 /usr/local/bin/configureBclToFastq.pl \
 --output-dir ./Unaligned \
 --input-dir /usr/local/share/bcl2fastq-1.8.4/examples/Validation/Default/Data/Intensities/BaseCalls

  make -C Unaligned  # use make -j <parallel jobs> to speedup the conversion

=back

=head1 DESCRIPTION

B<configureBclToFastq.pl> reads a given Base Call directory, a 
SampleSheet.csv file, and configures a process which produces the
hierarchy of folders where demultiplexed data is stored as compressed
fastq files.

=head1 DIAGNOSTICS

=head2 Exit status

=over 4

=item B<0:> successful completion

=item B<1:> abnormal completion

=item B<2:> fatal error

Z<>

=item B<Errors:> All error messages are prefixed with "ERROR: ".

Z<>

=item B<Warnings:> All warning messages generated by CASAVA are prefixed with "WARNING: ".

Z<>

=back

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Mauricio Varea

=cut


use warnings "all";
use strict;
use Getopt::Long;
use Pod::Usage;
use Pod::Text;
use File::Spec;
use File::Copy;
use File::Basename;
use File::Path qw(mkpath);
# realpath is not exported by default, so it has to be explicit
use Cwd 'realpath';
use Cwd;
use Carp;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Alignment;
use Casava::BaseCalls;
use Casava::Common::Utils qw(reallyRealPath);
use Casava::Common::Log qw(logWarning logInfo errorExit initLog);
use Casava::Demultiplex;
use Casava::Demultiplex::Dx;
use Casava::Demultiplex::SampleSheet;

initLog( undef, 3, 0, 1, 1);

my @CLI = @ARGV;

my $man = 0;
my $help = 0;
my @adapterSequenceFiles = ();
my $adapterStringency = undef;
#my $adapterTrimming=0;
my $minimumTrimmedReadLength = undef;
my $maskShortAdapterReads = undef;
my $inputDir = "_default_";
my $positionsFormat= undef;
my $positionsDir = undef;
my $filterDir = undef;
my $outputDir = "_default_";
my $intDir = "_default_";
my $sampleSheetFile = "_default_";
my $mask = undef;
my $noEamss = undef;
my $withFailedReads = undef;
my $detection = "0";
my $fastqClusterCount=4000000;
my $compression='gzip';
my $gzLevel="1";
my $ignoreMissingBcl=0;
my $ignoreMissingCtrl = 0;
my $ignoreMissingStats=0;
my $force = 0;

my $filterPerRead = undef;
my $separateControls = undef;
my $userFlowCellID;
my $tilesFilter;
my $skipVariableMetadata = 0;

# This string will have all shell escapes already translated don't expect it to be re-executable
my $cmdAndArgs = join (' ', $0, @ARGV);

my $result = GetOptions('help|?' => \$help,
                        'adapter-sequence=s@' => \@adapterSequenceFiles,
                        'adapter-stringency:s' => \$adapterStringency,
                        #'adapter-trimming' => \$adapterTrimming,
                        'minimum-trimmed-read-length:i' => \$minimumTrimmedReadLength,
                        'mask-short-adapter-reads:i' => \$maskShortAdapterReads,
                        'sample-sheet:s' => \$sampleSheetFile,
                        'use-bases-mask=s' => \$mask,
                        'no-eamss' => \$noEamss,
                        'with-failed-reads' => \$withFailedReads,
                        'input-dir=s' => \$inputDir,
                        'intensities-dir=s' => \$intDir,
                        'positions-dir=s' => \$positionsDir,
                        'positions-format=s' => \$positionsFormat,
                        'filter-dir=s' => \$filterDir,
                        'output-dir=s' => \$outputDir,
                        'mismatches=s' => \$detection,
                        'fastq-cluster-count=i' => \$fastqClusterCount,
                        'compression=s' => \$compression,
                        'gz-level=i' => \$gzLevel,
                        'ignore-missing-bcl' => \$ignoreMissingBcl,
                        'ignore-missing-stats' => \$ignoreMissingStats,
                        'ignore-missing-control'  => \$ignoreMissingCtrl,
                        'filter-per-read!' => \$filterPerRead,
                        'separate-controls!' => \$separateControls,
                        'force' => \$force,
                        'flowcell-id=s' => \$userFlowCellID,
                        'tiles=s' => \$tilesFilter,
                        'skip-variable-metadata!' => \$skipVariableMetadata,
                        man => \$man) or pod2usage(2);

$cmdAndArgs = '' if $skipVariableMetadata;

pod2usage(1) if $help;
pod2usage(-verbose => 2,  -input => $1) if ($man and $0 =~ /(.*)/);
errorExit("ERROR: Unrecognized command-line argument(s): @ARGV")  if (0 < @ARGV);

if ($inputDir eq "_default_")
{
    $inputDir = ".";
    reallyRealPath($inputDir);
    logInfo "--input-dir is unspecified. Will configure the conversion from $inputDir\n";
}
else
{
    if(reallyRealPath($inputDir)) { errorExit("ERROR: Can't interpret input directory argument as an existing BaseCalls folder path: '$inputDir'."); }
}

if ($outputDir eq "_default_")
{
    my ($vol, $dirs, $lastDir) = File::Spec->splitpath($inputDir);
    my @directories = File::Spec->splitdir($dirs);
    pop @directories if !$directories[-1];
    if ('BaseCalls' eq $lastDir && 
        'Intensities' eq (pop @directories) && 
        'Data' eq (pop @directories) && scalar(@directories))
    {
        $outputDir = File::Spec->catdir($vol, @directories, 'Unaligned');
        logInfo "--output-dir is unspecified. Will configure the conversion in $outputDir\n"
    }
    else
    {
        errorExit "ERROR: --output-dir is unspecified and impossible to guess\n"
    }
}
errorExit "ERROR: input dir does not exist\n   "  unless -d $inputDir;
errorExit "ERROR: --fastq-cluster-count unspecified\n   "  unless defined $fastqClusterCount;

$intDir = File::Spec->catdir($inputDir,File::Spec->updir())
    if ($intDir eq "_default_");


if(reallyRealPath($outputDir)) { errorExit("ERROR: Can't interpret output directory argument: '$outputDir'.\n"); }

if (-e $outputDir)
{
    errorExit "ERROR: output path '$outputDir' already exists. Please specify --force to overwrite." unless $force;
}
else
{
    logInfo "Creating directory '$outputDir'";
    mkdir( $outputDir ) or errorExit "ERROR: failed to create directory '$outputDir': $!";
}

$Casava::Demultiplex::sampleSheetFile   = $sampleSheetFile;
$Casava::Demultiplex::basecallsDir      = $inputDir;
$Casava::Demultiplex::demultiplexedDir  = $outputDir;
$Casava::Demultiplex::intensitiesDir    = $intDir;
$Casava::Demultiplex::tilesFilter       = $tilesFilter;

$Casava::Demultiplex::strippedSampleSheet = "SampleSheet.csv";
$Casava::Demultiplex::supportFile   = File::Spec->catfile($outputDir,"support.txt");
$Casava::Demultiplex::mappingFile   = File::Spec->catfile($outputDir,"SampleSheet.mk");
$Casava::Demultiplex::makeFile      = File::Spec->catfile($outputDir,"Makefile");
my $demultiplexConfigXml = File::Spec->catfile($outputDir,"DemultiplexConfig.xml");

my $demux = Casava::Demultiplex->new;
my $guessedMask = $demux->guessUseBasesMaskFromRunInfo($mask);
logInfo("Original use-bases mask: " . (defined $mask ? $mask : "undefined"), 0);
logInfo("Guessed use-bases mask: " . (defined $guessedMask ? $guessedMask : "undefined"), 0);
my $adapterReadIndex = 0;
foreach my $adapterSequenceFile (@adapterSequenceFiles)
{
    $demux->adapterSequenceFile(++$adapterReadIndex,$adapterSequenceFile);
}
$demux->adapterStringency($adapterStringency)  if defined $adapterStringency;
#$demux->adapterTrimming($adapterTrimming);
$demux->minimumTrimmedReadLength($minimumTrimmedReadLength)  if defined $minimumTrimmedReadLength;
$demux->maskShortAdapterReads($maskShortAdapterReads)        if defined $maskShortAdapterReads;
$demux->mask($guessedMask);
$demux->noEamss($noEamss);
$demux->withFailedReads($withFailedReads);
$demux->ignoreMissingBcl($ignoreMissingBcl);
$demux->ignoreMissingCtrl($ignoreMissingCtrl);
$demux->ignoreMissingStats($ignoreMissingStats);
$demux->errorDetection($detection);
$demux->fastqClusterCount($fastqClusterCount);
$demux->positionsDir($positionsDir)         if defined $positionsDir;
$demux->positionsFormat($positionsFormat)   if defined $positionsFormat;
$demux->filterDir($filterDir)               if defined $filterDir;
$demux->flowCellId($userFlowCellID)         if defined $userFlowCellID;
$demux->compression($compression);
$demux->gzLevel($gzLevel);
$demux->isFilterPerRead($filterPerRead)     if defined $filterPerRead;
$demux->needControlFile($separateControls)  if defined $separateControls;

$demux->logSupport('>',@CLI);

# my %tiles = getTilesFromDir($inputDir);

$demux->loadSampleSheet();
$demux->loadAdapterSequences($adapterReadIndex);

$demux->logSupport('>>');

$demux->run();

$demux->saveSampleSheet("mappingFile");
$demux->saveDemultiplexConfig($cmdAndArgs, $demultiplexConfigXml);

$demux->generateMakefile();
$demux->selfTest($outputDir);

1;

__END__

