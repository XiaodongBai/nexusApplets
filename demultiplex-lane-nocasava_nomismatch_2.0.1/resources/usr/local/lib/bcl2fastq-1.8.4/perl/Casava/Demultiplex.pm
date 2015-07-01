=head1 LICENSE

Copyright (c) 2009 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And VAriation
(CASAVA) software package.

=head1 NAME

Casava::Demultiplex - Utility library for indexing.

=head1 SYNOPSIS


=head1 DESCRIPTION

This script is a rewrite of the original GERALD library originally
implemented for the Pipeline.

Library for the creation, initialization and configuration of the
Alignment folder.

Provides a number of alignment-related services for CASAVA;

=head1 SUBROUTINES

=cut

# POD for for each subroutines is right before the implementation
# More POD after __END__

package Casava::Demultiplex;

use strict;
use warnings "all";
use Carp;
use Exporter 'import';
use Cwd;

use File::Spec;
use File::Path;
use File::Copy;
use File::Basename;

use XML::Simple;
#use Data::Dumper;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
my $DATA_DIR = '/usr/local/share/bcl2fastq-1.8.4';

use Casava::BaseCalls;
use Casava::Common::Log;
use Casava::Common::Utils qw(expandUseBases expandUseBasesString);
use Casava::Demultiplex::DemultiplexConfig;
use Casava::Demultiplex::Dx;
use Casava::Demultiplex::SampleSheet;
use Casava::Intensities;

our @EXPORT_OK = qw(&new $sampleSheetFile $demultiplexedDir $basecallsDir $intensitiesDir $supportFile $mappingFile $makeFile );

our $sampleSheetFile     = undef;
our $demultiplexedDir    = undef;
our $basecallsDir        = undef;
our $intensitiesDir      = undef;
our $tilesFilter         = undef;

our $supportFile         = undef; #"support.txt";
our $mappingFile         = undef; #"SampleSheet.mk";
our $makeFile            = undef; #"Makefile";
our $strippedSampleSheet = undef;
my $projectsRoot = "Unaligned";
my $defaultUseBasesMask = 'Y*';

sub new(;$)
{
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self  = {};

    $self->{NO_EAMSS}               = undef;
    $self->{WITH_FAILED_READS}      = undef;
    $self->{IGNORE_MISSING_BCL}     = undef;
    $self->{IGNORE_MISSING_STATS}   = undef;
    $self->{ADAPTER_SEQUENCE_FILES} = {};
    $self->{ADAPTER_SEQUENCES}      = {};
    $self->{ADAPTER_STRINGENCY}     = undef;
    #$self->{ADAPTER_TRIMMING}       = undef;
    $self->{MINIMUM_TRIMMED_READ_LENGTH} = undef;
    $self->{MASK_SHORT_ADAPTER_READS} = undef;
    $self->{MASK}                   = undef;
    $self->{ERROR_DETECTION}        = undef;
    $self->{POSITIONS_DIR}          = undef;
    $self->{POSITIONS_FORMAT}       = undef;
    $self->{FASTQ_CLUSER_COUNT}     = undef;
    $self->{COMPRESSION}            = undef;
    $self->{GZ_LEVEL}               = undef;
    $self->{COMPRESSED_BCL}         = undef;
    $self->{FILTER_DIR}             = undef;
    $self->{FILTER_PER_READ}        = 0;
    $self->{NEED_CONTROL_FILE}      = 0;
    $self->{BaseCalls}              = Casava::BaseCalls->new();
    $self->{SampleSheet}            = Casava::Demultiplex::SampleSheet->new();
    $self->{Support}                = Casava::Demultiplex::Dx->new();
    $self->{Intensities}            = Casava::Intensities->new();
    $self->{DEMULTIPLEXED}          = {};

    $self->{BaseCalls}->path($basecallsDir);
    $self->{BaseCalls}->fillTilesFromXml();
    $self->{BaseCalls}->filterTiles($tilesFilter) if defined $tilesFilter and length($tilesFilter);
    $self->{Intensities}->path($intensitiesDir);

    if (defined $self->{Intensities}->getCompressedBcl())
    {
        $self->{COMPRESSED_BCL} = 1  if $self->{Intensities}->getCompressedBcl() =~ m/true/i;
    } else {
        $self->{COMPRESSED_BCL} = 0;
    }

    my $baseCallsSoftwareName = $self->{BaseCalls}->getSoftwareName();
    if (defined $baseCallsSoftwareName)
    {
        my ($swVer,$swBuild) = $self->{BaseCalls}->getFormattedVersion();
        logInfo("Basecalling software: $baseCallsSoftwareName", 0);
        logInfo("             version: $swVer ".(defined $swBuild ? "(build $swBuild)" : ""), 0)
               if defined $swVer;

        if ($baseCallsSoftwareName eq 'RTA' && defined $swVer)
        {
            if ($swVer < 1.09)
            {
                $self->{FILTER_DIR} = $basecallsDir;
            }
            elsif ($swVer == 1.09)
            {
                $self->{FILTER_PER_READ} = defined $swBuild && $swBuild < 31;    # up to build 30
                $self->{NEED_CONTROL_FILE} = defined $swBuild && $swBuild > 30;  # from build 31
            }
            elsif ($swVer == 1.10)
            {
                $self->{FILTER_PER_READ} = defined $swBuild && $swBuild < 29;    # up to build 28
                $self->{NEED_CONTROL_FILE} = defined $swBuild && $swBuild > 28;  # from build 29
            }
            else
            {
                $self->{NEED_CONTROL_FILE} = 1;  # from build 29
            }
        }
        elsif ($baseCallsSoftwareName eq 'Bustard')
        {
            my $intensitiesSoftwareName = $self->{Intensities}->getSoftwareName();
            if (defined $intensitiesSoftwareName)
            {
                if ($intensitiesSoftwareName eq 'RTA')
                {
                    ($swVer,$swBuild) = $self->{Intensities}->getFormattedVersion();
                    logInfo("Image analysis software: $intensitiesSoftwareName, version: $swVer ".(defined $swBuild ? "(build $swBuild)" : ""), 0);
                }
                else
                {
                    logWarning ("Cannot auto-detect format/location of positions files when Image analyzer '$intensitiesSoftwareName' is in use.");
                    logWarning ("Auto-detection only works with Illumina's 'RTA' and 'OLB' softwares. Make sure these parameters are explicitly set.");
                    $swVer = undef;
                }
            }
        }
        else 
        {
            logWarning ("Cannot auto-detect format/location of filter files when Basecaller '$baseCallsSoftwareName' is in use.");
            logWarning ("Auto-detection only works with Illumina's 'RTA' and 'OLB' softwares. Make sure these parameters are explicitly set.");
            $swVer = undef;
        }

        if (defined $swVer)
        {
            if ($swVer < 1.09)
            {
                $self->{POSITIONS_DIR} = $intensitiesDir;
                $self->{POSITIONS_FORMAT} = '_pos.txt';
            }
            else
            {
                $self->{POSITIONS_FORMAT} = $self->{Intensities}->getLocationFileType();
                unless (defined $self->{POSITIONS_FORMAT})
                {
                    my $instrumentType = $self->{Intensities}->getInstrumentType();
                    $self->{POSITIONS_FORMAT} = ($instrumentType eq 'HiSeq' ? '.clocs' : '.locs')
                                                if defined $instrumentType;
                }
            }
        } else {
            my $intensitiesConfig = File::Spec->catfile($self->{Intensities}->path, $Casava::Intensities::configXml);
            my $baseCallsConfig = File::Spec->catfile($self->{BaseCalls}->path, $Casava::BaseCalls::configXml);
            logWarning ("Auto-detection failed to recognise software version. Please check '$baseCallsConfig' and '$intensitiesConfig'");
        }
    }
    
    bless ($self, $class);

    return $self;
}

sub dirs($)
{
    my $self = shift;
    return keys %{$self->{DEMULTIPLEXED}};
}

sub sampleSheet($)
{
    my $self = shift;
    return $self->{SampleSheet};    
}

sub guessUseBasesMaskFromRunInfo($$)
{
    my $self = shift;
    my ($mask) = @_;
    my $runInfo = $self->{Intensities}->getRunInfo();
    if (not $runInfo)
    {
        logWarning("Couldn't find RunInfo.xml for " . $self->{BaseCalls}->path);
        return $mask;
    }
    my $reads = $runInfo->{'Reads'}->{'Read'};
    if (not $reads)
    {
        logWarning("Couldn't find reads in RunInfo.xml");
        return $mask;
    }
    my @originalReadMasks = split(',', $mask) if defined $mask;
    my @readMasks;
    foreach my $read (@$reads)
    {
        my $readMask = shift @originalReadMasks;
        my $readLength = $read->{'NumCycles'};
        if ($readMask)
        {
            if ($readMask =~ /\*/)
            {
                push @readMasks, expandUseBases($readMask, $readLength);
            }
            else
            {
                push @readMasks, expandUseBasesString($readMask);
            }
        }
        else
        {
            my $selector = 'y';
            $selector = 'I' if exists $read->{'IsIndexedRead'} and  $read->{'IsIndexedRead'} eq 'Y';
            # By default, the mask uses all bases, except for the barcode where the last base is discarded
            my $lastSelector = ($selector eq 'I' ? 'n' : 'y');
            my $repeat = "$selector"x($readLength - 1);
            push @readMasks, "$selector"x($readLength - 1) . $lastSelector;
        }
    }
    return join(',', @readMasks);
}

sub adapterSequenceFile($$;$)
{
    my $self = shift;
    my $readNumber = shift;
    if (@_) { $self->{ADAPTER_SEQUENCE_FILES}->{$readNumber} = shift }

    if (!defined $self->{ADAPTER_SEQUENCE_FILES}->{$readNumber})
    {
        if (1 < $readNumber)
        {
            return $self->adapterSequenceFile($readNumber - 1);
        }
        return '';
    }
    return $self->{ADAPTER_SEQUENCE_FILES}->{$readNumber};
}

sub adapterSequence($$;$)
{
    my $self = shift;
    my $readNumber = shift;
    if (@_) { $self->{ADAPTER_SEQUENCES}->{$readNumber} = shift }

    if (!defined $self->{ADAPTER_SEQUENCES}->{$readNumber})
    {
        if (1 < $readNumber)
        {
            return $self->adapterSequence($readNumber - 1);
        }
        return [];
    }
    return $self->{ADAPTER_SEQUENCES}->{$readNumber};
}

sub adapterStringency($;$)
{
    my $self = shift;
    if (@_) {$self->{ADAPTER_STRINGENCY} = shift};
    return $self->{ADAPTER_STRINGENCY};
}

#sub adapterTrimming($;$)
#{
#    my $self = shift;
#    if (@_) {$self->{ADAPTER_TRIMMING} = shift};
#    return $self->{ADAPTER_TRIMMING};
#}

sub minimumTrimmedReadLength($;$)
{
    my $self = shift;
    if (@_) {$self->{MINIMUM_TRIMMED_READ_LENGTH} = shift};
    return $self->{MINIMUM_TRIMMED_READ_LENGTH};
}

sub maskShortAdapterReads($;$)
{
    my $self = shift;
    if (@_) {$self->{MASK_SHORT_ADAPTER_READS} = shift};
    return $self->{MASK_SHORT_ADAPTER_READS};
}

sub mask($;$)
{
    my $self = shift;
    ($self->{MASK}) = @_  if (@_);
    errorExit("ERROR: undefined mask")  unless( defined $self->{MASK} );
    return $self->{MASK};
}

sub errorDetection($;$)
{
    my $self = shift;
    if (@_) {$self->{ERROR_DETECTION} = shift};
    return $self->{ERROR_DETECTION};
}

sub noEamss($;$)
{
    my $self = shift;
    if (@_) {$self->{NO_EAMSS} = shift};
    return $self->{NO_EAMSS};
}

sub withFailedReads($;$)
{
    my $self = shift;
    if (@_) {$self->{WITH_FAILED_READS} = shift};
    return $self->{WITH_FAILED_READS};
}

sub ignoreMissingBcl($;$)
{
    my $self = shift;
    if (@_) {$self->{IGNORE_MISSING_BCL} = shift};
    return $self->{IGNORE_MISSING_BCL};
}

sub ignoreMissingStats($;$)
{
    my $self = shift;
    if (@_) {$self->{IGNORE_MISSING_STATS} = shift};
    return $self->{IGNORE_MISSING_STATS};
}

sub ignoreMissingCtrl($;$)
{
    my $self = shift;
    if (@_) {$self->{IGNORE_MISSING_CTRL} = shift};
    return $self->{IGNORE_MISSING_CTRL};
}

sub positionsDir($;$)
{
    my $self = shift;
    if (@_) {$self->{POSITIONS_DIR} = shift};
    return $self->{POSITIONS_DIR};
}

sub positionsFormat($;$)
{
    my $self = shift;
    if (@_) {$self->{POSITIONS_FORMAT} = shift};
    return $self->{POSITIONS_FORMAT};
}

sub filterDir($;$)
{
    my $self = shift;
    if (@_) {$self->{FILTER_DIR} = shift};
    return $self->{FILTER_DIR};
}

sub isFilterPerRead($;$)
{
    my $self = shift;
    if (@_) {$self->{FILTER_PER_READ} = shift};
    return $self->{FILTER_PER_READ};
}

sub needControlFile($;$)
{
    my $self = shift;
    if (@_) {$self->{NEED_CONTROL_FILE} = shift};
    return $self->{NEED_CONTROL_FILE};
}

sub compressedBcl($;$)
{
    my $self = shift;
    if (@_) {$self->{COMPRESSED_BCL} = shift};
    return $self->{COMPRESSED_BCL};
}

sub fastqClusterCount($;$)
{
    my $self = shift;
    if (@_) {$self->{FASTQ_CLUSTER_COUNT} = shift};
    return $self->{FASTQ_CLUSTER_COUNT};
}

sub flowCellId($;$)
{
    my $self = shift;
    if (@_) {$self->{FLOWCELL} = shift};
    return $self->{FLOWCELL};
}

sub compression($;$)
{
    my $self = shift;
    if (@_) {$self->{COMPRESSION} = shift};
    return $self->{COMPRESSION};
}

sub gzLevel($;$)
{
    my $self = shift;
    if (@_) {$self->{GZ_LEVEL} = shift};
    return $self->{GZ_LEVEL};
}

sub _getDefaultSampleSheet($$)
{
    my ($inputDirectory, $flowCellId) = @_;
    # create the content of the default sample sheet in-memory
    use Casava::BaseCalls;
    my $baseCalls = Casava::BaseCalls->new;
    $baseCalls->path($inputDirectory);
    $baseCalls->fillTilesFromXml;
    my $sampleSheet = "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n";
    foreach my $lane (sort $baseCalls->lanes())
    {
        my $laneNumber = $1 if $lane =~ /^s_(\d)$/ or errorExit("ERROR: failed to parse lane numbers from $inputDirectory");
        next unless($laneNumber);
        my $project = ( $flowCellId ? $flowCellId : 'default' );
        my $sampleName = "lane$laneNumber";
        $sampleSheet .="$flowCellId,$laneNumber,$sampleName,Unknown,,'DefaultSample',N,,,$project\n";
    }
    return \$sampleSheet;
}

sub loadAdapterSequences($$)
{
    my $self = shift;
    my $adapterFileCount = shift;

    foreach my $demultiplexedReadIndex (1 .. $adapterFileCount)
    {
        my $adapterSequenceFile = $self->adapterSequenceFile($demultiplexedReadIndex);
        return unless defined $adapterSequenceFile;

        my $transformFastaCmd = File::Spec->catfile('/usr/local/libexec/bcl2fastq-1.8.4', 'transformFasta.pl') .
                              " --targetPath - --chromNameSource=contigName --removeBreaks --removeHeaders $adapterSequenceFile";
        logDebug "executing: $transformFastaCmd", 2;
        open FASTA, '-|', "$transformFastaCmd" || errorExit("ERROR: Failed to load adapter sequences from $adapterSequenceFile: $!");

        my @adapterSequences = ();
        while (<FASTA>)
        {
            chomp;
            my $line = $_;
            push @adapterSequences, $line;
        }
        $self->adapterSequence($demultiplexedReadIndex, \@adapterSequences);
        logInfo "Loaded " . scalar(@adapterSequences) . " adapter sequences from $adapterSequenceFile";
        logDebug "  ... loaded adapter sequences are: '@adapterSequences'", 2;
    }
}

sub loadSampleSheet($;$)
{
    my $self = shift;
    if (@_) { $sampleSheetFile = shift };
    croak "ERROR: undefined Sample Sheet file\n   "  unless defined $sampleSheetFile;

    if ($sampleSheetFile eq "_default_")
    {
        my $defaultSampleSheetFile = File::Spec->catfile($self->{BaseCalls}->path, 'SampleSheet.csv');
        # simply use the the default sample sheet if it exists
        $sampleSheetFile = $defaultSampleSheetFile if -f $defaultSampleSheetFile;
    }

    my $flowCellId = $self->flowCellId();
    my $baseCallsFlowCellId = $self->{BaseCalls}->flowCellId();
    logWarning("Overriding flowcell ID '$baseCallsFlowCellId' from the BaseCalls configuration with '$flowCellId'") 
        if $baseCallsFlowCellId && $flowCellId && ($flowCellId ne $baseCallsFlowCellId);

    $flowCellId = $baseCallsFlowCellId unless $flowCellId;

    if (!$sampleSheetFile or '_default_' eq $sampleSheetFile)
    {
        $self->{SampleSheet} = Casava::Demultiplex::SampleSheet::create("SampleSheet.csv");
        $sampleSheetFile = _getDefaultSampleSheet($self->{BaseCalls}->path, $flowCellId) ;
    }
    else
    {
        $self->{SampleSheet} = Casava::Demultiplex::SampleSheet::create($sampleSheetFile)
             or errorExit "ERROR: sample sheet provided is of unknown extension: ".
                      "only *.csv or *.xml are supported\n   ";
    }
    open( my $sampleSheetHandle, "<", $sampleSheetFile )
        or errorExit("ERROR: Failed to open file '$sampleSheetFile', for reading: $!\n   ");
    $self->{SampleSheet}->load($sampleSheetHandle);
    errorExit("ERROR: Did not retrieve any barcodes from '$sampleSheetFile'\n") unless $self->{SampleSheet}->laneNumbers;
    close( $sampleSheetHandle );

    my $sampleSheetFlowCellId = $self->{SampleSheet}->flowCellId();
    logWarning("SampleSheet flowcell ID '$sampleSheetFlowCellId' does not match. Expected: '$flowCellId'\n")
               if($sampleSheetFlowCellId && $flowCellId && ($flowCellId ne $sampleSheetFlowCellId));

    $flowCellId = $sampleSheetFlowCellId unless $flowCellId;
    $self->flowCellId($flowCellId);
    $self->{SampleSheet}->flowCellId($flowCellId);

    # check that the barcodes have the expected length
    my $expectedBarcodeLength = $self->mask =~ tr/iI//;
    foreach my $lane ($self->{SampleSheet}->lanes)
    {
        foreach my $barcode ($self->{SampleSheet}->barcodes($lane))
        {
            $barcode =~ s/$Casava::Demultiplex::SampleSheet::barcodeComponentSeparator//g;
            next unless $barcode and length($barcode) > 0;
            next if 'Undetermined' eq $barcode;
            next if 'NoIndex' eq $barcode;
            errorExit("ERROR: barcode $barcode for lane $lane has length " .
                      length($barcode) .
                      ": expected barcode lenth (including delimiters) is $expectedBarcodeLength") if length($barcode) != $expectedBarcodeLength;
        }
    }
}

sub saveSampleSheet($$)
{
    my $self = shift;
    my $varname = "\$".(shift);  # has to be one of *our* exported variables
    croak "ERROR: undefined Sample Sheet file\n   "  unless defined $sampleSheetFile;
    my $file = eval $varname;
    croak "ERROR: unknown Sample Sheet type (don't know anything about '$varname')\n   "
          unless defined $file;
    my $ext;
    (undef,undef,$ext) = fileparse($file,@Casava::Demultiplex::SampleSheet::FILE_FORMATS);
    croak "ERROR: unknown extension for Sample Sheet: ".
          "only \{".join(',',@Casava::Demultiplex::SampleSheet::FILE_FORMATS)."\} are supported.\n   "
          unless defined $ext;

    my $ss;
    # different extension => needs cloning
    if ($sampleSheetFile =~ /\\$ext$/)
    {
        $ss = $self->{SampleSheet};
    } else {
        $ss = Casava::Demultiplex::SampleSheet::create($file);
        $ss->clone($self->{SampleSheet});
    }
    open( my $handle, '>', $file )
        or croak "ERROR: Failed to open file '$file', for writing: $!\n   ";
    $ss->save($handle);
    close( $handle );    
}

sub saveDemultiplexConfig($$$)
{
    my ($self, $cmdAndArgs, $file) = @_;

    my $dc = Casava::Demultiplex::DemultiplexConfig->new();
    $dc->software({Name => (File::Spec->splitpath($0))[2], Version => 'bcl2fastq-1.8.4'});
    $dc->cmdAndArgs($cmdAndArgs);
    $dc->baseCallsSoftware($self->{BaseCalls}->getSoftware());
    $dc->intensitiesSoftware($self->{Intensities}->getSoftware());
    $dc->controlSoftware($self->{Intensities}->getControlSoftware());
    $dc->clone($self->{SampleSheet});

    open( my $handle, '>', $file ) or errorExit "ERROR: Failed to open file '$file', for writing: $!";
    $dc->save($handle);
    close( $handle );
}


sub logSupport($$;@)
{
    my $self = shift;
    my $mode = shift;
    my @args = @_;
    croak "ERROR: undefined output dir\n   "    unless defined $demultiplexedDir;
    croak "ERROR: undefined Support file\n   "  unless defined $supportFile;
    croak "ERROR: undefined mode: Only '>' and '>>' are allowed\n   "
          unless defined $mode && ($mode eq '>' || $mode eq '>>');
    
    open( my $supportHandle, $mode, $supportFile )
        or croak "ERROR: Failed to open file '$supportFile', for writing: $!\n   ";
    $self->{Support}->generalInfo($supportHandle, \@args)                              if $mode eq '>';
    $self->{Support}->specificInfo($supportHandle, $self->mask, $self->{SampleSheet})  if $mode eq '>>';
    close( $supportHandle );
}



sub run()
{
    my $self = shift;
    croak "ERROR: undefined output dir\n   "                unless defined $demultiplexedDir;
    croak "ERROR: undefined Sample Sheet output file\n   "  unless defined $strippedSampleSheet;

    my $isDemux = $self->{SampleSheet}->isDemux;

    foreach my $lane (map {$_ =~ /_(\d+)/; $1} sort $self->{BaseCalls}->lanes)
    {
        my $tiles = $self->{BaseCalls}->tiles($lane);
        croak "ERROR: No tiles in lane '$lane'. Please review your sample sheet.\n   "
              unless defined $tiles;

        unless (defined $self->{SampleSheet}->barcodes($lane))
        {
            logWarning "No sample sheet data for lane '$lane'. Skipping lane.";
            next;
        }

        for my $barcode ($self->{SampleSheet}->barcodes($lane))
        {

            my $projectSampleDirName = $self->{SampleSheet}->projectSampleDirName($lane,$barcode);
            if (defined $projectSampleDirName && !exists $self->{DEMULTIPLEXED}->{$projectSampleDirName})
            {
                my $samplePath=File::Spec->catdir($demultiplexedDir, $projectSampleDirName);

                $self->{DEMULTIPLEXED}->{$projectSampleDirName} = Casava::BaseCalls->new();
                $self->{DEMULTIPLEXED}->{$projectSampleDirName}->path($samplePath);
                $self->{DEMULTIPLEXED}->{$projectSampleDirName}->initialize();

                my $ssFile = File::Spec->catfile($samplePath,$strippedSampleSheet);
                my $ss = Casava::Demultiplex::SampleSheet::create($strippedSampleSheet);
                my %cond = (projectId => $self->{SampleSheet}->projectId($lane,$barcode),
                            sampleId =>  $self->{SampleSheet}->sampleId($lane,$barcode));
                $ss->clone( $self->{SampleSheet}, \%cond );
                open( my $strippedHandle, '>', $ssFile)
                  or croak "ERROR: Failed to open file '$ssFile', for writing: $!\n   ";
                $ss->save($strippedHandle);
                close( $strippedHandle );
                $self->{DEMULTIPLEXED}->{$projectSampleDirName}->addTiles($lane,@{$tiles});
            }
        }
    }
}

sub generateMakefile($)
{
    my $self = shift;
    open( my $makefileHandle, ">$Casava::Demultiplex::makeFile" )
        or croak "ERROR: Failed to open file '$Casava::Demultiplex::makeFile', for writing: $!\n   ";
    $self->makefile($makefileHandle);
    close( $makefileHandle );
}

sub makefile($$)
{
    my $self = shift;
    my $handle = shift;
    # only the lanes mentioned in sample sheet must go through
    my $lanesFilter = join('|', $self->{SampleSheet}->laneNumbers());
    my @lanes = grep(/$lanesFilter/, (map {$1 if $_ =~ /^.+_(\d)$/} keys %{$self->{BaseCalls}->{TILES}}));
    my $baseCalls = $self->{BaseCalls}->path;
    my $intensities = $self->{Intensities}->path;
    my $isDemux = $self->{SampleSheet}->isDemux;

    print $handle "\# This makefile was automatically generated by ".$self->{Support}->id().".\n";
    print $handle "\# Please do not edit.\n\n";

#    print $handle "POSITIONS_DIR:=$self->positionsDir()" if($self->positionsDir());
#    print $handle "POSITIONS_FORMAT:=$self->positionsFormat()" if($self->positionsFormat());
    print ($handle 'FILTER_DIR:=' . $self->filterDir() . "\n") if($self->filterDir());
    print ($handle "IGNORE_MISSING_BCL:=yes\n") if($self->ignoreMissingBcl());
    print ($handle "IGNORE_MISSING_STATS:=yes\n") if($self->ignoreMissingStats());
    print ($handle "IGNORE_MISSING_CTRL:=yes\n") if($self->ignoreMissingCtrl());
    
    print $handle <<EOF;
\# first target needs to be defined in the beginning. Ohterwise includes such as
\# Log.mk cause unexpected behavior
firsttarget: all

BASECALLS_DIR:=$baseCalls
INTENSITIES_DIR:=$intensities
MAKEFILES_DIR:=$DATA_DIR/makefiles

\# Import the global configuration
include \$(MAKEFILES_DIR)/Config.mk

include \$(MAKEFILES_DIR)/Sentinel.mk

\# Import the debug functionalities
include \$(MAKEFILES_DIR)/Debug.mk

\# Import the logging functionalities
include \$(MAKEFILES_DIR)/Log.mk

\# Setup:
EOF
    # this needs to be changed accordingly, once BF-363 is in
    #print $handle "DEMUX_READS:=1".(($self->mask =~ /y/) ? " 2\n" : "\n");
    my $barcodeCycles = $self->{BaseCalls}->barcodeCycles($self->{MASK});
    unless (defined $barcodeCycles)
    {
        print STDERR "WARNING: Unable to detect barcode cycles from config.xml. Assuming non-multiplexed data.\n";
        my @emptyArray;
        $barcodeCycles = \@emptyArray;
    }
    print $handle "BARCODE_CYCLES:=". join(' ',@{$barcodeCycles}) ."\n";

    my $readCycles = $self->{BaseCalls}->readCycles($self->{MASK}, 0);
    my $readDemuxCycles = $self->{BaseCalls}->readCycles($self->{MASK}, 1);
    errorExit "ERROR: Unable to detect read cycles from config.xml when the following mask applied: $self->{MASK}"
              unless defined $readCycles;
# This will renumber the reads so that they are sequential.
    my @excludedReads=();
    my @includedReads=();# = (1..$self->{BaseCalls}->getReads($self->{BaseCalls}->path));
    my $demuxRead = 1;
    foreach my $read (sort {$a<=>$b} keys(%$readCycles))
    {
        print $handle "r${demuxRead}_CYCLES:=". join(' ',@{$readCycles->{$read}}) ."\n";
        print $handle "r${demuxRead}_DEMUX_CYCLES:=". join(' ',@{$readDemuxCycles->{$read}}) ."\n";
        if (scalar(@excludedReads) + $demuxRead != $read)
        {
            my $lastExcluded = $excludedReads[-1];
            push @excludedReads, ( ++$lastExcluded .. $read - 1);
        }
        push @includedReads, $read;
        print $handle "r${demuxRead}_ADAPTER:=". join(' ',@{$self->adapterSequence($demuxRead)}) . "\n";
        ++$demuxRead;
    }
    my @demuxReads = (1 .. keys(%$readCycles));
    my @allCycles = $self->{BaseCalls}->allCycles();
    # Original reads and demux reads must correspond as the read renumbering in xsl uses translate
    # function to convert from old to new numbers
    print $handle "ALL_ORIGINAL_CYCLES:=@allCycles\n";
    print $handle "INCLUDED_ORIGINAL_READS:=@includedReads\n";
    print $handle "DEMUX_READS:=@demuxReads\n";
#    print $handle "EXCLUDED_ORIGINAL_READS:=@excludedReads\n";
# the section below preserves original read numbers. At the moment it is believed the post-alignment
# will be unable to cope with read 3 being read2. 
#    foreach my $read (keys(%$readCycles))
#    {
#        print $handle "r${read}_CYCLES:=". join(' ',@{$readCycles->{$read}}) ."\n";
#        push @allCycles, @{$readCycles->{$read}};
#    }
#    print $handle "DEMUX_READS:=". join(' ', keys(%$readCycles)) ."\n";


    print $handle 'COMPRESSION:=' . $self->compression() . "\n";
    print $handle 'FLOWCELL:=' . $self->flowCellId() . "\n";
    
    print $handle 'COMPRESSIONSUFFIX:=' . {undef=>'', none=>'', gzip=>'.gz', bzip2=>'.bz2'}->{$self->compression()} .  "\n";

    print $handle "FILTER_PER_READ:=1\n" if $self->isFilterPerRead();
#    print $handle "QSEQ_MASK:=".$self->mask()."\n";
    print $handle "DEMUX_OPTIONS:=" . 
                  ($self->noEamss() ? " --no-eamss " : "" ).
                  ($self->compressedBcl() ? " --compressed-bcl " : "" ).
                  ($self->withFailedReads() ? " --with-failed-reads " : "" ).
                  ($self->adapterStringency() ? " --adapter-stringency=" . $self->adapterStringency() : "" ).
                  #($self->adapterTrimming() ? " --adapter-trimming " : "" ).
                  ($self->minimumTrimmedReadLength() ? " --minimum-trimmed-read-length=" . $self->minimumTrimmedReadLength() : "" ).
                  ($self->maskShortAdapterReads() ? " --mask-short-adapter-reads=" . $self->maskShortAdapterReads() : "" ).
                  ($self->isFilterPerRead() ? " --filter-per-read" : "" ).
                  ($self->needControlFile() ? " --need-separate-controls" : "" ).
                  ($self->positionsDir() ? " --positions-dir=" . $self->positionsDir() : "" ).
                  ($self->positionsFormat() ? " --positions-format=" . $self->positionsFormat() : "" ).
                  (defined $self->fastqClusterCount() ? " --fastq-cluster-count=" . $self->fastqClusterCount() : "" ).
                  ($self->gzLevel() ? " --gz-level=" . $self->gzLevel() : "" ).
                  (" --mismatches='" . $self->errorDetection()."'") .
                  (" --instrument-name='" . $self->{BaseCalls}->instrumentName()).
                  ("' --run-number='" . $self->{BaseCalls}->runNumber()).
                  "'\n\n";

    print $handle "STATS_TO_SIGNAL_MEANS_OPTIONS:=" . 
                  ($self->needControlFile() ? " --need-separate-controls" : "" ).
                  "\n\n";

    print $handle <<EOF;
\# Information extracted from SampleSheet.csv
include ./SampleSheet.mk
\# this list gets populated by DemultiplexLaneRead.mk
ALL_DEMUX_SUMMARIES:=
EOF

#READS:=\$(words \$(ORIGINAL_READS))

    my @laneTargets = (map {"l${_}.done"} @lanes);
    print $handle "post_run.done: @laneTargets DemultiplexedBustardSummary.xml DemultiplexedBustardConfig.xml Basecall_Stats_\$(FLOWCELL)/Demultiplex_Stats.htm\n";
    print $handle "\t\$(POST_RUN_COMMAND)\n\n";
    print $handle "all: clean_intermediate.done\n";
    print $handle "\t\@\$(LOG_INFO)  \$\@ completed successfully.\n\n";
    
    for my $read (1 .. keys(%$readCycles))
    {
        my @readTargets = map {"l${_}_r${read}.done"} @lanes;
        print $handle "r${read}_post_run.done: @readTargets\n";
        print $handle "\t\$(POST_RUN_COMMAND_R$read)\n\n";
        print $handle "r${read}: r${read}_post_run.done\n";
        print $handle "\t\@\$(LOG_INFO)  \$\@ completed successfully.\n\n";
    }

    print $handle "\# DEMULTIPLEXING:\n\n";

    $tilesFilter = '' unless defined $tilesFilter;
    print $handle "TILES_FILTER := $tilesFilter\n\n";

    foreach my $ln (@lanes)
    {
        print $handle "l${ln}_TILES:=".join(' ', sort {$a <=> $b} @{$self->{BaseCalls}->tiles($ln)})."\n";
        print $handle "lane:=$ln\n";
        print $handle "include \$(MAKEFILES_DIR)/DemultiplexLane.mk\n\n";
    }
    print $handle <<"EOF";

include \$(MAKEFILES_DIR)/DemultiplexFlowCell.mk

.PHONY: self_test
self_test: DemultiplexedBustardConfig.xml
\t\@\$(LOG_INFO) \$\@ completed successfully.;\n\n
EOF

}

sub selfTest
{
    my ($self, $demuxDir) = @_;
    use Cwd;
    my $currDir = getcwd();
    chdir $demuxDir;
    my $selfTestCommand = "make self_test";
    logInfo("Running self tests: '$selfTestCommand'", 0);
    my $ENV_REF = \%ENV;
    $ENV_REF->{'PATH'} = $1 if $ENV_REF->{'PATH'} =~ /^(.*)$/;    # untaint PATH variable
    my $selfTestOutput = `$selfTestCommand 2>&1`;
    my $selfTestResult = $?;
    my $exitValue = ($selfTestResult >> 8);
    if ( $exitValue != 0 )
    {
        my $signal = ($selfTestResult & 127);
        logWarning("\noutput of '$selfTestCommand':\n\n$selfTestOutput\n\n" .
             "Self test command exited with error $exitValue (signal $signal)\n\n" .
             "Investigate and fix errors, then retry\n");
        chdir $currDir;
        errorExit("Self test command exited with error $exitValue (signal $signal)");
    }
    chdir $currDir;
    logInfo("Running self tests on $demuxDir completed with no problems\n", 0);
}


1;
__END__

=pod

=head1 DIAGNOSTICS

=head2 Exit status

=over 4

=item 0: succesful completion

=item 1: abnormal completion

=item 2: fatal error

=back

=head2 Errors

All error messages are prefixed with "ERROR: ". 

=head2 Warnings

All warning messages generated by CASAVA are prefixed with "WARNING: ".

=head1 CONFIGURATION AND ENVIRONMENT

Name and location of configuration files, with a complete description of the properties that can be set.

Name and description of the relevant environment variables

=head1 DEPENDENCIES

=over 4

=item Standard perl modules

strict, warnings, Carp,

=item External perl modules

XML::Simple

=item Casava perl modules

Casava::Demultiplex

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Mauricio Varea

=cut
