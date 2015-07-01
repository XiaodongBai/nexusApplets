=head1 LICENSE

Copyright (c) 2007-2009 Illumina, Inc.

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

Casava::Intensities - Utilities for CASAVA for the intensities directory

Library for the extraction of various information from the intensities directory.

=cut

=head1 AUTHOR

Mauricio Varea

=head1 SYNOPSIS

  use Casava::Intensities;
  ...

=head1 DESCRIPTION

=head2 Overview

Provides a number of intensities-related services for CASAVA;

=head2 Dependencies

  strict, warnings, Carp, Exporter, XML::Simple, File::Spec

=head2 Methods

=cut

package Casava::Intensities;

use strict;
use warnings "all";
use Carp;
use Exporter 'import';
use File::Spec;
use XML::Simple;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install

use Casava::Common::Log;

our @EXPORT_OK = qw(&new $tilesTxt $configXml);
our $configXml = "config.xml";

sub new(;$)
{
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self  = {};

    $self->{PATH}  = undef;

    bless ($self, $class);
    return $self;
}

sub path($;$)
{
    my $self = shift;
    if (@_) { $self->{PATH} = shift }
    return $self->{PATH};
}

sub getConfig
{
    my $self = shift;
    errorExit "ERROR: undefined Configuration file\n   "  unless defined $configXml;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    errorExit "ERROR: $configFile: file does not exist\n   " unless -e $configFile;
    my $configRef = XMLin($configFile, SuppressEmpty => 1) 
                    or croak "ERROR: couldn't load $configFile: $!\n   ";
    errorExit "ERROR: 'Run' element missing from $configFile\n   " unless exists $configRef->{Run};
    my $run = $configRef->{Run};
}


=pod

=item getRunInfo($)

Extract the "Reads" information from the RunInfo.xml if available.

Note that the structure and content of the RunInfo has changed over the RTA versions.
This version returns returns some sort of "union" of these version. Each 'Read' has:
- 1-based Index and Number
- FirstCycle
- LastCycle
- NumCycles
- IsIndexedRead

  use Casava::Intensities;
  my $baseCalls = Casava::BaseCalls->new();
  $baseCalls->path('/run-folder/Data/Intensities/BaseCalls');
  my $runInfo = $baseCalls->getRunInfo();
  foreach $read (@{$runInfo->{Reads}->{Read}}) {print "$read->{Number}: $read->{FirstCycle} - $read->{LastCycle}: $read->{NumCycles} cycles\n";}

B<Parameters:>

  $self - reference to the BaseCalls instance

B<Returns:>

    the number of reads (integer > 0).

=cut

sub getRunInfo
{
    my $self = shift;
    my $runInfoPath = File::Spec->catfile($self->path, (File::Spec->updir(),'') x 2, 'RunInfo.xml');
    if (not -e $runInfoPath)
    {
        logWarning("Couldn't find run info in $runInfoPath");
        return undef;
    }
    my $runInfoRef = XMLin($runInfoPath, ForceArray => ['Read']) or errorExit("ERROR: couldn't load $runInfoPath: $!");
    errorExit("'Run' element not found in $runInfoPath") unless exists $runInfoRef->{'Run'} and defined $runInfoRef->{'Run'};
    my $runRef = $runInfoRef->{'Run'};
    errorExit("'Read' element not found in $runInfoPath") unless exists $runRef->{'Reads'}->{'Read'} and defined $runRef->{'Reads'}->{'Read'};
    my $readsRef = $runRef->{'Reads'}->{'Read'};
    # Sort the reads either by Index or by FirstCycle
    my $sortKey = exists $readsRef->[0]->{'Number'} ? 'Number' : exists $readsRef->[0]->{'Index'} ? 'Index' : exists $readsRef->[0]->{'FirstCycle'} ? 'FirstCycle' : undef;
    my $lastCycle = 0;
    my $readIndex = 0;
    my @sortedReads = defined $sortKey ? (sort { $a->{$sortKey} > $b->{$sortKey} } @$readsRef) : @$readsRef;
    for my $read (@sortedReads)
    {
        ++$readIndex;
        my $firstCycle = $lastCycle + 1;
        if (exists $read->{'FirstCycle'})
        {
            # In older versions of the run info, the "Index" is an element that indicated that the read contains a barcode.
            # In that case, the attributes "Index", "IsIndexedRead" and "NumCycles" do not exist.
            $read->{'IsIndexedRead'} = (exists $read->{'Index'} ? 'Y' : 'N');
            $read->{'Index'} = $readIndex;
            $read->{'Number'} = $readIndex;
            logWarning("Read FirstCycle = $read->{'FirstCycle'}: expected $firstCycle (Read Index $readIndex)") if $read->{'FirstCycle'} != $firstCycle;
            $read->{'NumCycles'} = $read->{'LastCycle'} - $read->{'FirstCycle'} + 1;
            $lastCycle = $read->{'LastCycle'};
        }
        else
        {
            # In newer versions of the run info, the FirstCycle and LastCycle attributes do not exist
            $read->{'FirstCycle'} = $firstCycle;
            $lastCycle += $read->{'NumCycles'};
            $read->{'LastCycle'} = $lastCycle;
            logWarning("Found Read Index = $read->{'Index'}: expected $readIndex") if exists $read->{'Index'} && $read->{'Index'} != $readIndex;
            logWarning("Found Read Number = $read->{'Number'}: expected $readIndex") if exists $read->{'Number'} && $read->{'Number'} != $readIndex;
            $read->{'Index'} = $readIndex;
            $read->{'Number'} = $readIndex;
        }
    }
    $runRef->{'Reads'}->{'Read'} = \@sortedReads;
    return $runRef;
}

sub getRunParameters
{
    my $self = shift;
    my $run = $self->getConfig;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    errorExit "ERROR: 'RunParameters' element missing from $configFile" unless exists $run->{RunParameters};
    my $runParameters = $run->{RunParameters};
    return $runParameters;
}

sub getSoftware
{
    my $self = shift;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    if (!-e $configFile)
    {
        logWarning("Couldn't determine image processing software version due to missing $configFile");
        return undef;
    }
    my $run = $self->getConfig;
    errorExit "ERROR: 'RunParameters' element missing from $configFile" unless exists $run->{Software};
    my $software = $run->{Software};
    return $software;
}

sub getRtaConfiguration
{
    my $self = shift;
    my $rtaConfigurationPath = File::Spec->catfile($self->path, 'RTAConfiguration.xml');
    logWarning("Missing $rtaConfigurationPath")  unless (-e $rtaConfigurationPath);
    return (-e $rtaConfigurationPath ? $rtaConfigurationPath : undef);
}

sub getCompressedBcl
{
    my $self = shift;
    my $compressedBcl;

    my $rtaConfigurationPath = $self->getRtaConfiguration();
    if (defined $rtaConfigurationPath)
    {
        my $rtaConfigurationRef = XMLin($rtaConfigurationPath)  or errorExit("ERROR: couldn't load $rtaConfigurationPath: $!");
        $compressedBcl = $rtaConfigurationRef->{CompressBCLs};
        logWarning("CompressBCLs' element not found in $rtaConfigurationPath")  unless defined $compressedBcl;
    }

    return $compressedBcl;
}

sub getLocationFileType
{
    my $self = shift;
    my $locationFileType;

    my $rtaConfigurationPath = $self->getRtaConfiguration();
    if (defined $rtaConfigurationPath)
    {
        my $rtaConfigurationRef = XMLin($rtaConfigurationPath) or errorExit("ERROR: couldn't load $rtaConfigurationPath: $!");
        $locationFileType = $rtaConfigurationRef->{LocationFileType};
        logWarning("'LocationFileType' element not found in $rtaConfigurationPath") unless defined $locationFileType;
    }

    return $locationFileType;
}

sub getInstrumentType
{
    my $self = shift;
    my $instrumentType;

    my $rtaConfigurationPath = $self->getRtaConfiguration();
    if (defined $rtaConfigurationPath)
    {
        my $rtaConfigurationRef = XMLin($rtaConfigurationPath) or errorExit("ERROR: couldn't load $rtaConfigurationPath: $!");
        $instrumentType = $rtaConfigurationRef->{InstrumentType};
        if (!defined $instrumentType)
        {
            # an attempt to determine whether it is HiSeq or not in older RTAs
            my $swath = $rtaConfigurationRef->{IsSwathData};
            if (defined $swath)
            {
                $instrumentType = "GA"     if $swath eq 'false';
                $instrumentType = "HiSeq"  if $swath eq 'true';
            }
        }
        logWarning("Couldn't determine control software platform due to missing 'InstrumentType' and 'IsSwathData' in $rtaConfigurationPath")
                  unless defined $instrumentType;
    }

    return $instrumentType;
}

sub getControlSoftware
{
    my $self = shift;
    my $controlSoftware = {};
    
    my $runParametersPath = File::Spec->catfile($self->path, (File::Spec->updir(),'') x 2, 'runParameters.xml');
    if (not -e $runParametersPath)
    {
        logWarning("Couldn't determine control software version due to missing $runParametersPath");
        $controlSoftware->{Name} = 'SCS';
    }
    else
    {
        my $runParametersRef = XMLin($runParametersPath) or errorExit("ERROR: couldn't load $runParametersPath: $!");
        my $applicationName = $runParametersRef->{Setup}->{ApplicationName};
        my $applicationVersion = $runParametersRef->{Setup}->{ApplicationVersion};
        errorExit("'ApplicationName' element not found in $runParametersPath") unless defined $applicationName;
        errorExit("'ApplicationVersion' element not found in $runParametersPath") unless defined $applicationVersion;
        $controlSoftware->{Name} = $applicationName;
        $controlSoftware->{Version} = $applicationVersion;
    }

    $controlSoftware->{Platform} = $self->getInstrumentType();

    return $controlSoftware;
}


sub getSoftwareName
{
    my $self = shift;
    my $software = $self->getSoftware;
    return $software->{Name};
}

sub getSoftwareVersion
{
    my $self = shift;
    my $software = $self->getSoftware;
    return $software->{Version};
}

sub getFormattedVersion
{
    my $self = shift;
    my $versionString = $self->getSoftwareVersion;
    if ($versionString =~ m/(\d+)\.(\d+)\.?(\d*).*/)
    {
    	my @fmtVersion =( sprintf("%d.%02d", $1, $2) );
    	push @fmtVersion, $3  if defined $3;
        return @fmtVersion;
    }
    return undef;
}

sub getRunFolderFromXml($)
{
    my $self = shift;
    my $runParameters = $self->getRunParameters;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    errorExit "ERROR: 'RunFolder' element missing from $configFile" unless exists $runParameters->{RunFolder};
    my $runFolder = $runParameters->{RunFolder};
    return $runFolder;
}

sub getRunParameter($$;$)
{
    my $self = shift;
    my ($parameter, $optional) = @_;
    my $runParameters = $self->getRunParameters;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    unless (exists $runParameters->{$parameter})
    {
        errorExit "ERROR: '$parameter' element missing from $configFile" unless $optional;
        return '';
    }
    return $runParameters->{$parameter};
}

sub instrumentName($)
{
    my $self = shift;
    return $self->getRunParameter('Instrument');
}

sub flowCellId($)
{
    my $self = shift;
    return getFlowCellId(File::Spec->catfile($self->path, $configXml));
}

sub runNumber($)
{
    my $self = shift;
    return $self->getRunParameter('RunFolderId');
}



=pod

=item getReads($)

Extract the number of reads from the 'config.xml' in the base call directory.

  use Casava::BaseCall;
  Casava::Alignment::getReads($baseCallDir);

B<Parameters:>

  $baseCallDir - full path to the base call folder

B<Returns:>

    the number of reads (integer > 0).

=cut

sub getReads($)
{
    my $self = shift;
    croak "ERROR: undefined Configuration file\n   "  unless defined $configXml;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    croak "ERROR: $configFile: file does not exist\n   " unless -e $configFile;
    my $configRef = XMLin($configFile, ForceArray => 1, SuppressEmpty => 1) or croak "ERROR: couldn't load $configFile: $!";
    croak "ERROR: 'Run' element missing from $configFile" unless exists $configRef->{Run};
    my $run = $configRef->{Run}->[0];
    croak "ERROR: 'Run/RunParameters' element missing from $configFile" unless exists $run->{RunParameters};
    my $runParameters = $run->{RunParameters}->[0];
    croak "ERROR: 'Run/RunParameters/Reads' element missing from $configFile" unless exists $runParameters->{Reads};
    return scalar @{$runParameters->{Reads}}
}


1;
__END__

