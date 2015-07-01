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

Casava::BaseCalls - Utilities for CASAVA for the base call directory

Library for the extraction of various information from the base call directory.

=cut

=head1 DESCRIPTION

This script is a rewrite of the original GERALD library originally
implemented for the Pipeline.

=cut

=head1 AUTHOR

Come Raczy

=head1 SYNOPSIS

  use Casava::BaseCalls;
  Casava::BaseCalls::getTiles($baseCallDir);

=head1 DESCRIPTION

=head2 Overview

Provides a number of base call-related services for CASAVA;

=head2 Exports

  getTiles($baseCallDir)

Global variables:

=head2 Dependencies

  strict, warnings, Carp, Exporter, XML::Simple, File::Spec

=head2 Methods

=cut

package Casava::BaseCalls;

use strict;
use warnings "all";
use Carp;
use Exporter 'import';
use File::Basename;
use File::Copy;
use File::Path qw(mkpath);
use File::Spec;
use XML::Simple;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install

use Casava::Common::Log;
use Casava::Common::Utils qw(expandUseBases expandUseBasesString getFlowCellId);

our @EXPORT_OK = qw(&new $tilesTxt $configXml);
our $tilesTxt  = "tiles.txt";
our $configXml = "config.xml";

sub new(;$)
{
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self  = {};

    $self->{PATH}  = undef;
    $self->{TILES} = {};

    bless ($self, $class);
    return $self;
}

sub path($;$)
{
    my $self = shift;
    if (@_) { $self->{PATH} = shift }
    return $self->{PATH};
}

sub lanes($)
{
    my $self = shift;
    return keys %{$self->{TILES}};
}

sub tiles($$;$)
{
    my $self = shift;
    my $lane = shift;
    $lane =~ s/s_(\d)/$1/;
    if (@_) { $self->{TILES}->{"s_$lane"} = shift };
    return undef unless exists $self->{TILES}->{"s_$lane"};
    return $self->{TILES}->{"s_$lane"};
}

sub addTiles($$@)
{
    my $self = shift;
    my $lane = shift;
    $lane =~ s/s_(\d)/$1/;
    my @tiles = @_;
    $self->{TILES}->{"s_$lane"} = []  unless exists $self->{TILES}->{"s_$lane"};
    push @{$self->{TILES}->{"s_$lane"}}, @tiles;
}

sub fillTilesFromDir($)
{
    my $self = shift;

    opendir (DIR,$self->path)
         or croak "ERROR: Cannot open basecall directory '".$self->path()."'";

    my %tilesMap = map{ $_ =~ /(s_\d)/; {$1 => []} }
                   grep{ /s_\d_\d+_\d+_qseq\.txt/ }
                   readdir( DIR );

    foreach my $lane (keys %tilesMap)
    {
        rewinddir(DIR);
        my %tiles = map{ $_ =~ /(\d+)_qseq\.txt/; {$1 => 1} }
                    grep{ /${lane}_\d+_\d+_qseq\.txt/ }
                    readdir( DIR );

        $tilesMap{$lane} = [ sort keys(%tiles) ];
    }

    closedir(DIR);

    $self->{TILES} = \%tilesMap;
}

sub fillTilesFromXml($)
{
    my $self = shift;

    croak "ERROR: undefined Configuration file\n   "  unless defined $configXml;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    croak "ERROR: $configFile: file does not exist\n   " unless -e $configFile;
    my $configRef = XMLin($configFile, SuppressEmpty => 1, ForceArray => ['Lane', 'SelectedTiles', 'Tile', 'TileRange']) 
                    or croak "ERROR: couldn't load $configFile: $!\n   ";
    croak "ERROR: 'Run' element missing from $configFile\n   " unless exists $configRef->{Run};
    my $run = $configRef->{Run};

    # Check if individual tiles have been selected.
    # Needed because only first tile range is recorded in some PL versions.
    croak "ERROR: 'RunParameters' element missing from $configFile"
          unless exists $run->{RunParameters};
    my $runParameters = $run->{RunParameters};

    my %laneSelectedTilesMap;

    if (exists($runParameters->{SelectedTiles}))
    {
        my @selectedTiles = @{$runParameters->{SelectedTiles}};

        foreach my $selectedTile (@selectedTiles)
        {
            next if (!($selectedTile =~ /.*_(\d)_(\d{4})/));

            my $laneNum = $1;
            my $tileNum = $2;
            push(@{$laneSelectedTilesMap{$laneNum}}, $tileNum);
        }
    }

    croak "ERROR: 'Run/TileSelection' element missing from $configFile" unless exists $run->{TileSelection};
    my $tileSelection = $run->{TileSelection};
    croak "ERROR: 'Run/TileSelection/Lane' element missing from $configFile" unless exists $tileSelection->{Lane};
    my %tilesMap;
    foreach my $lane (@{$tileSelection->{Lane}})
    {
        croak "ERROR: lane without 'Index'" unless exists $lane->{Index};
        my $index = $lane->{Index};
        croak "ERROR: missing 'Sample' in lane $index" unless exists $lane->{Sample};
        my $sample = $lane->{Sample};
        my $key = $sample."_${index}";
        croak "ERROR: multiple occurrences of lane Index $index" if exists $tilesMap{$key};
        my @tiles;

        # In some Pipeline versions, only the first tile range was stored,
        # so use the SelectedTiles if present for this lane instead of
        # concatenating ranges.

        if (exists($laneSelectedTilesMap{$index}))
        {
            @tiles = sort(@{$laneSelectedTilesMap{$index}});
        }
        else
        {
            @tiles = @{$lane->{Tile}} if exists $lane->{Tile};
            foreach my $tileRange (@{$lane->{TileRange}})
              {
                  croak "ERROR: missing 'Min' in tile range for lane $key" unless exists $tileRange->{Min};
                  croak "ERROR: missing 'Max' in tile range for lane $key" unless exists $tileRange->{Max};
                  my $min = $tileRange->{Min};
                  my $max = $tileRange->{Max};
                  push @tiles, ($min..$max);
              }
            @tiles = map {sprintf "%04d", $_} @tiles;
        }

        $tilesMap{$key} = \@tiles;
    }
    $self->{TILES} = \%tilesMap;
}

sub filterTiles($$)
{
    my $self = shift;
    my ($tilesFilter) = @_;
    my $tilesRegex = $tilesFilter;
    $tilesRegex =~ s/,/|/g;
    my @allTiles;
    foreach my $lane (keys %{$self->{TILES}})
    {
        my $tiles = $self->{TILES}->{$lane};
        push @allTiles, map("${lane}_$_", @$tiles);
    }
    my @selectedTiles = grep(/$tilesRegex/, @allTiles);
    errorExit("ERROR: tiles filter '$tilesFilter' resulted into an empty tile selection") unless @selectedTiles;
    my %tilesMap;
    foreach my $tile (@selectedTiles)
    {
        next unless $tile =~ /^(.+)_(\d{4})$/;
        $tilesMap{$1} = [] unless exists $tilesMap{$1};
        push @{$tilesMap{$1}}, $2;
    }
    $self->{TILES} = \%tilesMap;
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
    my $run = $self->getConfig;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    errorExit "ERROR: 'Software' element missing from $configFile" unless exists $run->{Software};
    my $software = $run->{Software};
    return $software;
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
    # First try to get the RunParameters/RunFlowcellId element
    my $runParameters = $self->getRunParameters;
    my $runFlowcellId = $runParameters->{'RunFlowcellId'} if exists $runParameters->{'RunFlowcellId'};
    return $runFlowcellId if defined $runFlowcellId and 0 < length($runFlowcellId);
    # Otherwise fallback to the old method
    return getFlowCellId(File::Spec->catfile($self->path, $configXml));
}

sub runNumber($)
{
    my $self = shift;
    return $self->getRunParameter('RunFolderId');
}

sub writeTiles($@)
{
    my $self = shift;
    my @lanes = @_;
    croak "ERROR: undefined Tiles file\n   "  unless defined $tilesTxt;
    my $tilesFile = File::Spec->catfile($self->path, $tilesTxt);

    open( my $tilesHandle, ">$tilesFile" )  or croak "ERROR: Failed to open file $tilesFile: $!\n    ";

    my @sortedLanes = sort @lanes;
    foreach my $lane (@sortedLanes)
    {
        next unless defined $self->tiles($lane);
        my @tiles = @{$self->tiles($lane)};
        my @sortedTiles = sort @tiles;
        foreach my $tile (@sortedTiles)
        {
            print $tilesHandle "s_${lane}_${tile}";
            if ($lane == $sortedLanes[$#lanes] && $tile == $sortedTiles[$#tiles]) {
                print $tilesHandle "\n";
            } else {
                print $tilesHandle " ";
            }
        }
    }
    close( $tilesHandle );
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


=pod

=item guessUseBases()

Make an educated guess about what the use bases should be, based on the 'config.xml' $self->path.

  use Casava::BaseCall;
  Casava::Alignment::guessUseBases();

B<Parameters:>

  $baseCallDir - full path to the base call folder

B<Returns:>

    guessed Use Bases parameter.

=cut

sub guessUseBases($)
{
    my $self = shift;
    my $useBases = undef;
    croak "ERROR: undefined Configuration file\n   "  unless defined $configXml;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    croak "ERROR: $configFile: file does not exist\n   " unless -e $configFile;
    my $configRef = XMLin($configFile, SuppressEmpty => 1, ForceArray => 1, ForceContent => 1) 
       or croak "ERROR: couldn't load $configFile: $!   ";
    croak "ERROR: 'Run' element missing from $configFile\n   " unless exists $configRef->{Run};
    my $run = $configRef->{Run}->[0];
    croak "ERROR: 'Run/RunParameters' element missing from $configFile\n   " unless exists $run->{RunParameters};
    my $runParameters = $run->{RunParameters}->[0];
    croak "ERROR: 'Run/RunParameters/Reads' element missing from $configFile\n   " unless exists $runParameters->{Reads};
    my @reads = @{$runParameters->{Reads}};
    my $barcode = $runParameters->{Barcode}->[0];

    my @cycle = (defined $barcode) ? @{$barcode->{'Cycle'}} : ({content=>0});
    my $indexMask;
    my $indexPos = -1; #assume no barcode read initially
    for (my $i = 0; $i < scalar(@reads); $i++)
    {
        if( ($reads[$i]->{'FirstCycle'}->[0]->{'content'} == $cycle[0]->{'content'}) &&
            ($reads[$i]->{'LastCycle'}->[0]->{'content'}  == $cycle[$#cycle]->{'content'}) )
        {
            $indexPos = $reads[$i]->{'Index'};
            $indexMask = join ('', 
                               map { ( exists $_->{'Use'} 
                                    && $_->{'Use'} eq "false" ) ? 'n' : 'I' } @cycle 
                         );
            my $len = length($indexMask);
            if ($len < (scalar @cycle))
            {
                while ($len > 1)
                {
                    my $ii = 'I' x $len;
                    $indexMask =~ s/$ii/I$len/g; 
                    $len--;
                };
            } else {
                $indexMask = "I".$#cycle."n";
            }
            last;
        }
    }
    my @readLength = map{ $_->{'LastCycle'}->[0]->{'content'} - $_->{'FirstCycle'}->[0]->{'content'} + 1 } @reads;
    $useBases = join(',', 
                     map { ($_->{'Index'} != $indexPos) ? ("Y".$readLength[$_->{'Index'} - 1]) : $indexMask } @reads
                );
    return $useBases;
}

=pod

=item allCycles($)

Retrive the list of all cycles cycles from 'config.xml' in the base call directory

  use Casava::BaseCall;
  Casava::Alignment::barcodeCycles($mask);

B<Returns:>

    sorted array of cycle numbers

=cut

sub allCycles($)
{
    my ($self) = @_;
    croak "ERROR: undefined Configuration file\n   "  unless defined $configXml;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    croak "ERROR: $configFile: file does not exist\n   " unless -e $configFile;
    my $configRef = XMLin($configFile, SuppressEmpty => 1, ForceArray => 1, ForceContent => 1) 
       or croak "ERROR: couldn't load $configFile: $!   ";
    croak "ERROR: 'Run' element missing from $configFile\n   " unless exists $configRef->{Run};
    my $run = $configRef->{Run}->[0];
    croak "ERROR: 'Run/RunParameters' element missing from $configFile\n   " unless exists $run->{RunParameters};
    my $runParameters = $run->{RunParameters}->[0];
    croak "ERROR: 'Run/RunParameters/Reads' element missing from $configFile\n   " unless exists $runParameters->{Reads};
    my @reads = @{$runParameters->{Reads}};

    my @allCycles=();
    for (my $i = 0; $i < scalar(@reads); $i++)
    {
        push @allCycles, ($reads[$i]->{'FirstCycle'}->[0]->{'content'} .. $reads[$i]->{'LastCycle'}->[0]->{'content'});
    }
    return sort {$a<=>$b} @allCycles;
}

=pod

=item barcodeCycles($$)

Identify the barcode cycles from the use base mask. If the mask contains fillers,
retrive the list of barcode cycles from 'config.xml' in the base call directory
given the use base mask.

  use Casava::BaseCall;
  Casava::Alignment::barcodeCycles($mask);

B<Parameters:>

  $mask - use base mask with I corresponding to cycle numbers that will be returned

B<Returns:>

    reference to cycle numbers array or undef.

=cut

sub barcodeCycles($$)
{
    my ($self, $mask) = @_;
    my $barcodeCycles = [];
    # If all index cycles indicated in the mask are before any filler, use the mask
    if ($mask and $mask =~ /^([0-9yYnNiI,]*)[^iI]*$/ and $mask !~ /[iI]\*/)
    {
        my @tmp = split('\*', $mask);
        my $currentCycle = 0;
        for my $read (split(',', $tmp[0]))
        {
            foreach my $cycleSelector (split('', expandUseBasesString($read)))
            {
                ++$currentCycle;
                push @$barcodeCycles, $currentCycle if $cycleSelector eq 'I' or $cycleSelector eq 'i';
            }
        }
    }
    else
    {
        # The mask doesn't provide enough information. Fall back to the config file.
        croak "ERROR: undefined Configuration file\n   "  unless defined $configXml;
        my $configFile = File::Spec->catfile($self->path, $configXml);
        croak "ERROR: $configFile: file does not exist\n   " unless -e $configFile;
        my $configRef = XMLin($configFile, SuppressEmpty => 1, ForceArray => 1, ForceContent => 1) 
        or croak "ERROR: couldn't load $configFile: $!   ";
        croak "ERROR: 'Run' element missing from $configFile\n   " unless exists $configRef->{Run};
        my $run = $configRef->{Run}->[0];
        croak "ERROR: 'Run/RunParameters' element missing from $configFile\n   " unless exists $run->{RunParameters};
        my $runParameters = $run->{RunParameters}->[0];
        croak "ERROR: 'Run/RunParameters/Reads' element missing from $configFile\n   " unless exists $runParameters->{Reads};
        my @reads = @{$runParameters->{Reads}};
        my $barcode = $runParameters->{Barcode}->[0];
        
        my @useBases = split(',', $mask);
        push @useBases, 'n*' while scalar(@useBases) < scalar(@reads);
        for (my $i = 0; $i < scalar(@reads); $i++)
        {
            my $length = 1 + $reads[$i]->{'LastCycle'}->[0]->{'content'} - $reads[$i]->{'FirstCycle'}->[0]->{'content'};
            my $readMask = shift @useBases;
            my @expandedUseBases = split //, expandUseBases($readMask, $length);
            my @thisReadBarcodeCycles = grep 
            {'i' eq $expandedUseBases[$_ - $reads[$i]->{'FirstCycle'}->[0]->{'content'}]} 
            ($reads[$i]->{'FirstCycle'}->[0]->{'content'} .. $reads[$i]->{'LastCycle'}->[0]->{'content'});
            push @$barcodeCycles, @thisReadBarcodeCycles;
        }
    }
    return $barcodeCycles;
}

=pod

=item readCycles($)

Retrive the list of regular read cycles from 'config.xml' in the base call directory.

  use Casava::BaseCall;
  Casava::Alignment::readCycles($mask);

B<Parameters:>

  $mask - use base mask with 'Y' or 'y' corresponding to cycle numbers that will be returned
  $demuxCycles - if set, the cycle numbers after masking are returned. Ohterwise original cycle
                 numbers are returned.

B<Returns:>

    reference to hash of cycle numbers array refs or undef. Only cycle arrays that stayed
    non-empty after masking are returned.

=cut

sub readCycles($$$)
{
    my ($self, $mask, $demuxCycles) = @_;
    
    my $nonBarcodeCycles = {};
    croak "ERROR: undefined Configuration file\n   "  unless defined $configXml;
    my $configFile = File::Spec->catfile($self->path, $configXml);
    croak "ERROR: $configFile: file does not exist\n   " unless -e $configFile;
    my $configRef = XMLin($configFile, SuppressEmpty => 1, ForceArray => 1, ForceContent => 1) 
       or croak "ERROR: couldn't load $configFile: $!   ";
    croak "ERROR: 'Run' element missing from $configFile\n   " unless exists $configRef->{Run};
    my $run = $configRef->{Run}->[0];
    croak "ERROR: 'Run/RunParameters' element missing from $configFile\n   " unless exists $run->{RunParameters};
    my $runParameters = $run->{RunParameters}->[0];
    croak "ERROR: 'Run/RunParameters/Reads' element missing from $configFile\n   " unless exists $runParameters->{Reads};
    my @reads;
    # Ensure that the reads are in the right order (assume contiguous index, starting at 1)
    foreach my $read (@{$runParameters->{Reads}})
    {
        $reads[$read->{'Index'} - 1] = $read;
    }
    
    my @useBases = split(',', $mask);
    push @useBases, 'n*' while scalar(@useBases) < scalar(@reads);
    my $cyclesMasked = 0;
    #for (my $i = 0; $i < scalar(@reads); $i++)
    my $currentReadIndex = 0;
    my $currentCycle = 1;
    for my $readMask (@useBases)
    {
        ++$currentReadIndex;
        my $read = shift @reads;
        my $firstCycle = undef;
        $firstCycle = $read->{'FirstCycle'}->[0]->{'content'} if $read;
        logWarning("Read $currentReadIndex: current cycle in mask is $currentCycle: first cycle in read is $firstCycle")
            if defined $firstCycle and $firstCycle != $currentCycle;
        my $length = undef;
        my $expectedLength = undef;
        $expectedLength = 1 + $read->{'LastCycle'}->[0]->{'content'} - $read->{'FirstCycle'}->[0]->{'content'} if $read;
        if ($readMask !~ /\*/)
        {
            $readMask = expandUseBasesString($readMask);
            $length = length($readMask);
            logWarning("use-base-mask length is $length for read $currentReadIndex:" .
                       " expected length from config file is $expectedLength") if $read and $length != $expectedLength;
        }
        else
        {
            $length = $expectedLength;
        }
        logInfo("Read $currentReadIndex: length = $length: mask = $readMask", 0)  unless $demuxCycles;
        my @expandedUseBases = split //, expandUseBases($readMask, $length);
        #my @thisReadCycles = grep 
        #    {'y' eq $expandedUseBases[$_ - $reads[$i]->{'FirstCycle'}->[0]->{'content'}]} 
        #    ($reads[$i]->{'FirstCycle'}->[0]->{'content'} .. $reads[$i]->{'LastCycle'}->[0]->{'content'});
        #my $indexPos = $reads[$i]->{'Index'};
        my @thisReadCycles = grep {'y' eq $expandedUseBases[$_ - $currentCycle]} ($currentCycle .. ($currentCycle + $length - 1));
        my $indexPos = $currentReadIndex;
        if (@thisReadCycles)
        {
            croak "ERROR: Masking out cycles in the middle of the read is not allowed: read $indexPos mask '$readMask' in '$mask' is bad\n   " 
                unless scalar(@thisReadCycles) == 1 + ($thisReadCycles[-1] - $thisReadCycles[0]);
            my @shiftedThisReadCycles = map {$_ - $cyclesMasked} @thisReadCycles;
            $nonBarcodeCycles->{$indexPos} = 
                $demuxCycles ? \@shiftedThisReadCycles : \@thisReadCycles if (scalar(@thisReadCycles));
        }
        $currentCycle += $length;
        $cyclesMasked += $length - scalar(@thisReadCycles);
    }
    
    
#    my $barcode = $runParameters->{Barcode}->[0];
#
#    # non-multiplexed data does not have Barcode element
#    my @cycle = ({'content'=>0});
#    if (defined $barcode)
#    {
#        @cycle = @{$barcode->{'Cycle'}};
#    }
#    my $indexMask;
#    my $indexPos;
#    for (my $i = 0; $i < scalar(@reads); $i++)
#    {
#        if( ($reads[$i]->{'FirstCycle'}->[0]->{'content'} != $cycle[0]->{'content'}) &&
#            ($reads[$i]->{'LastCycle'}->[0]->{'content'}  != $cycle[$#cycle]->{'content'}) )
#        {
#            $indexPos = $reads[$i]->{'Index'};
#            my @cyclesArray = ($reads[$i]->{'FirstCycle'}->[0]->{'content'} .. $reads[$i]->{'LastCycle'}->[0]->{'content'});
#            $nonBarcodeCycles->{$indexPos} = \@cyclesArray;
#        }
#    }
    
    return keys(%$nonBarcodeCycles) ? $nonBarcodeCycles : undef;
}


sub initialize($;$)
{
    my ($self,$scall) = @_;

    my @emptyDirs = ($self->path);

    # these are filled up at 'make' stage
    foreach my $dir (@emptyDirs)
    {
        croak "ERROR: file '$dir' is supposed to be a directory\n   "
              if -f $dir;
        unless (-e $dir)
        {
            print STDERR "Creating directory '$dir'\n";
            mkpath $dir;
            croak "ERROR: failed to create '$dir'\n" unless -d $dir;
        }
    }

#    # these are copied across at 'perl' run-time
#    if (defined $subdirs)
#    {
#        foreach my $dir (values %$subdirs)
#        {
#            my $target = File::Spec->catdir($self->path,basename($dir));
#            croak "ERROR: '$target' already exists.\n   "  if -e $target;
#            my @files = glob(File::Spec->catfile($dir,'*'));
#            print "Copying '$dir' into '$target'...\n";
#            mkdir $target;
#            copy($_,$target)  foreach (@files);
#        }
#    }
}

1;
__END__

