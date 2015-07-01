#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2010-2011 Illumina, Inc.

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

copyConfig.pl

Copy the BaseCalls "config.xml" while applying some modifications.

=cut

=head1 SYNOPSIS

=over 12

=item B<copyConfig.pl>

[S<B<--input-file> I<input_config>>]
[S<B<--output-file> I<output_config>>]
[S<B<--tiles> I<tile_selection>>]

=item B<copyConfig.pl>

B<--help> or B<--man>

=back

=head1 OPTIONS AND ARGUMENTS

=over 8

=item B<--input-file> I<input_config>

Full path to a valid BaseCalls XML configuration file.

=item B<--output-file> I<output_config>

Full path to the modified file. Prints to the standard output if unspecified.

=item B<--tiles> I<regex>[[I<,regex>]...]

Comma-separated list of regular expressions to select only a subset of the tiles available in the flow-cell.

 - to select all the tiles ending with "5" in all lanes: --tiles [0-9][0-9][0-9]5
 - to select tile 2 in lane 1 and all the tiles in the other lanes: --tiles s_1_0002,s_[2-8]

=item B<--help> Print a brief help message and exit.

Z<>

=item B<--man> View this help formatted in 'man' style.

Z<>

=back

=head1 DESCRIPTION

B<copyConfig.pl> reads a given Base Call XML configuration file, applies the
required modification and saved the resulting configuration file.

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

Come Raczy

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
use XML::Simple;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Alignment;
use Casava::BaseCalls;
use Casava::Common::Utils qw(reallyRealPath);
use Casava::Common::Log qw(logWarning logInfo errorExit initLog);

sub filterTiles($$);

initLog( undef, 3, 0, 1, 1);

my @CLI = @ARGV;

my $man = 0;
my $help = 0;
my $inputFile;
my $outputFile;
my $tilesFilter;

my $result = GetOptions('help|?' => \$help,
                        'input-file=s' => \$inputFile,
                        'output-file=s' => \$outputFile,
                        'tiles=s' => \$tilesFilter,
                        man => \$man) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2,  -input => $1) if ($man and $0 =~ /(.*)/);

errorExit "ERROR: --input-file is unspecified\n" unless defined $inputFile and length($inputFile);

my $configRef = XMLin($inputFile,
                      SuppressEmpty => 1,
                      ForceArray => ['Lane', 'SelectedTiles', 'Tile', 'TileRange', 'Sample',
                                     'BaseCallFolder', 'ImageAnalysisFolder', 'InstrumentFolder',
                                     'IparFlag', 'Compression', 'CompressionSuffix',
                                     'ImageHeight', 'ImageWidth', 'NumMatrixCycles', 'QhgCompression',
                                     'SeqCompression', 'Sig2Compression',
                                     'PhasingRate', 'PrephasingRate', 'RunFlowcellId',
                                     'AutoFlag', 'AutoLane', 'Cycle', 'CycleOffset',
                                     'AutoCycleFlag', 'BasecallFlag', 'Deblocked', 'DebugFlag',
                                     'ChastityThreshold', 'PureBases', 'SmtFilter', 'SmtRelation', 'SmtThreshold',
                                     'FirstRunOnlyFlag', 'Instrument', 'IterativeMatrixFlag', 'MakeFlag',
                                     'MaxCycle', 'MinCycle', 'QTableVersion', 'RunFolderDate', 'RunFolderId',
                                     'FirstCycle', 'LastCycle', 'RunFolder', 'Read'])
or errorExit "ERROR: couldn't load config file $inputFile: $!\n   ";

if(defined $tilesFilter and length($tilesFilter))
{
    filterTiles($configRef, $tilesFilter);
}

$outputFile = '&STDOUT' unless defined $outputFile and length($outputFile);

open my $output, ">$outputFile" or errorExit("ERROR: couldn't open output file $outputFile: $!");

my $xml = XMLout($configRef,
                 RootName => 'BaseCallAnalysis',
                 XMLDecl => '<?xml version="1.0" encoding="utf-8"?>',
                 OutputFile => $output) or errorExit "ERROR: couldn't write config file $outputFile: $!\n   ";

close $output;

1;

sub filterTiles($$)
{
    my ($configRef, $tilesFilter) = @_;
    my $tilesRegex = $tilesFilter;
    $tilesRegex =~ s/,/|/g;
    if (not exists($configRef->{'Run'}->{'TileSelection'}->{'Lane'}))
    {
        logWarning("WARNING: config file doesn't have any 'TileSelection/Lane' element");
        return;
    }
    my @laneList = @{$configRef->{'Run'}->{'TileSelection'}->{'Lane'}};
    my @allTiles;
    foreach my $lane (@laneList)
    {
        next unless exists $lane->{Sample} and exists $lane->{Index};
        my $sample = $lane->{Sample}->[0];
        my $index = $lane->{Index};
        my @tiles;
        push @tiles, @{$lane->{Tile}} if exists $lane->{Tile};
        if (exists $lane->{TileRange})
        {
            foreach my $tileRange (@{$lane->{TileRange}})
            {
                next unless exists $tileRange->{Min} and exists $tileRange->{Max};
                push @tiles, $tileRange->{Min}..$tileRange->{Max};
            }
        }
        push @allTiles, map(sprintf("${sample}_${index}_%04d", $_), @tiles);
    }
    my @selectedTiles = grep(/$tilesRegex/, @allTiles);
    errorExit("ERROR: tiles filter '$tilesFilter' resulted into an empty tile selection") unless @selectedTiles;
    my %tileMap;
    foreach my $tile (@selectedTiles)
    {
        next unless $tile =~ /^(.+_\d)_(\d{4})$/;
        $tileMap{$1} = [] unless exists $tileMap{$1};
        push @{$tileMap{$1}}, sprintf("%d", $2);
    }
    my @tileSelection;
    foreach my $lane (sort keys %tileMap)
    {
        next unless $lane =~ /^(.*)_(\d)$/;
        my $sample = $1;
        my $index = $2;
        push @tileSelection, {Sample => [$sample], Index => $index, Tile => $tileMap{$lane}};
    }
    $configRef->{'Run'}->{'TileSelection'}->{'Lane'} = \@tileSelection;
}

__END__

