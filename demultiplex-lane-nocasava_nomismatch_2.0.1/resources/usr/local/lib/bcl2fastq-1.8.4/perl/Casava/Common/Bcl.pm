package Casava::Common::Bcl;

use File::Spec;
use Data::Dumper;
use Carp;
use Casava::Common::Eamss qw(maskQvalsByEamss);

BEGIN {
    use Exporter();
    use vars qw (@ISA @EXPORT @EXPORT_OK);
    @ISA    = qw(Exporter);
    @EXPORT = qw(%bclDef $BCL_FIELD_SEPARATOR);
    @EXPORT_OK = qw(&new &close &getRead &getFilter &getClusterCount);
}

# PROJECT: Pipeline
# MODULE:  Bcl.pm
# AUTHOR:  Come Raczy
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

# _bcl file API functions

=pod

=head1 NAME

Casava::Common::Bcl.pm - _bcl file API functions

=head2 SYNOPSIS

use Casava::Common::Bcl.pm qw();  

=head2 AUTHORSHIP

Copyright (c) 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

=head1 DESCRIPTION

Bcl Read file (_bcl file API functions)

=head2 Overview

Definition of the _qseq data format

Currently tab separated text format.

    * Machine name: (hopefully) unique identifier of the sequencer. [:alnum:]
    * Run number: (hopefully) unique number to identify the run on the sequencer
    * Lane number: (strictly) positive integer (currently 1-8)
    * Tile number: (strictly) positive integer
    * X: x coordinate of the spot. Integer (can be negative)
    * Y: y coordinate of the spot. Integer (can be negativE)
    * Index: (strictly) positive integer. No indexing should have a value of 1.
    * Read Number: 1 for single reads; 1 or 2 for paired ends
    * Sequence: the sequence (ACGT.)
    * Quality: the calibrated quality string
    * Filter: Did the read passed filter (0 - No, 1 - Yes)

1 file for 1 "virtual lane" (concatenation of the lane number on 1 digit, with the index number on 4 digits), names s_lllll_seq.txt, where "lllll" is the id of the virtual lane.

Note: to avoid problems with concurrent access to the files, intermediary "tile-based _qseq files" are created. When they are all completed, they are dispatched by index number and concatenated to get the final _qseq file.
Messaging Functions 

=head2 Exports

    Casava::Common::Bcl::new($$$$$$$$$$$$$);
    Casava::Common::Bcl::closeQseq($);
    Casava::Common::Bcl::getRead($);

Global variables:

=head2 Depends on
    
    Casava::Common::QseqRead
    warnings, strict, POSIX, constant, IO::File

=cut

use warnings;
use strict;
use POSIX;
use IO::File;
use constant DEBUG => 0;    # set to 1 to get debug info
use Casava::Common::QseqRead;

sub openBcl($$$$$$$$$$$$);
sub close($);
sub getRead($);

sub openPositions($$$$$);
sub getClusterCount($$);

=pod

=head1 FUNCTIONS

=head2 General Functions

=over 4

=cut

=pod

=item new($$$$$$$$$$$$$)

The procedure opens bcl file and returns ref to a structure destribing 
    bcl file. 

B<Parameters:>

    $machineName          - name of the instrument
    $runNumber            - run number on that instrument
    $laneNumber           - lane number (1-8)
    $tileNumber           - tile number (> 0)
    $inputDirectory       - the BaseCalls directory
    $cyclesRef            - arrayref to the list of cycles
    $positionsDirectory   - directory containing the positions file
    $positionsFileFormat  - "txt", "locs" or "clocs"
    $filtersDirectory      - directory containing the filters file

B<Returns:> 

    HASH MAP Ref to Bcl structure:
    qseqRead
    bclFiles
    positionsFile
    filtersFile

=cut

sub new($$$$$$$$$$$$$)
{
    my $class = shift;
    my ($machineName, $runNumber, $laneNumber, $tileNumber, $inputDirectory,
        $cyclesRef, $barcodeCyclesRef, $read, $positionsDirectory,
        $positionsFileFormat, $filtersDirectory, $withEamss) = @_;
    croak "Input directory $inputDirectory does not exist" unless -d $inputDirectory;
    croak "Filters directory $filtersDirectory does not exist" unless -d $filtersDirectory;
    croak "Positions directory $positionsDirectory does not exist" unless -d $positionsDirectory;
    croak "Empty list of cycles" unless @$cyclesRef;
    croak "Invalid read $read" unless 0 < $read;
    croak "Support for barcodes is not implemented yet" if @$barcodeCyclesRef;
    croak "Unknown positions file format $positionsFileFormat" unless grep(/^$positionsFileFormat$/, 'txt', 'locs');

    my $laneDirectory = File::Spec->catfile($inputDirectory, sprintf("L%03u", $laneNumber));
    my $bclName = sprintf("s_%u_%u.bcl", $laneNumber, $tileNumber);
    my @bclFiles;
    # TODO: implement support for barcodes
    # my @barcodeFiles;
    my $clusterCount;
    foreach my $cycle (@$cyclesRef)
    {
        my $bclFilePath = File::Spec->catfile($laneDirectory, sprintf("C%u.1", $cycle), $bclName);
        my $handle = new IO::File $bclFilePath, "r";
        croak "Couldn't open BCL file $bclFilePath: $!" unless $handle;
        binmode($handle);
        my $tmpCount = getClusterCount($handle, $bclFilePath);
        $clusterCount = $tmpCount unless defined $clusterCount;
        croak "Unexpected cluster count in BCL file $bclFilePath: expected $clusterCount: got $tmpCount" unless $tmpCount==$clusterCount;
        push @bclFiles, $handle;
    }
    my $positionsFile = openPositions($positionsDirectory, $positionsFileFormat, $laneNumber, $tileNumber, $clusterCount);
    my $filtersFilePath = File::Spec->catfile($filtersDirectory, sprintf("s_%u_%04u.filter", $laneNumber, $tileNumber));
    my $filtersFile = new IO::File $filtersFilePath;
    croak "Couldn't open filters file $filtersFilePath: $!" unless $filtersFile;
    binmode($filtersFile);
    my $tmpCount = getClusterCount($filtersFile, $filtersFilePath);
    croak "Unexpected cluster count in filter file $filtersFilePath: expected $clusterCount: got $tmpCount" unless $tmpCount==$clusterCount;
    my $qseqRead = Casava::Common::QseqRead->new([$machineName, $runNumber, $laneNumber, $tileNumber, 0, 0, 0, $read, '', '', 0]);
    my $self = {qseqRead=>$qseqRead, bclFiles=>\@bclFiles, positionsFile=>$positionsFile, 
                filtersFile=>$filtersFile, clusterCount=>$clusterCount, clustersRead=>0, withEamss=>$withEamss};
    bless $self, $class;
    return $self;
}    # sub openBcl

=pod

=item close($)

    The procedure closes bcl and associated file

B<Parameters:>

    $self             - HASH MAP Ref to Bcl structure:

B<Returns:> 

    status (0 or 1)

=cut

sub close($) {
    my ($self) = @_;
    foreach my $bclFile (@{$self->{bclFiles}})
    {
        CORE::close($bclFile);
    }
    CORE::close($self->{positionsFile}->{handle});
    CORE::close($self->{filtersFile});
    return 1;
}    # sub close

=pod

=item getQseqRead($)

The procedure reads one read (ARRAY Ref) from bcl file. If end of file
then returns undef.
To iterate

while ( defined( my $readTmpRef = getQseqRead( %{$bclFileRef} ) ) ) {
    my $tile = $readTmpRef->tile();
}

B<Parameters:>

    $self             - HASH MAP Ref to Bcl structure see openFile

B<Returns:> 

    QseqRead object Ref (undef if EOF)

=cut

our @baseLiterals = ('A', 'C', 'G', 'T');

sub getRead($) {
    my ($self) = @_;
    my @bases;
    my @qualities;
    if ($self->{clusterCount} <= $self->{clustersRead})
    {
        return undef;
    }
    foreach my $bclFile(@{$self->{bclFiles}})
    {
        my $c;
        read($bclFile, $c, 1) or croak "Failed to read bcl file $bclFile: $!";
        my $v = unpack("C", $c);
        if ($v)
        {
            push @bases, $baseLiterals[$v & 3];
            push @qualities, chr(64 + ($v >> 2));
        }
        else
        {
            push @bases, 'N';
            push @qualities, chr(66);
        }
    }
    # TODO: apply the EAMSS filter if needed
    if ($self->{withEamss})
    {
        maskQvalsByEamss(\@qualities, \@bases);
    }
    my $coordinates = $self->{positionsFile}->{getCoordinates}();
    my $filter = $self->getFilter();
    my $qseqRead = $self->{qseqRead};
    $qseqRead->X(sprintf("%.0f", 1000.0 + 10.0*$coordinates->[0]));
    $qseqRead->Y(sprintf("%.0f", 1000.0 + 10.0*$coordinates->[1]));
    # TODO: implement support for barcodes
    #$qseqRead->index = 0;
    $qseqRead->seq(join("", @bases));
    $qseqRead->qualityString(join("", @qualities));
    $qseqRead->filter($filter);
    $self->{clustersRead}++;
    return $qseqRead;
}    # sub getQseqRead

=pod

=item writeQseqRead($\%)

The procedure writes one qseq read (ARRAY Ref) to qseq file. 

B<Parameters:>

    $qseqReadRef         - ARRAY MAP Ref to qseq read structure see getQseqRead
    $qseqRef             - QseqRead object Ref

B<Returns:> 

    Nothing

=cut

sub getClusterCount($$)
{
    my ($handle, $filePath) = @_;
    my $buffer;
    my $r = read($handle, $buffer, 4);
    croak "Couldn't read cluster count from $filePath: $!" unless $r and 4==$r;
    my $clusterCount = unpack("L", $buffer);
    return $clusterCount;
}

sub openTxtPositions($$$$);
sub openLocsPositions($$$$);

sub openPositions($$$$$)
{
    my ($positionsDirectory, $positionsFileFormat, $laneNumber, $tileNumber, $clusterCount) = @_;
    if ('txt' eq $positionsFileFormat)
    {
        return openTxtPositions($positionsDirectory, $laneNumber, $tileNumber, $clusterCount)
    }
    elsif ('locs' eq $positionsFileFormat)
    {
        return openLocsPositions($positionsDirectory, $laneNumber, $tileNumber, $clusterCount)
    }
    else
    {
        croak "Unsupported position file format: $positionsFileFormat";
    }
}

sub openTxtPositions($$$$)
{
    my ($positionsDirectory, $laneNumber, $tileNumber, $clusterCount) = @_;
    my $positionsFilePath = File::Spec->catfile($positionsDirectory, sprintf("s_%u_%04u_pos.txt", $laneNumber, $tileNumber));
    my $positionsFile = new IO::File $positionsFilePath, "r";
    croak "Couldn't open positions file $positionsFilePath: $!" unless $positionsFile;
    my $getCoordinates = sub
    {
        croak "Unexpected error in positions file $positionsFilePath" unless $positionsFile;
        my $coordinatesLine;
        $coordinatesLine = <$positionsFile>;
        croak "Failed to read coordinates from positions file $positionsFilePath: $!" unless $coordinatesLine;
        $coordinatesLine =~ /^\s*(-?\d+([.]\d*)?)\s+(-?\d+([.]\d*)?)(\s|\n|\r)*$/ 
        or croak "Failed to read coordinates from positions file $positionsFilePath: $coordinatesLine";
        return [$1, $3];
    };
    return {getCoordinates=>$getCoordinates, handle=>$positionsFile};
}

sub openLocsPositions($$$$)
{
    my ($positionsDirectory, $laneNumber, $tileNumber, $clusterCount) = @_;
    my $positionsFilePath = File::Spec->catfile($positionsDirectory, sprintf("s_%u_%04u.locs", $laneNumber, $tileNumber));
    my $positionsFile = new IO::File $positionsFilePath, "r";
    croak "Couldn't open positions file $positionsFilePath: $!" unless $positionsFile;
    binmode($positionsFile);
    my $tmpCount = getClusterCount($positionsFile, $positionsFilePath);
    croak "Unexpected cluster count in positions file $positionsFilePath: expected $clusterCount: got $tmpCount" unless $tmpCount==$clusterCount;
    my $getCoordinates = sub
    {
        my $c;
        my $r = read($positionsFile, $c, 4);
        croak "Couldn't read X coordinate from positions file $positionsFilePath: $!" unless $r and 4 == $r;
        my $x = unpack("l", $c);
        $r = read($positionsFile, $c, 4);
        croak "Couldn't read Y coordinate from positions file $positionsFilePath: $!" unless $r and 4 == $r;
        my $y = unpack("l", $c);
        return [$x, $y];
    };
    return {getCoordinates=>$getCoordinates, handle=>$positionsFile};
}

sub getFilter($)
{
    my ($self) = @_;
    my $buffer;
    my $r = read($self->{filtersFile}, $buffer, 1);
    croak "Failed to read filter information: $!" unless $r and 1 == $r;
    # TODO: add full support for the controls
    my $fullFilter = unpack('C', $buffer);
    my $passFilter = ($fullFilter & 1);
    my $isControl = ($fullFilter > 1);
    return ($isControl ? 0 : $passFilter);
}

1;                                  # says use was ok
__END__

