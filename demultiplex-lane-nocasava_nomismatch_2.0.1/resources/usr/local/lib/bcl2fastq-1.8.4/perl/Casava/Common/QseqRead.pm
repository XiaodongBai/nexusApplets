package Casava::Common::QseqRead;

BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT = qw($QSEQ_FIELD_SEPARATOR %qseqDef 
      &machineName &runNumber &lane &tile &X &Y
      &Index &seq &qualityString &filter &arrayRef &toString);
    @EXPORT_OK = qw(&new);
}

# PROJECT: Pipeline
# MODULE:  Qseq.pm
# AUTHOR:  Lukasz Szajkowski
#

# Qseq Class functions

=pod

=head1 NAME

Common::QseqRead.pm - Messaging Functions 

=head2 SYNOPSIS

use Common::QseqRead.pm qw();  

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

Read file used between Bustard and Gerald

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

1 file for 1 "virtual lane" (concatenation of the lane number on 1 digit, with the index number on 4 digits), names s_lllll_seq.txt, where "lllll" is the id of the virtual lane.

Note: to avoid problems with concurrent access to the files, intermediary "tile-based _qseq files" are created. When they are all completed, they are dispatched by index number and concatenated to get the final _qseq file.
Messaging Functions 

=head2 Exports


Global variables:

=head2 Depends on

    warnings, strict, POSIX, constant

=cut

use warnings FATAL => 'all';
use strict;
use POSIX;
use IO::File;
use constant DEBUG => 0;    # set to 1 to get debug info

## _qseq format fields (try to keep it as sub class of _expot.txt format)
our $QSEQ_FIELD_SEPARATOR = "\t";
our %qseqDef = (
    MachineName   => 0,
    RunNumber     => 1,
    Lane          => 2,
    Tile          => 3,
    X             => 4,
    Y             => 5,
    Index         => 6,
    ReadNum       => 7,
    Seq           => 8,
    QualityString => 9,
    Filter        => 10,
);

=pod

=head1 METHODS

=head2 General Methods

=over 4

=cut

=pod

=head1 NAME

Common::Read.pm - Read class

=head2 SYNOPSIS

use Common::Read.pm qw();  

=head1 DESCRIPTION

=head2 Overview

    Classs with fieles
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
    
Messaging Functions 

=head2 Exports


Global variables:

=head2 Depends on

    warnings, strict, POSIX, XML::Simple, constant

=cut

=pod

=head1 METHODS

=head2 General Methods

=over 4

=cut

#constructor

=pod

=item new

construtor

=cut

sub new {
    my ( $class, $readArrayRef ) = @_;
    my $self;
    if ( defined $readArrayRef ) {
        $self = { _readArrayRef => $readArrayRef };
    }
    else {
        my @readArray = ();
        foreach my $key ( keys %qseqDef ) {
            push( @readArray, undef );
        }    # foreach
        $self = { _readArrayRef => \@readArray };
    }

    bless $self, $class;
    return $self;
}

=pod

=item machineName($value)

Accessor method for machine name field

=cut

sub machineName {
    my ( $self, $machine ) = @_;
    $self->{_readArrayRef}[ $qseqDef{MachineName} ] = $machine
      if defined($machine);
    return $self->{_readArrayRef}[ $qseqDef{MachineName} ];
}

=pod

=item runNumber($value)

Accessor method for run number

=cut

sub runNumber {
    my ( $self, $runNumber ) = @_;
    $self->{_readArrayRef}[ $qseqDef{RunNumber} ] = $runNumber
      if defined($runNumber);
    return $self->{_readArrayRef}[ $qseqDef{RunNumber} ];
}

=pod

=item lane($value)

Accessor method for line number

=cut

sub lane {
    my ( $self, $lane ) = @_;
    $self->{_readArrayRef}[ $qseqDef{Lane} ] = $lane if defined($lane);
    return $self->{_readArrayRef}[ $qseqDef{Lane} ];
}

=pod

=item tile($value)

Accessor method for tile number

=cut

sub tile {
    my ( $self, $tile ) = @_;
    $self->{_readArrayRef}[ $qseqDef{Tile} ] = $tile if defined($tile);
    return $self->{_readArrayRef}[ $qseqDef{Tile} ];
}

=pod

=item x($value)

Accessor method for x

=cut

sub X {
    my ( $self, $x ) = @_;
    $self->{_readArrayRef}[ $qseqDef{X} ] = $x if defined($x);
    return $self->{_readArrayRef}[ $qseqDef{X} ];
}

=pod

=item Y($value)

Accessor method for y

=cut

sub Y {
    my ( $self, $y ) = @_;
    $self->{_readArrayRef}[ $qseqDef{Y} ] = $y if defined($y);
    return $self->{_readArrayRef}[ $qseqDef{Y} ];
}

=pod

=item index($value)

Accessor method for index

=cut

sub Index {
    my ( $self, $index ) = @_;
    $self->{_readArrayRef}[ $qseqDef{Index} ] = $index if defined($index);
    return $self->{_readArrayRef}[ $qseqDef{Index} ];
}

=pod

=item readNum($value)

Accessor method for read number

=cut

sub readNum {
    my ( $self, $readNum ) = @_;
    $self->{_readArrayRef}[ $qseqDef{ReadNum} ] = $readNum if defined($readNum);
    return $self->{_readArrayRef}[ $qseqDef{ReadNum} ];
}

=pod

=item seq($value)

Accessor method for sequence

=cut

sub seq {
    my ( $self, $seq ) = @_;
    $self->{_readArrayRef}[ $qseqDef{Seq} ] = $seq if defined($seq);
    return $self->{_readArrayRef}[ $qseqDef{Seq} ];
}

=pod

=item qualityString($value)

Accessor method for quality

=cut

sub qualityString {
    my ( $self, $qualityString ) = @_;
    $self->{_readArrayRef}[ $qseqDef{QualityString} ] = $qualityString
      if defined($qualityString);
    return $self->{_readArrayRef}[ $qseqDef{QualityString} ];
}

=pod

=item filter($value)

Accessor method for Filter

=cut

sub filter {
    my ( $self, $filter ) = @_;
    $self->{_readArrayRef}[ $qseqDef{Filter} ] = $filter
      if defined($filter);
    return $self->{_readArrayRef}[ $qseqDef{Filter} ];
}

=pod

=item filter($value)

Accessor method for readArray

=cut

sub arrayRef {
    my ( $self, $readArrayRef ) = @_;
    $self->{_readArrayRef} = $readArrayRef
      if defined($readArrayRef);    
    return $self->{_readArrayRef};
}

=pod

=item toString

The method prints Qseqread to string

=cut

sub toString {
    my ( $self ) = @_;    
    return join ($QSEQ_FIELD_SEPARATOR, @{$self->{_readArrayRef}});
}
1;                                  # says use was ok
__END__
