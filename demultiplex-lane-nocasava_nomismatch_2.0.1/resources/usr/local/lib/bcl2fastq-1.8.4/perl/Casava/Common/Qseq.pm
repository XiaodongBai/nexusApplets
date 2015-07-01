package Casava::Common::Qseq;

BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT = qw(%qseqDef $QSEQ_FIELD_SEPARATOR);
    @EXPORT_OK = qw(&new &close &getRead &writeQseqRead);
}

# PROJECT: Pipeline
# MODULE:  Qseq.pm
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

# _qseq file API functions

=pod

=head1 NAME

Casava::Common::Qseq.pm - _qseq file API functions

=head2 SYNOPSIS

use Casava::Common::Qseq.pm qw();  

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

Qseq Read file used between Bustard and Gerald (_qseq file API functions)

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

    new($$$);
    close($);
    getRead($);
    writeQseqRead($$);

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

sub new($$$);
sub close($);
sub getRead($);
sub writeQseqRead($$);

=pod

=head1 FUNCTIONS

=head2 General Functions

=over 4

=cut

=pod

=item new($$$)

The procedure opens qseq file and returns ref to a structure destribing 
    qseq file. 

B<Parameters:>

    $filePath             - path to _qseq.txt file

B<Returns:> 

    HASH MAP Ref to Qseq structure:
    file
    filePath

=cut

sub new($$$) {
    my $class = shift;
    my ( $filePath, $mode ) = @_;
    die "Undefined file path" unless defined $filePath;
    die "Undefined open-mode" unless defined $mode;
    my $file = IO::File->new("$mode$filePath") or die "$0:ERROR: Couldn't open $mode$filePath $!\n";
    my $self = {mode=>$mode, readCount=>0, fieldCount=>scalar( keys %qseqDef ), file=>$file, filePath=>$filePath};
    bless $self, $class;
    return $self;
}    # sub new

=pod

=item close($) 

    The procedure closes qseq file

B<Parameters:>

    $qseqRef             - HASH MAP Ref to Qseq structure:

B<Returns:> 

    status (0 or 1)

=cut

sub close($) {
    my ($qseqRef) = @_;
    my $file = $qseqRef->{file};
    CORE::close($file);
    return 1;
}    # sub closeQseq

=pod

=item getRead($)

The procedure reads one read (ARRAY Ref) from qseq file. If end of file
then returns undef.
To iterate

while ( defined( my $readTmpRef = getQseqRead( %{$qseqFileRef} ) ) ) {
    my $tile = $readTmpRef->tile();
}

B<Parameters:>

    $qseqRef             - HASH MAP Ref to Qseq structure see openFile

B<Returns:> 

    QseqRead object Ref

=cut

sub getRead($) {
    my ($qseqRef) = @_;
    my $file      = $qseqRef->{file};
    my $line      = <$file>;
    if ( !defined $line || $line eq '' ) {
        return undef;
    }
    chomp($line);
    my @read = split $QSEQ_FIELD_SEPARATOR, $line;
    if ( scalar(@read) != $qseqRef->{fieldCount} ) {
        my $fieldsCount = scalar(@read);
        die
"ERROR: getQseqRead wrong number of fields [$fieldsCount] expected [$qseqRef->{fieldCount}] in [$line]\n";
    }
    $qseqRef->{readCount}++;
    return Casava::Common::QseqRead->new( \@read );
}    # sub getQseqRead

=pod

=item writeQseqRead($$)

The procedure writes one qseq read (ARRAY Ref) to qseq file. 

B<Parameters:>

    $qseqReadRef         - ARRAY MAP Ref to qseq read structure see getQseqRead
    $qseqRef             - QseqRead object Ref

B<Returns:> 

    Nothing

=cut

sub writeQseqRead($$) {
    my ( $qseqReadRef, $qseqRef ) = @_;
    my $file = $qseqRef->{file};
    print $file $qseqReadRef->toString() . "\n";    
}    # sub writeQseqRead
1;                                  # says use was ok
__END__

