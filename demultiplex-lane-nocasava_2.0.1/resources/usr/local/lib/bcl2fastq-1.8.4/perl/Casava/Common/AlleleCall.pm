package Casava::Common::AlleleCall;

BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT = qw($FIELD_SEPARATOR %def &new &arrayRef &toString);
    @EXPORT_OK = qw();
}

# PROJECT: CASAVA
# MODULE:  AlleCall.pm
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

# AlleCall Class functions

=pod

=head1 NAME

Common::AlleCallRead.pm - Messaging Functions 

=head2 SYNOPSIS

use Casava::Common::AlleCall.pm qw();  

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

Read sort.count file 

=head2 Overview

Definition of the _qseq data format

Currently tab separated text format.
 Position  A  C  G  T  calledBase  Total_count
                # Called_bases_count next_best_call score
                
    * Position: position on the chromosome/sccaffold
    * A count: number of A bases
    * C count: number of C bases
    * G count: number of G bases
    * T count: number of T bases
    * base: called base 
    * total: Total number of bases
    * called: Number of called bases
    * Score: ":" seperated score

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

## format fields

our $FIELD_SEPARATOR = '\s+';

our %def             = (
    position => 0,
    Acount   => 1,
    Ccount   => 2,
    Gcount   => 3,
    Tcount   => 4,
    base     => 5,
    total    => 6,
    called   => 7,
    score    => 8,
);

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
        foreach my $key ( keys %def ) {
            push( @readArray, undef );
        }    # foreach
        $self = { _readArrayRef => \@readArray };
    }

    bless $self, $class;
    return $self;
}

=pod

=item position($value)

Accessor method for position field

=cut

sub position($;\@) {
    my ( $self, $position ) = @_;
    $self->{_readArrayRef}[ $def{position} ] = $position
      if defined($position);
    return $self->{_readArrayRef}[ $def{position} ];
}

=pod

=item Acount($value)

Accessor method for number of A bases

=cut

sub Acount {
    my ( $self, $acount ) = @_;
    $self->{_readArrayRef}[ $def{Acount} ] = $acount    
      if defined($acount);
    return $self->{_readArrayRef}[ $def{Acount} ];
}

=pod

=item Acount($value)

Accessor method for number of C bases

=cut

sub Ccount {
    my ( $self, $ccount ) = @_;
    $self->{_readArrayRef}[ $def{Ccount} ] = $ccount    
      if defined($ccount);
    return $self->{_readArrayRef}[ $def{Ccount} ];
}
=pod

=item Gcount($value)

Accessor method for number of G bases

=cut

sub Gcount {
    my ( $self, $gcount ) = @_;
    $self->{_readArrayRef}[ $def{Gcount} ] = $gcount    
      if defined($gcount);
    return $self->{_readArrayRef}[ $def{Gcount} ];
}

=pod

=item Tcount($value)

Accessor method for number of T bases

=cut

sub Tcount {
    my ( $self, $tcount ) = @_;
    $self->{_readArrayRef}[ $def{Tcount} ] = $tcount    
      if defined($tcount);
    return $self->{_readArrayRef}[ $def{Tcount} ];
}

=pod

=item base($value)

Accessor method for base called

=cut
sub base {
    my ( $self, $base ) = @_;
    $self->{_readArrayRef}[ $def{base} ] = $base    
      if defined($base);
    return $self->{_readArrayRef}[ $def{base} ];
}

=pod

=item total($value)

Accessor method for total number of bases

=cut

sub total {
    my ( $self, $total ) = @_;
    $self->{_readArrayRef}[ $def{total} ] = $total    
      if defined($total);
    return $self->{_readArrayRef}[ $def{total} ];
}

=pod

=item called($value)

Accessor method for number of called bases

=cut

sub called {
    my ( $self, $called ) = @_;
    $self->{_readArrayRef}[ $def{called} ] = $called    
      if defined($called);
    return $self->{_readArrayRef}[ $def{called} ];
}

=pod

=item score($value)

Accessor method for allele call score

=cut

sub score {
    my ( $self, $score ) = @_;
    $self->{_readArrayRef}[ $def{score} ] = $score    
      if defined($score);
    return $self->{_readArrayRef}[ $def{score} ];
}

=pod

=item arrayRef($value)

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

The method prints AlleCall to string

=cut

sub toString {
    my ($self) = @_;
    return join( $FIELD_SEPARATOR, @{ $self->{_readArrayRef} } );
}
1;    # says use was ok
__END__

