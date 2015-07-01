package Casava::Common::Message;

BEGIN {
    use Exporter();
    @ISA       = qw(Exporter);
    @EXPORT = qw(&info &warning &error);    
    @EXPORT_OK = qw();
}
# PROJECT: GERALD
# MODULE:  Message.pm
# AUTHOR:  A. J. Cox
#
# Copyright (c) 2007 Solexa, 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

# Functions common to the GERALD module

=pod

=head1 NAME

Common::Message.pm - Messaging Functions 

=head2 SYNOPSIS

use Common::Message.pm qw();  

=head2 AUTHORSHIP

Copyright (c) 2007 Solexa, 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

=head1 DESCRIPTION

=head2 Overview

Messaging Functions 

=head2 Exports

    info($)
    warning($)
    error($)

Global variables:

=head2 Depends on

    warnings, strict, POSIX, XML::Simple, constant

=cut

use warnings FATAL => 'all';
use strict;
use POSIX;
use XML::Simple;
use constant DEBUG => 0;    # set to 1 to get debug info

sub info($);
sub warning($);
sub error($);

=pod

=head1 FUNCTIONS

=head2 General Functions

=over 4

=cut

=pod

=item info($msg)

The procedure prints message 

B<Parameters:>

    $msg            - Message

B<Returns:> 

    Nothing

=cut

sub info($) {
    my $msg = shift;
    print STDERR $msg;
} # info

=pod

=item warning($msg)

The procedure prints warning message 

B<Parameters:>

    $msg            - warning message

B<Returns:> 

    Nothing

=cut
sub warning($) {
    my $msg = shift;
    print STDERR "Warning: ", $msg;
} # warning

=pod

=item error($msg)

The procedure prints error message 

B<Parameters:>

    $msg            - error message

B<Returns:> 

    Nothing

=cut
sub error($) {
    my $msg = shift;
    die "Error: ", $msg;
} # sub error
1;    # says use was ok
__END__

