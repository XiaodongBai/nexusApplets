package  Casava::Common::TabFile;

BEGIN {
    use Exporter();
    @ISA       = qw(Exporter);
    @EXPORT    = qw(&openFile &closeFile &getNext &append2);
    @EXPORT_OK = qw();
}

# PROJECT: CASAVA
# MODULE:  TabFile.pm
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

# Generic functions to reaf/write AlleleCall and SNP (sort.count and snp.txt) files

=pod

=head1 NAME

Casava::Common::TabFile.pm - Generic API functions for tab delimited files

=head2 SYNOPSIS

use Casava::Common::TabFile.pm qw();  

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

Generic functions to reaf/write AlleleCall and SNP (sort.count and snp.txt) files

=head2 Overview


=head2 Exports

    append2(\%;\%);
    getNext(\%);
    closeFile(\%);
    openFile($;$;$);

Global variables:

=head2 Depends on
    
    warnings, strict, POSIX, constant, IO::File

=cut

use warnings FATAL => 'all';
use strict;
use POSIX;
use IO::File;
use Carp;
use Data::Dumper;
use constant DEBUG => 0;    # set to 1 to get debug info
use Casava::Common::Log;

sub append2(\%;\%);
sub getNext(\%);
sub closeFile(\%);
sub openFile($;$;$);

=pod

=head1 FUNCTIONS

=head2 General Functions

=over 4

=cut

=pod

=item openFile($;$)

The procedure opens file type of $class and returns 
    ref to a structure destribing file handle. 

B<Parameters:>

    $filePath             - path to a file
    $mode                 - > or < (read or write)
    $class                - File class (
                                Casava::Common::AlleleCall | 
                                Casava::Common::SNP)

B<Returns:> 

    HASH MAP Ref to Qseq structure:
    file
    filePath

=cut

sub openFile($;$;$) {
    croak("ERROR: openFile wrong number of parameters")
      unless ( @_ == 3 );
    my ( $filePath, $mode, $class ) = @_;
    my %handleMap = ();
    print "$filePath\n";
    my $file      = IO::File->new("$mode$filePath")
      or errorExit "$0:ERROR: Couldn't open $mode$filePath $!\n";
    $handleMap{mode}           = $mode;
    $handleMap{class}          = $class;
    $handleMap{readCount}      = 0;
    $handleMap{fieldCount}     = scalar( eval( 'keys %' . $class . '::def' ) );
    $handleMap{fieldSeparator} = eval( '$' . $class . '::FIELD_SEPARATOR' );
    $handleMap{file}     = $file;       
    $handleMap{filePath} = $filePath;
    return \%handleMap;
}    # sub openQseq

=pod

=item closeFile(\%) 

    The procedure closes the file

B<Parameters:>

    $handleMapRef         - HASH MAP Ref to structure destribing file handle

B<Returns:> 

    status (0 or 1)

=cut

sub closeFile(\%) {
    croak("ERROR: closeFile wrong number of parameters")
      unless ( @_ == 1 );
    my ($handleMapRef) = @_;
    my $file = $handleMapRef->{file};
    close($file);
    return 1;
}    # sub close

=pod

=item getNext(\%)

The procedure reads one row (ARRAY Ref) from file. If end of file
then returns undef.
To iterate

while ( defined( my $itemRef = getNext( %{$handleMapRef} ) ) ) {
    my $content = $itemRef->toString();
}

B<Parameters:>

    $handleMapRef         - HASH MAP Ref to structure destribing file handle

B<Returns:> 

    Class object Ref

=cut

sub getNext(\%) {
    croak("ERROR: getQseqRead wrong number of parameters")
      unless ( @_ == 1 );
    my ($handleMapRef) = @_;
    my $file           = $handleMapRef->{file};
    my $line           = <$file>;
    while ( defined $line && $line ne '' && $line =~ m/^#/) {
        $line = <$file>;
    }

    if ( !defined $line || $line eq '' ) {
        return undef;
    }

    chomp($line);
#    print "$line\n$handleMapRef->{fieldSeparator}";
    my @item = split( $handleMapRef->{fieldSeparator}, $line );
    if ( scalar(@item) != $handleMapRef->{fieldCount} ) {
        my $fieldsCount = scalar(@item);
        croak
"ERROR: getNext wrong number of fields [$fieldsCount] expected [$handleMapRef->{fieldCount}] in [$line]\n";
    }
    $handleMapRef->{count}++;
    return $handleMapRef->{class}->new( \@item );
}    # sub getNext

=pod

=item append2(\%)

The procedure append one item as a string (ARRAY Ref) to a file. 

B<Parameters:>

    $handleMapRef        - HASH MAP Ref to structure destribing file handle
    $itemRef             - Item object Ref

B<Returns:> 

    Nothing

=cut

sub append2(\%;\%) {
    croak("ERROR: append2 wrong number of parameters")
      unless ( @_ == 2 );
    my ( $qseqReadRef, $itemRef ) = @_;
    my $file = $qseqReadRef->{file};
    print $file $itemRef->toString() . "\n";
}     # sub append
1;    # says use was ok
__END__
