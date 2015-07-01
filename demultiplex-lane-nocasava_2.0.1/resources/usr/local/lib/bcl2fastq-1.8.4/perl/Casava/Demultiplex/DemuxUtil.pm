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

Casava::Demultiplex::DemuxUtil - Simple utilities to synchronize
demux/fastq converter routines.

=cut

# POD for for each subroutines is right before the implementation
# More POD after __END__

package Casava::Demultiplex::DemuxUtil;


BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT =
      qw($fastqSuffix $compressedFastqSuffix checkMakeDir defaultSampleName defaultSampleDirName defaultFastqPrefix);
}



use strict;
use warnings "all";

use Carp;
use File::Path qw(mkpath);


our $fastqSuffix = ".fastq";
our $compressedFastqSuffix = ".fastq.gz";


sub checkMakeDir($) {
    my $dir = shift;
    unless (-e $dir) {
        mkpath($dir) || croak("ERROR: Can't create directory '$dir'\n");
    } else {
        croak "ERROR: Path is not a directory '$dir'\n"  unless -d $dir;
    }
}



sub defaultSampleName($) {
    my $lane = int(shift);
    return "lane".$lane;
}



my $sampleDirPrefix="Sample_";

sub defaultSampleDirName($) {
    my $lane = shift;
    return $sampleDirPrefix . defaultSampleName($lane);
}



sub defaultFastqPrefix($$) {
    my $lane = shift;
    my $read = shift;

    my $sampleName = defaultSampleName($lane);
    my $padLane = sprintf "%03d", $lane;
    return "$sampleName\_NoIndex_L$padLane\_R$read\_";
}



1;
__END__

=pod

=cut
