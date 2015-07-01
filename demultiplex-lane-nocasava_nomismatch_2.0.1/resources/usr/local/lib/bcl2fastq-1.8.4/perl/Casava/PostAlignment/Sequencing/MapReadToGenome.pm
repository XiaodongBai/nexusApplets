package Casava::PostAlignment::Sequencing::MapReadToGenome;

BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT =
      qw(&mapReadToGenome $VERSION);
}
# PROJECT: CASAVA
# MODULE:  $RCSfile: MapReadToGenome.pm,v $
# File: mapReadToGenome.pl
#
# Copyright (c) 2006,2007 Solexa
# Copyright (c) 2007-2010 Illumina, Inc.
#
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).
#
# This file is part of the Consensus Assessment of Sequence And VAriation
# (CASAVA) software package.
#
# Author: A. J. Cox
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

# Description:
# Convert a read into a set of bases that can be mapped to a reference
# file, accounting for:
# 1. orientation of its match relative to reference
# 2. the way ELAND handles N characters
# Convert the reads associated fastq string to match

# Usage:
# mapReadToGenome( read, fastq, nDescription, fastqString

use strict;
use English;

sub mapReadToGenome
{
    my ($read, $strand, $nDescription, $fastq)=@_;
    my ($r2, $q2);
    my ($leadingNs, $trailingNs);
    my ($leadingQs, $trailingQs);

    $fastq='0' x length($read) unless (defined($fastq));

    die "length of read $read (".length($read).
    ") not equal to length of quality string $fastq (".length($fastq).
    ")\n" unless (length($fastq)==length($read));

    # This code changed by TC 22.4.8 to fix issue IMP-15
    # Function needs to retain leading Ns

    # Trim off trailing Ns
#    $read=~s/[Nn]+$//;
#    $fastq=substr($fastq,0,length($read));

    # Trim off leading Ns
#    $read=~s/^[Nn]+//;
#    $fastq=substr($fastq,length($fastq)-length($read));

    # separate leading and trailing Ns from read and split
    # corresponding quality string in the same way
    if ($read=~/(^[Nn]*)(.*?)([Nn]*$)/)
    {
    ($leadingNs, $read, $trailingNs)=($1,$2,$3);
    $leadingQs=substr($fastq,0,length($leadingNs));
    $trailingQs=substr($fastq,-length($trailingNs),length($trailingNs));
    $fastq=substr($fastq,
              length($leadingNs),
              length($fastq)-length($leadingNs)-length($trailingNs));
#   print "$leadingNs $read $trailingNs $leadingQs $fastq $trailingQs\n";
    } # if

    
    # if an N is an I, delete the base (insertion in read relative to genome)
    # if an N is a D, leave it in (base read but not detected)
    for my $nType (split(//,$nDescription))
    {
    next unless ($nType=~/[DI]/);
    if ($read=~/([Nn])/)
    {
        $r2.=$PREMATCH;
        $q2.=substr($fastq,0,length($PREMATCH));
        if ($nType eq 'D')
        {
        $r2.=$1;
        $q2.=substr($fastq,length($PREMATCH),1);
        } # if
        $read=$POSTMATCH;
        $fastq=substr($fastq,length($fastq)-length($read));
    }
    else
    {
        print STDERR "Expected to find N";
    } # else
    } # for

    $r2.=$read;
    $q2.=$fastq;
    
    # Finally, change strand if required

    if ($strand eq 'R')
    {
    # Append trailing Ns, these become leading Ns when translated to 
    # forward strand
    $r2.=$trailingNs;
    $q2.=$trailingQs;
    $r2=reverse($r2);
    $r2=~tr/ACGTacgt/TGCAtgca/;
    $q2=reverse($q2);
    } # if
    else
    {
    # Add in leading Ns
    $r2=$leadingNs.$r2;
    $q2=$leadingQs.$q2;
    } # else

    die "length of mapped read $r2 (".length($r2).
    ") not equal to length of quality string $q2 (".length($q2).
    ")\nOriginal read = $read original quality=$fastq\n" 
    unless (length($q2)==length($r2));
    return ($r2,$q2);
} # sub mapReadToGenome


# Comment out next 2 lines to use as standalone 
1;
__END__
# Simple main program to test

my ($read, $strand, $nDescription, $fastq)=@ARGV;

die unless ((@ARGV==3)||(@ARGV==4));

print "$read $fastq\n";

my ($readOut, $fastqOut)=mapReadToGenome($read, $strand, $nDescription, $fastq);
print "$readOut $fastqOut\n";

