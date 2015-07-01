#!/usr/bin/env perl

# To get the sequence information from BAM and output in FASTQ

use strict;
use warnings;

die("\nUSAGE: $0 <BAM file> <basename>
BAM file needs to be sorted by read names.\n\n") unless (@ARGV == 2);

my ($bamfile,$basename) = @ARGV;

open(I,"samtools view $bamfile |");
open(Q1,">$basename.R1.fq");
open(Q2,">$basename.R2.fq");
while(<I>) {
    chomp;
    my @t = split /\t/;
    my $bitflag = $t[1];
    next if ($bitflag & 256);

    my $id = $t[0];
    my $seq = $t[9];
    my $qual = $t[10];
    ($seq,$qual) = &reverse($seq,$qual), if ($bitflag & 16);

    if ($bitflag & 64) {
	print Q1 "@",$id," 1:N:0\n",$seq,"\n+\n",$qual,"\n";
    }
    if ($bitflag & 128) {
	print Q2 "@",$id," 2:N:0\n",$seq,"\n+\n",$qual,"\n";
    }
}
close I;
close Q1;
close Q2;
    
sub reverse {
    my ($seq,$qual) = @_;
    $seq = reverse($seq);
    $seq =~ tr/ACGTN/TGCAN/;
    $qual = reverse($qual);
    return ($seq,$qual);
}
