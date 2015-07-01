#!/usr/bin/env perl

# To calculate the base coverage on Y chromosome

use strict;
use warnings;

my $gvcffile = shift;
open(G,"$gvcffile");
my $total = 0;
while(<G>) {
    chomp;
    my @t = split /\t/;
    my @info = split(":",$t[8]);
    my $dpindex = 0;
    for my $i (0..$#info) {
	$dpindex = $i, if ($info[$i] eq "DP");
    }
    my @gpinfo = split(":",$t[9]);
    my $depth = $gpinfo[$dpindex];
    next if ($depth == 0);
    my $end = $t[1];
    if ($t[7] =~ /END/) {
	$t[7] =~ /END=(\d+)$/;
	$end = $1;
    }
    my $blocksize = $end-$t[1]+1;
    $total += $depth*$blocksize;
}
close G;

if ($total > 500000) {
    print "Male";
} else {
    print "Female";
}
