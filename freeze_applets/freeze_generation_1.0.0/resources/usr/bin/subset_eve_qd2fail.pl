#!/usr/bin/env perl

# To get the subset of the eve_qd2fail file

use strict;
use warnings;

die("\nUSAGE: $0 <file> <sample_code>\n\n") unless (@ARGV == 2);

my ($file,$code) = @ARGV;

open(F,"$file");
while(<F>) {
    chomp;
    next unless (/$code/);
    my @t = split /\t/;
    my @samples = ();
    for my $i (4..$#t) {
	push @samples, $t[$i], if ($t[$i] =~ /^$code/);
    }
    print join("\t",(@t[0..3],@samples)),"\n";
}
close F;
