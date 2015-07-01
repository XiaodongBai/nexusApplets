#!/usr/bin/env perl

# To process the information file generated by vcf_eve_qdfilter

use strict;
use warnings;
use vars qw(%sites);

die("\nUSAGE: $0 <original>\n\n") unless (@ARGV == 1);

my ($orifile) = @ARGV;

if ($orifile =~ /\.gz/) {
    open(O,"gunzip -c $orifile |");
} else {
    open(O,"$orifile");
}
while(<O>) {
    chomp;
    s/\:$//;
    my @t = split /\t/;
    push @{$sites{$t[0]}}, $t[1];
}
close O;

foreach my $key (sort keys %sites) {
    print join("\t",($key,join(":",@{$sites{$key}}))),"\n";
}
