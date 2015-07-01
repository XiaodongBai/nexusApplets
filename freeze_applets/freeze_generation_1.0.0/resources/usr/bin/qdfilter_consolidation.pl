#!/usr/bin/env perl

# To consolidate the information of QDfilter

use strict;
use warnings;
use vars qw(%hash %sites %blacklist);

die("\nUSAGE: $0 <header> <QD> <blacklist>\n\n") unless (@ARGV == 3);

my ($headerfile,$qdfile,$blacklistfile) = @ARGV;

open(B,"$blacklistfile");
while(<B>) {
    chomp;
    $blacklist{$_} = 1;
}
close B;

open(H,"$headerfile");
chomp(my $line = <H>);
my @t = split(/\t/,$line);
for my $i (9..$#t) {
    my $index = $i-9;
    $hash{$t[$i]} = $index;
}
close H;

open(Q,"$qdfile");
while(<Q>) {
    chomp;
    my @t = split /\t/;
    if ($t[0] eq "X") {
	$t[0] = 23;
    } elsif ($t[0] eq "Y") {
	$t[0] = 24;
    }
    my $name = join(":",@t[0..1]);
    for my $j (4..$#t) {
	next if ($blacklist{$t[$j]});
	$sites{$name}{$hash{$t[$j]}} = 1;
    }
}
close Q;

foreach my $site (sort keys %sites) {
    my $samplesites = $sites{$site};
    print join("\t",($site,join(":",sort {$a <=> $b} keys %$samplesites))),"\n";
}
