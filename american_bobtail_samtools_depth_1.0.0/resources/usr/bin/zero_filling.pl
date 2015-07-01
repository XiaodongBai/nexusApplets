#!/usr/bin/env perl

use strict;
use warnings;
use vars qw(%depth);

die("\nUSAGE: $0 <target bed> <nonzero depth>\n\n") unless (@ARGV == 2);

my ($targetbed,$non0depth) = @ARGV;

open(N,"$non0depth");
while(<N>) {
    chomp;
    my @t = split /\t/;
    my $name = join("\t",@t[0..1]);
    $depth{$name} = $t[2];
}
close N;

open(T,"$targetbed");
while(<T>) {
    chomp;
    my @t = split /\t/;
    my $start = $t[1]+1;
    my $end = $t[2];
    for my $i ($start..$end) {
	my $name = join("\t",($t[0],$i));
	if ($depth{$name}) {
	    print join("\t",($name,$depth{$name})),"\n";
	} else {
	    print join("\t",($name,0)),"\n";
	}
    }
}
close T;
