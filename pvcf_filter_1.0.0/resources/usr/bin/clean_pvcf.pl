#!/usr/bin/env perl

use strict;
use warnings;
use vars qw(@headers);

die("\nUSAGE: $0 <pVCF>\n\n") unless (@ARGV == 1);

my $pvcf = shift;

open(O,">filter.log");
open(P,"gunzip -c $pvcf |");
while(<P>) {
    if (/^\#/) {
	print;
	if (/^\#CHROM/) {
	    chomp;
	    @headers = split /\t/;
	}
	next;
    }
    chomp;
    my $adindex = -1;
    my $gqindex = -1;
    my $dpindex = -1;
    my @t = split /\t/;
    my @formats = split(":",$t[8]);
    for my $i (0..$#formats) {
	$adindex = $i, if ($formats[$i] eq "AD");
	$gqindex = $i, if ($formats[$i] eq "GQ");
	$dpindex = $i, if ($formats[$i] eq "DP");
    }

    my @lineculprits = ();
    for my $j (9..$#t) {
	my @culprits = ();
	my @fields = split(":",$t[$j]);
	next unless (scalar(@fields) == scalar(@formats));

	if ($adindex != -1) {
            my @genotypes = split("/",$fields[0]);
            my @adfields = split(",",$fields[$adindex]);
            if (scalar(@adfields) == 2) {
		next if ($genotypes[0] eq ".");
                if ($genotypes[0] != $genotypes[1]) { # heterozygous call
                    if ($adfields[0] != 0 && $adfields[1] != 0) {
                        my $allelefrequency = $adfields[0]/$adfields[1];
                        if ($allelefrequency > 4 || $allelefrequency < 0.25) { # heterozygous call needs to have a lower ratio than 80:20
                            push @culprits, "AD:$adfields[0],$adfields[1]:".sprintf("%.3f",$allelefrequency);
                        }
                    } else { # one AD is zero for a heterozygous call is an offending case
                        push @culprits, "AD:$adfields[0],$adfields[1]";
                    }
                } else { # homozygous call
                    if ($adfields[0] != 0 && $adfields[1] != 0) {
                        my $allelefrequency = $adfields[0]/$adfields[1];
                        if ($allelefrequency < 4 || $allelefrequency > 0.25) { # homozygous call needs to have a higher ratio than 80:20
                            push @culprits, "AD:$adfields[0],$adfields[1]:".sprintf("%.3f",$allelefrequency);
                        }
                    }
                }
            }
	}

	if ($gqindex != -1 && $fields[$gqindex] < 30) {
	    push @culprits, "GQ:".$fields[$gqindex];
	}
	if ($dpindex != -1 && ($fields[$dpindex] eq "." || $fields[$dpindex] < 8)) {
	    push @culprits, "DP:".$fields[$dpindex];
	}
	
	if (scalar(@culprits) > 0) {
	    $fields[0] = "./.";
	    push @lineculprits, join("; ",($headers[$j],$j+1,@culprits));
	    $t[$j] = join(":",@fields);
	}
    }
    print O join("\t",(@t[0..4],scalar(@lineculprits),@lineculprits)),"\n", if (@lineculprits > 0);

    print join("\t",@t),"\n";
}
close P;
close O; 
