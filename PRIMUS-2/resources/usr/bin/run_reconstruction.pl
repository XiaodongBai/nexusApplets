#! /usr/bin/perl

use lib '/usr/lib/perl_modules';
use lib '/home/ubuntu/applets/PRIMUS/resources/usr/lib/perl_modules';
use strict;
use PRIMUS::reconstruct_pedigree_v7;

#$ENV{'PERL5LIB'} =~ s/5\.14\.2/5\.10\.1/g;	

PRIMUS::reconstruct_pedigree_v7::reconstruct_pedigree(@ARGV);
