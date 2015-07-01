#!/usr/bin/env perl

# To generate a master GVCF

use strict;
use warnings;
#use lib '/usr/local/lib';
use Bio::SeqIO;
use vars qw(%seqs %chrinfo);

die("\nUSAGE: $0 <EVE_file> <freeze_project_name> <seqfile> <VCRome_region>\n\n") unless (@ARGV == 4);

my ($evefile,$projectname,$seqfile,$vcrregionfile) = @ARGV;

#my $seqfile = "/home/ubuntu/databases/hs37d5/genome.fa";
my $seqin = Bio::SeqIO->new(-file => $seqfile, -format => 'fasta');
while(my $seqobj = $seqin->next_seq) {
    my $id = $seqobj->id;
    $seqs{$id} = $seqobj;
}
$seqin->close;

my $tmpvcf = "$projectname.mock.gvcf";
if ($evefile =~ /\.gz$/) {
    open(E,"gunzip -c $evefile |");
} else {
    open(E,"$evefile");
}
open(T,">$tmpvcf");
while(<E>) {
    chomp;
    next if /^\#/;
    my @t = split /\t/;
    $t[2] = ".";
    $t[4] = $t[4].",<NON_REF>";
    $t[7] = ".";
    print T join("\t",(@t,"GT:AD:DP:GQ:PL:SB","./.:0,16,0:16:.:476,48,0,476,48,476:0,0,16,0")),"\n";
}
close E;
close T;

system("cat $vcrregionfile $projectname.mock.gvcf | sort -k1,1 -k2,2n -k3,3 > $projectname.vcrb37.combined.txt");

my $infile = "$projectname.vcrb37.combined.txt";
open(I,"$infile");
my $buffer;
read(I,$buffer,-s $infile);
my @lines = split("\n",$buffer);

foreach my $i (0..$#lines) {
    my $line1 = $lines[$i];
    my @line1elems = split("\t",$line1);
    next if ($line1elems[2] eq "#"); # Not to process the two lines across two regions
    my $line2 = $lines[$i+1];
    my @line2elems = split("\t",$line2);
    next if ($line1elems[0] ne $line2elems[0]); # Not to process the two lines across two chromosomes

    my $coord1 = $line1elems[1];
    my $coord2 = $line2elems[1];
    my $coord1next= $coord1+1;
    if ($coord1 == $coord2 || $coord2 == $coord1next) { # Not to process two lines with the same or immediate adjacent coordinates
	push @{$chrinfo{$line2elems[0]}}, join("\t",@line2elems), if (@line2elems > 3);
	next;
    }

    my $regionend = $coord2 - 1;
    if ($line2elems[2] eq "#") {
	$regionend = $coord2;
    }

    my $chrobj = $seqs{$line1elems[0]};

    # output the gvcf line with the reference base and the blocks without an alternate variant
    if ($line1elems[2] eq "!") {
	my $refbase = $chrobj->subseq($coord1,$coord1);
	push @{$chrinfo{$line2elems[0]}}, join("\t",($line1elems[0],$coord1,".",$refbase,"<NON_REF>",".",".","END=$regionend","GT:DP:GQ:MIN_DP:PL","./.:100:99:100:0,0,1000"));
    } elsif ($regionend > $coord1+1) {
	my $refbase = $chrobj->subseq($coord1+1,$coord1+1);
	push @{$chrinfo{$line2elems[0]}}, join("\t",($line1elems[0],$coord1+1,".",$refbase,"<NON_REF>",".",".","END=$regionend","GT:DP:GQ:MIN_DP:PL","./.:100:99:100:0,0,1000"));
    }
    
    next unless ($line2elems[3]);
    push @{$chrinfo{$line2elems[0]}}, join("\t",@line2elems); # output the gvcf line with an alternate variant
}

open(O,">$projectname.master.gvcf");
print O "##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
";
print O join("\t",("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","template")),"\n";
my @chromosomes = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
for my $chr (@chromosomes) {
    print O join("\n",@{$chrinfo{$chr}}),"\n";
}
close O;
