#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2007-2011 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And
VAriation (CASAVA) software package.

=cut

=head1 SYNOPSIS

finishFastq.pl

Helper function for batched Qseq To Fastq conversion which splits
per-tile fastq files into files with a maximum record count, gzip
compresses these and finally moves them from their tempoaray to final
location to indicate successful completion.

=cut

use warnings "all";
use strict;

use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Copy qw(move);
use Carp;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Demultiplex::DemuxUtil;


my $scriptName = (File::Spec->splitpath($0))[2];

my $argCount=scalar(@ARGV);
my $help;

my $sampleTmpDir;
my $sampleOutputDir;
my $lane;
my $read;
my $fastqClusterCount;

my $result = GetOptions('help|h' => \$help,
                        'sample-tmp-dir=s' => \$sampleTmpDir,
                        'sample-output-dir=s' => \$sampleOutputDir,
                        'lane=s' => \$lane,
                        'read=s' => \$read,
                        'fastq-cluster-count=i' => \$fastqClusterCount) or pod2usage(2);

if($help or (not $argCount)) {
    pod2usage(2);
}

croak "ERROR: sample output directory argument: '$sampleOutputDir' does not exist or is not a directory\n"
  unless (defined $sampleOutputDir) and (-d $sampleOutputDir);
croak "ERROR: sample temp directory argument: '$sampleTmpDir' does not exist or is not a directory\n"
  unless (defined $sampleTmpDir) and (-d $sampleTmpDir);
croak "ERROR: --fastq-cluster-count must be greater than 0\n"
  unless (defined $fastqClusterCount) &&  ($fastqClusterCount>0);
croak "ERROR: --lane must be defined.\n"
  unless (defined $lane);
croak "ERROR: --read must be defined.\n"
  unless (defined $read);

my @fastqs=@ARGV;

croak "ERROR: no fastqs specifed for splitting\n"
  unless (scalar(@fastqs));

# assume project directories are already created:
#

my $samplePath= $sampleOutputDir;
my $tempSamplePath= $sampleTmpDir;

my $fastqPrefix=defaultFastqPrefix($lane,$read);
my $fastqPathPrefix=File::Spec->catfile($tempSamplePath,$fastqPrefix);
my $finalPathPrefix=File::Spec->catfile($samplePath,$fastqPrefix);


checkMakeDir($tempSamplePath);
unlink(glob("$fastqPathPrefix*"));

my $splitNumber=0;
my $SPLITFH;

sub splitFileSuffix($) {
    return sprintf("%03i",$_[0]) . $compressedFastqSuffix;
}

sub updateSplitFH {
    if(defined $SPLITFH) {
        close($SPLITFH)
          or croak("ERROR: error in split file write\n");
    }
    $splitNumber++;
    my $cmd="| gzip -c --fast > ". $fastqPathPrefix . splitFileSuffix($splitNumber);
    open($SPLITFH,$cmd)
      or croak("ERROR: can't open process '$cmd'\n");
}

updateSplitFH();
my $count=0;

# Note fastq parse is setup to be relatively general (ie. do not
# require the '+' header-line to be empty and do not really on 2 w/s
# separated fields in the '@' header line). This means that we're
# obliged to transition between three states to parse the file
# correctly because both the '@' and '+' character can occur in the
# quality string. Thus we loop through neutral->read->qual states,
# taking advantage of two constraints: (1) No '+' can occur in the
# read and (2) the qual string must be the same length as the read
# string.
#
my ($is_read,$is_qual,$read_length,$qual_length) = (0,0,0,0); # use these to track where we are in the fastq
for my $fastq (@fastqs) {
    my $fastqPath = File::Spec->catfile($sampleTmpDir,$fastq);
    open(my $FASTQFH,"<$fastqPath") or croak("ERROR: can't open file: $fastqPath\n");
    while(<$FASTQFH>) {
        chomp;
        if     ($is_qual) {
            $qual_length += length($_);
            if     ($qual_length == $read_length) { # qual->neutral state transition
                $is_qual=0;
            } elsif($qual_length > $read_length) {
                croak("ERROR: Can't parse fastq file: '$fastq' $!\n");
            }
        } elsif(/^\+/){ #read->qual state transition
            $is_read=0;
            $is_qual=1;
            $qual_length=0;
        } elsif($is_read) {
            $read_length += length($_);
        } elsif(/^@/) { # neutral->read state transition
            $count++;
            if($count > $fastqClusterCount) {
                $count=1;
                updateSplitFH();
            }
            $is_read=1;
            $read_length=0;
        }
        print $SPLITFH $_ . "\n";
    }
    close($FASTQFH);
}
if(defined $SPLITFH) {
    close($SPLITFH)
      or croak("ERROR: error in split file write\n");
}

for(my $i=$splitNumber;$i;--$i) {
    my $suffix = splitFileSuffix($i);
    my $old=$fastqPathPrefix . $suffix;
    my $new=$finalPathPrefix . $suffix;
    move($old,$new)
      or croak("ERROR: failed to move '$old' to '$new'\n");
}

# remove input fastqs after completing split & move:
#
for my $fastq (@fastqs) {
    my $fastqPath = File::Spec->catfile($sampleTmpDir,$fastq);
    unlink($fastqPath);
}

1;

__END__

