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

This file is part of the Consensus Assessment of Sequence And VAriation
(CASAVA) software package.

=cut

=head1 NAME

configureDataset.pl

Create the analysis-specific configuration file for each dataset.

=cut

=head1 DESCRIPTION

Uses the content of the config.xml to generate the initialization of all the
Makefile variables required for the analysis of the dataset.

=cut

=head1 AUTHOR

Come Raczy

=head1 SYNOPSIS

    $ configureDataset.pl --project <project> --sample <sample> --lane <lane> --barcode <barcode> <config.xml>

Runs configureDataset.pl on <config.xml> and generate the resulting
Makefile configuration on the standard output.

=head1 DETAIL

The top level operations are:

* parse the config file
* generate the list of variables
* assign a value to each variable
* produce the makefile assignment

Example output:<<EOF;

dataset := $(project)_$(sample)_$(lane)_$(barcode)

$(dataset)_ANALYSIS := none

EOF

=cut

use warnings "all";
use strict;
use Getopt::Long;
use Pod::Usage;
my $ENV_REF = \%ENV; # protect from CMake substitution
$ENV_REF->{"PATH"} = $1 if $ENV_REF->{"PATH"} =~ /^(.*)$/; # untaint the path

my $man = 0;
my $help = 0;
my $project;
my $sample;
my $lane;
my $barcode;
my $reference;
my $configFile = "config.xml";

GetOptions('help|?' => \$help,
           man => \$man,
           'project=s' => \$project,
           'sample=s' => \$sample,
           'lane=s' => \$lane,
           'barcode=s' => \$barcode,
           'reference=s' => \$reference) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2,  -input => $1) if ($man and $0 =~ /(.*)/);
pod2usage("$0: The path to the configuration file is mandatory.\n")  if (1 != @ARGV);
$configFile = $ARGV[0] if (1 == @ARGV);

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl';

use Casava::Alignment;
use Casava::Common::Log qw(initLog);

initLog( undef, 3, 0, 0, 0);
Casava::Alignment::configureDataset($configFile, $project, $sample, $lane, $barcode, $reference);

1;

__END__
