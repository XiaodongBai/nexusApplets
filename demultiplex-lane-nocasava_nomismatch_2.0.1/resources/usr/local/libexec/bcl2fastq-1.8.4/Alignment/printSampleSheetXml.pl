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

printSampleSheetXml.pl

=cut

=head1 DESCRIPTION

Converts the given sample sheet into an Xml sample sheet

=cut

=head1 AUTHOR

Roman Petrovski

=head1 SYNOPSIS

    $ printSampleSheetXml.pl --sample-sheet <sample-sheet>

Converts the given sample sheet into an Xml sample sheet

=head1 DETAIL

All output is on STDOUT.

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
my $inputSampleSheet = undef;

GetOptions('help|?' => \$help, man => \$man, 'sample-sheet=s' => \$inputSampleSheet) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2,  -input => $1) if ($man and $0 =~ /(.*)/);
pod2usage("$0: Unexpected parameter(s): @ARGV.\n")  if (0 < @ARGV);

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl';

use Casava::Common::Log qw(initLog);
use Casava::Demultiplex::SampleSheet;

my $ss = Casava::Demultiplex::SampleSheet::create($inputSampleSheet);
open(my $sampleSheetFH,"<$inputSampleSheet") 
    or errorExit("ERROR: Failed to open input file '$inputSampleSheet' $!");
$ss->load($sampleSheetFH);
close($sampleSheetFH);

my $sx = Casava::Demultiplex::SampleSheet::Xml->new();
$sx->clone($ss);

open(my $ssXmlFH,">-") or errorExit("ERROR: Failed to open STDOUT for output $!");
$sx->save($ssXmlFH);
close($ssXmlFH);

1;

__END__
