=head1 LICENSE

Copyright (c) 2007-2009 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And VAriation
(CASAVA) software package.

=head1 NAME

Casava::SampleSheet - Utility library for the management of sample sheets.

=head1 SYNOPSIS

  use Casava::Demultiplex::SampleSheet;
  $sampleSheet = Casava::SampleSheet->new()
  $sampleSheet->loadCsv("SampleSheet.csv") 
      or die "Couldn't load SampleSheet.csv");
  $recipe = $sampleSheet->recipe($lane, $barcode);
  $operator = $sampleSheet->operator($lane, $barcode);
  foreach $lane ($sampleSheet->laneNumbers)
  {
    foreach $barcode ($sampleSheet->barcodes($lane))
    {
        $sample = $sampleSheet->sample($lane, $barcode);
        $species = $sampleSheet->species($lane, $barcode);
        $description = $sampleSheet->description($lane, $barcode);
        $control = $sampleSheet->control($lane, $barcode);
    }
  }

=head1 DESCRIPTION

The sample sheet is used to describe the list of samples, and their associated
barcodes, available in each lane of a multiplexed run.

The default format for the sample sheet is a csv file, with the following
columns:

  Flow Cell Id
  Lane Id
  Sample Id
  Species
  Barcode
  Description
  Control
  Recipe
  Operator

The following optional column is added for the directory mapping after 
demultiplexing:

  Directory

It is assumed that there is only one flowcell id for the whole sample sheet.
Each lane can have several barcodes. Each line in the sample sheed is uniquely
identified by the couple (lane id, barcode). There are no constraints on the
other columns of the file.

=head1 SUBROUTINES

=cut

# POD for for each subroutines is right before the implementation
# More POD after __END__

package Casava::Demultiplex::SampleSheet::Csv;
use base 'Casava::Demultiplex::SampleSheet';

use strict;
use warnings "all";
use Exporter 'import';
use Carp;
use Text::ParseWords;

use Casava::Common::Utils;

=pod

=head2 load

Load the sample sheet as cvs data.

I<Parameters:>

=over 4

=item *

$handle

open handle to the actual data

=back

I<Returns:>

Nothing

I<Exceptions:>

=over 4

=item *

Incorrect number of columns

=back

=cut

sub load($$;$)
{
    my $self = shift;
    my $handle = shift;
    my $generateUndeterminedSamples = shift;
    $generateUndeterminedSamples = 1 if !defined $generateUndeterminedSamples;
    my $eol = getEOL($handle);
    local $/ = $eol; # tell chomp the EOL from .csv 

    my $header = <$handle>; # discard the first line
    my @allLanes;
    my %laneLastBarcode;
    while (<$handle>)
    {
        next if $_ =~ /^#|^\s*$/; 
        chomp;
        ##$_ =~ s/[\n\r]+//g; RP: as eol is supposed to be correct now, chomp should do the proper job.
        my @words = map {defined $_ ? $_ : ""} quotewords( ",", 1, $_);
        croak "ERROR: Wrong number of fields in sample sheet (expected: 10, got ".scalar(@words).": ".join(',', @words).")\n"
              unless scalar(@words) == 10;
        $self->flowCellId($words[0]) unless defined $self->flowCellId();
        croak "ERROR: FlowCell ID is inconsistent across Sample Sheet lines. Expected: '$self->flowCellId()', got $words[0]\n"
              unless $self->flowCellId() eq $words[0];

        $self->recipe($words[7]);
        $self->operator($words[8]);
        my $laneNumber = $words[1];
        my $barcode = $words[4];
        $barcode = 'NoIndex' if $barcode =~ /^\s*$/;
        my $projectId = ($words[9] ? $words[9] : ( $words[0] ? $words[0] : 'default' ) );
        my $sampleId = $words[2];
        # empty sample id in csv means the lane-barcode must not produce the output
        my $sample = {SAMPLE_ID=>$sampleId,
                      SPECIES=>$words[3],
                      DESCRIPTION=>$words[5],
                      CONTROL=>$words[6],
                      PROJECT_ID=>$projectId};
        $self->sample($laneNumber, $barcode, $sample);
        push @allLanes, $laneNumber;
        croak "ERROR: Conflicting sample sheet definitions for lane lane $laneNumber. Sample sheet line: $.. Existing:"
            .$laneLastBarcode{$laneNumber} . " New: $barcode. Context: ".join(',', @words)."\n"
              if exists $laneLastBarcode{$laneNumber} && 
                 ('NoIndex' eq $laneLastBarcode{$laneNumber} || 'NoIndex' eq $barcode);
        $laneLastBarcode{$laneNumber} = $barcode;
    }
    
    if ($generateUndeterminedSamples)
    {
        for my $laneNumber (@allLanes)
        {
            if ('NoIndex' ne $laneLastBarcode{$laneNumber} && !defined $self->sample($laneNumber, 'Undetermined'))
            {
                my $sample = {SAMPLE_ID=>'lane'.$laneNumber,
                              SPECIES=>'unknown',
                              DESCRIPTION=>"Clusters with unmatched barcodes for lane $laneNumber",
                              CONTROL=>'N',
                              PROJECT_ID=>'Undetermined_indices'};
                $self->sample($laneNumber, 'Undetermined', $sample);
            }
        }
    }
    return $self;
}

=pod

=head2 save

Dump sample sheet as csv data.

I<Parameters:>

=over 4

=item *

$handle

open handle to the actual data

=back

I<Returns:>

Nothing

I<Exceptions:>

=over 4

=item *

Incorrect number of columns

=back

=cut

sub save($$)
{
    my $self = shift;
    my $handle = shift;

    my @header = ("FCID","Lane","SampleID","SampleRef","Index","Description","Control","Recipe","Operator","SampleProject");
    my ($aLane)   = $self->laneNumbers;
    if (defined $aLane)
    {
        my ($anIndex) = $self->barcodes($aLane);
    }
    print $handle (join ',',@header),"\n";

    my $estr='';

    for my $lane (sort $self->laneNumbers)
    {
        for my $barcode (sort $self->barcodes($lane))
        {
            my @words = ();
            push @words, ($self->flowCellId() or $estr);
            push @words, $lane;
            push @words, ($self->sampleId($lane,$barcode) or $estr);
            push @words, ($self->species($lane,$barcode) or $estr);
            push @words, 'NoIndex' eq $barcode ? '' : $barcode;
            push @words, ($self->description($lane,$barcode) or $estr);
            push @words, ($self->control($lane,$barcode) or $estr);
            push @words, ($self->recipe or $estr);
            push @words, ($self->operator or $estr);
            push @words, ($self->projectId($lane,$barcode) or $estr);
            print $handle (join ',',@words),"\n";
        }
    }
}


1;
__END__

=pod

=head1 DEPENDENCIES

=over 4

=item Standard perl modules

strict, warnings, Exporter, Carp, Text::ParseWords

=item External perl modules

XML::Simple

=item Casava perl modules

None

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Come Raczy

=cut


