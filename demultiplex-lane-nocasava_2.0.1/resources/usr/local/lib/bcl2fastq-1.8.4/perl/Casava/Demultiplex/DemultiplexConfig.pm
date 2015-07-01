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
  $flowCellId = $sampleSheet->flowCellId;
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

package Casava::Demultiplex::DemultiplexConfig;
use base 'Casava::Demultiplex::SampleSheet';

use strict;
use warnings "all";
use Exporter 'import';
use Carp;
use XML::Simple;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Common::Log qw(errorExit);


=pod

=head2 load

Load the sample sheet as XML data.

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

sub load($$)
{
    my $self = shift;
    my $handle = shift;
    my $xml = XMLin($handle, SuppressEmpty => 1, ForceArray => ['Lane', 'Sample'])
              or errorExit "ERROR: couldn't load sample sheet: $!";

#    die Data::Dumper->Dumper($xml);
    $self->software($xml->{Software});
    
    my $flowCellInfo = $xml->{FlowcellInfo};
    foreach (qw(Lane ID Operator Recipe Desc))
    {
        errorExit "ERROR: no '$_' in sample sheet FlowcellInfo." unless exists $flowCellInfo->{$_};
    }
    $self->flowCellId($flowCellInfo->{ID});
    $self->operator($flowCellInfo->{Operator});
    $self->recipe($flowCellInfo->{Recipe});
    $self->summary($flowCellInfo->{Desc});
    foreach my $lane (@{$flowCellInfo->{Lane}})
    {
        logWarning ("No 'Sample' element in <Lane Number=\"$lane->{Number}\"> of sample sheet")
             unless exists $lane->{Sample};
        foreach my $sample (@{$lane->{Sample}})
        {
            errorExit "ERROR: no 'SampleId' property found in Sample element\n"
                  unless exists $sample->{SampleId};
            $self->sample($lane->{Number}, $sample->{Index},
                          {SAMPLE_ID => $sample->{SampleId},
                           SPECIES => $sample->{Ref},
                           DESCRIPTION => $sample->{Desc},
                           CONTROL => $sample->{Control},
                           PROJECT_ID => $sample->{ProjectId}});
        }
    }
    
    return $self;
}

sub cmdAndArgs($;$)
{
    my $self = shift;
    if (@_) { $self->{Software}->{CmdAndArgs} = shift }
    return $self->{Software}->{CmdAndArgs};
}

sub baseCallsSoftware($;$)
{
    my $self = shift;
    if (@_) { $self->{Software}->{Software} = shift }
    return $self->{Software}->{Software};
}

sub intensitiesSoftware($;$)
{
    my $self = shift;
    if (@_) { $self->{Software}->{Software}->{Software} = shift }
    return $self->{Software}->{Software}->{Software};
}

sub controlSoftware($;$)
{
    my $self = shift;
    if (@_) { $self->{Software}->{Software}->{Software}->{Software} = shift }
    return $self->{Software}->{Software}->{Software}->{Software};
}

sub software($;$)
{
    my $self = shift;
    if (@_) { $self->{Software} = shift }
    return $self->{Software};
}

=pod

=head2 save

Dump sample sheet as xml data.

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

    my $xml = {};
    $xml->{Software} = $self->software;

    my $flowCellInfo = {};
    $flowCellInfo->{ID} = $self->flowCellId;
    $flowCellInfo->{Recipe} = $self->recipe;
    $flowCellInfo->{Operator} = $self->operator;
    $flowCellInfo->{Desc} = $self->summary;
    for my $lane (sort $self->laneNumbers)
    {
        my @samples = ();
        for my $barcode (sort $self->barcodes($lane))
        {
            my $sample = {Ref => $self->species($lane,$barcode),
                          SampleId => $self->sampleId($lane,$barcode),
                          Index => $barcode, 
                          Desc => $self->description($lane,$barcode),
                          Control => $self->control($lane,$barcode)};
            $sample->{ProjectId} = $self->projectId($lane,$barcode)
                        if defined $self->projectId($lane,$barcode);
            push @samples, $sample;
        }
        push @{$flowCellInfo->{Lane}}, {Number => $lane,
                               Sample => [@samples] };
    }
    $xml->{FlowcellInfo} = $flowCellInfo;
    print $handle XMLout($xml, KeyAttr=>[], RootName => "DemultiplexConfig", NoSort => 1);
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


