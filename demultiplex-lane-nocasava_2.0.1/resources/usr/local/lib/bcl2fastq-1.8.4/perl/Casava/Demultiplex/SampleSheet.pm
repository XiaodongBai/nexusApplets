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
  my $sampleSheet = Casava::SampleSheet::create('SampleSheet.csv'); # Create Csv instance
  open my $handle, '<', 'SampleSheet.csv' or die "Couldn't open SampleSheet.csv: $!";
  $sampleSheet->load($handle) or die "Couldn't load SampleSheet.csv");
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
  Project

The following optional column is added for the directory mapping after 
demultiplexing:

  Directory

It is assumed that there is only one flowcell id for the whole sample sheet.
Each lane can have several barcodes. Each line in the sample sheet is uniquely
identified by the couple (lane id, barcode). There are no constraints on the
other columns of the file.

=head1 SUBROUTINES

=cut

# POD for for each subroutines is right before the implementation
# More POD after __END__

package Casava::Demultiplex::SampleSheet;

use strict;
use warnings "all";
use Exporter 'import';
use Carp;

use Data::Dumper;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Demultiplex::SampleSheet::Csv;
use Casava::Demultiplex::SampleSheet::Make;
use Casava::Demultiplex::SampleSheet::Xml;

our @EXPORT_OK = qw(&new &create @FILE_FORMATS $barcodeComponentSeparator);

our @FILE_FORMATS = qw(csv xml mk);

our $barcodeComponentSeparator = '-';

=pod

=head2 new

Create a new empty samplesheet, using the derived class

I<Parameters:>

Filename (to figure out type from extension)

I<Returns:>

Nothing

I<Exceptions:>

None

=cut

sub create($)
{
    my ($filename) = @_;
    return Casava::Demultiplex::SampleSheet::Xml->new()  if $filename =~ /\.xml$/i;
    return Casava::Demultiplex::SampleSheet::Csv->new()  if $filename =~ /\.csv$/i;
    return Casava::Demultiplex::SampleSheet::Make->new() if $filename =~ /\.mk$/i;
}

=pod

=head2 new

Create a new empty samplesheet

I<Parameters:>

None

I<Returns:>

Nothing

I<Exceptions:>

None

=cut

sub new(;$)
{
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self  = {};
    $self->{FLOW_CELL_ID} = undef;
    $self->{RECIPE}       = undef;
    $self->{OPERATOR}     = undef;
    $self->{SUMMARY}      = undef;
    $self->{LANES}        = {};
    bless ($self, $class);
    return $self;
}

sub dump($;$)
{
    my $self = shift;
    my $lane = undef;
    if (@_)
    {
        $lane = shift;
        return undef unless exists $self->{LANES}->{$lane}; 
    }
    my $data = Data::Dumper->new( ($lane) ? [$self->{LANES}->{$lane}] : [$self->{LANES}] );
    $data->Purity(1); 
    return $data->{todump}->[0];
}

# filter is used to specify a particular project/sample combination to
# be included, while everything else is skipped over:
#
sub clone($$;\%)
{
    my $self         = shift;
    my $sampleSheet  = shift;
    my $filter      = shift;
    $filter = {} if !defined $filter;
    $self->flowCellId($sampleSheet->flowCellId);
    $self->recipe($sampleSheet->recipe);
    $self->operator($sampleSheet->operator);
    $self->summary($sampleSheet->summary);

    foreach my $lane (sort $sampleSheet->laneNumbers)
    {
        foreach my $barcode ($sampleSheet->barcodes($lane))
        {
            my $skip = 0;
            foreach my $fkey (keys %$filter)
            {
                if ($filter->{$fkey} ne $sampleSheet->$fkey($lane,$barcode))
                {
                    $skip = 1;
                    last;
                }
            }
            $self->sample($lane,$barcode,{ %{$sampleSheet->sample($lane,$barcode)} }) unless $skip;
        }
    }
}

sub flowCellId($;$)
{
    my $self = shift;
    if (@_) { $self->{FLOW_CELL_ID} = shift }
    return $self->{FLOW_CELL_ID};
}

sub recipe($;$)
{
    my $self = shift;
    if (@_) { $self->{RECIPE} = shift }
    return $self->{RECIPE};
}

sub operator($;$)
{
    my $self = shift;
    if (@_) { $self->{OPERATOR} = shift }
    return $self->{OPERATOR};
}

sub summary($;$)
{
    my $self = shift;
    if (@_) { $self->{SUMMARY} = shift }
    return $self->{SUMMARY};
}

=pod

=head2 laneNumbers

Return the list of lane numbers matching a given filtering criteria.

I<Parameters:>

=over 4

=item self

=item filter (optional)

a {key, value}, where the key is a field of the sample (targetted use case is the 'SAMPLE_ID').

=back

I<Returns:>

The list of lanes in the sample sheet that match the filter (all lanes
if the filter is not present). Note that depending on the filter and
on the content of the sample sheet, some lanes might appear several
times.

I<Exceptions:>

None

=cut
sub laneNumbers($;\%)
{
    my $self  = shift;
    my @lanes = keys %{$self->{LANES}};
    return @lanes  unless @_;
    my $filter         = shift; 
    my ($fkey,$fvalue) = each %$filter;
    my @filtered = ();
    foreach my $lane (@lanes)
    {
        foreach my $barcode ($self->barcodes($lane))
        {
            push @filtered, $lane
                 if !defined $fkey || !defined $fvalue 
                 || $fvalue eq $self->$fkey($lane,$barcode);
        }
    }
    return @filtered;
}

sub lanes($)
{
    my $self = shift;
    return sort keys %{$self->{LANES}};
}

sub barcodes($$)
{
    my $self = shift;
    my $laneNumber = shift;
    croak "Casava::SampleSheet::barcodes: Unspecified lane"  unless $laneNumber;
    return undef unless exists $self->{LANES}->{$laneNumber};
    return keys %{$self->{LANES}->{$laneNumber}};
}

sub isDemux($)
{
    my $self = shift;
    for my $lane ($self->lanes) {
        return 1 if($self->isLaneDemux($lane));
    }
    return 0;
}

sub isLaneDemux($$)
{
    my $self = shift;
    my $lane = shift;
    for my $barcode ($self->barcodes($lane)) {
        return 1 if($barcode ne '');
    }
    return 0;
}

sub sample($$$;$)
{
    my $self = shift;
    my $laneNumber = shift;
    my $barcode = shift;
    if (@_) { $self->{LANES}->{$laneNumber}->{$barcode} = shift }
    return undef unless exists $self->{LANES}->{$laneNumber};
    return $self->{LANES}->{$laneNumber}->{$barcode};
}

sub sampleId($$$;$)
{
    my $self = shift;
    my $laneNumber = shift;
    my $barcode = shift;
    if (@_) { $self->{LANES}->{$laneNumber}->{$barcode}->{SAMPLE_ID} = shift }
    return undef unless exists $self->{LANES}->{$laneNumber};
    return undef unless exists $self->{LANES}->{$laneNumber}->{$barcode};
    return $self->{LANES}->{$laneNumber}->{$barcode}->{SAMPLE_ID};
}

sub species($$$;$)
{
    my $self = shift;
    my $laneNumber = shift;
    my $barcode = shift;
    if (@_) { $self->{LANES}->{$laneNumber}->{$barcode}->{SPECIES} = shift }
    return undef unless exists $self->{LANES}->{$laneNumber};
    return undef unless exists $self->{LANES}->{$laneNumber}->{$barcode};
    return $self->{LANES}->{$laneNumber}->{$barcode}->{SPECIES};
}

sub description($$$;$)
{
    my $self = shift;
    my $laneNumber = shift;
    my $barcode = shift;
    if (@_) { $self->{LANES}->{$laneNumber}->{$barcode}->{DESCRIPTION} = shift }
    return undef unless exists $self->{LANES}->{$laneNumber};
    return undef unless exists $self->{LANES}->{$laneNumber}->{$barcode};
    return $self->{LANES}->{$laneNumber}->{$barcode}->{DESCRIPTION};
}

sub control($$$;$)
{
    my $self = shift;
    my $laneNumber = shift;
    my $barcode = shift;
    if (@_) { $self->{LANES}->{$laneNumber}->{$barcode}->{CONTROL} = shift }
    return undef unless exists $self->{LANES}->{$laneNumber};
    return undef unless exists $self->{LANES}->{$laneNumber}->{$barcode};
    return $self->{LANES}->{$laneNumber}->{$barcode}->{CONTROL};
}

sub projectId($$$;$)
{
    my $self = shift;
    my $laneNumber = shift;
    my $barcode = shift;
    if (@_) { $self->{LANES}->{$laneNumber}->{$barcode}->{PROJECT_ID} = shift }
    return undef unless exists $self->{LANES}->{$laneNumber};
    return undef unless exists $self->{LANES}->{$laneNumber}->{$barcode};
    return $self->{LANES}->{$laneNumber}->{$barcode}->{PROJECT_ID};
}

sub projectDirName($$$)
{
    my $self = shift;
    my $laneNumber = shift;
    my $barcode = shift;
    return undef unless exists $self->{LANES}->{$laneNumber};
    return undef unless exists $self->{LANES}->{$laneNumber}->{$barcode};
    my $projectId = $self->projectId($laneNumber,$barcode);
    return "Project_$projectId" unless 'Undetermined_indices' eq $projectId;
    return $projectId;
}

sub sampleDirName($$$)
{
    my $self = shift;
    my $laneNumber = shift;
    my $barcode = shift;
    return undef unless exists $self->{LANES}->{$laneNumber};
    return undef unless exists $self->{LANES}->{$laneNumber}->{$barcode};
    return "Sample_" . $self->sampleId($laneNumber,$barcode);
}

sub projectSampleDirName($$$)
{
    my $self = shift;
    my $laneNumber = shift;
    my $barcode = shift;
    return undef unless exists $self->{LANES}->{$laneNumber};
    return undef unless exists $self->{LANES}->{$laneNumber}->{$barcode};
    return undef unless length($self->sampleId($laneNumber,$barcode));
    return File::Spec->catdir(projectDirName($self, $laneNumber, $barcode), sampleDirName($self, $laneNumber, $barcode));
}


sub forall($$;\@)
{
    my $self = shift;
    my $function = shift;
    foreach my $lane (sort $self->laneNumbers)
    {
        foreach my $barcode ($self->barcodes($lane))
        {
            $self->$function($lane,$barcode,@_);
        }
    }
    return $self;
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


