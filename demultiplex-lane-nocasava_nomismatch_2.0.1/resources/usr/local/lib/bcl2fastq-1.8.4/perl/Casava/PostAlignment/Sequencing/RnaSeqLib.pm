package Casava::PostAlignment::Sequencing::RnaSeqLib;

# PROJECT: CASAVA
# MODULE:  $RCSfile: RnaSeqLib.pm,v $
# AUTHOR:  Durinck Steffen, Lukasz Szajkowski, Richard Carter
#
# Copyright (c) 2008, 2009 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#
# The library contains procedures and variables specific to RNA-Seq analysis.
#

=pod

=head1 NAME

Casava::PostAlignment::Sequencing::RnaSeqLib.pm - Perl utility library for running BullForg 
    analysis.

=head1 SYNOPSIS

The library contains procedures and variables specific to CASAVA
application. 
 
use Casava::PostAlignment::Sequencing::RnaSeqLib.pm qw();  

=head1 DESCRIPTION
    
# Global variable


=head1 AUTHORSHIP

Copyright (c) 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

Created by Durinck, Steffen <sdurinck@illumina.com>,
Lukasz Szajkowski <lszajkowski@illumina.com>

=cut

BEGIN {
	use Exporter();
	@ISA       = qw(Exporter);
	@EXPORT    = qw();
	@EXPORT_OK = qw(&getGenesFromExons &countSpliceJunctions $tmpSpliceCountFileName);
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);
use Carp;

use Casava::Common::Log;
use Casava::PostAlignment::Sequencing::GenomicIOLib qw( %exportFields parseSpliceName);
sub getGenesFromExons(\%\%);
sub countSpliceJunctions($$\%\%);

# filename for temp storage of splice junction counts:
#
our $tmpSpliceCountFileName = "tmpSpliceJunctionCounts.txt";

=pod

=head1 Calculates the number of reads covering each splice site

=over 4

=item countSpliceJunctions($inFileHandleRef, $chrom, $srasMin, $Position, $junstionStorage)

The procedure calculates the number of read starting at any given postion.

Parameters:
    $inFileHandleRef - handle to reads input stream
    $srasMin         - Reads with alignment score below that will be skipped   
    $Position        - Hash that will receive number of reads starting at a particular position
    $junctionStorage - Hash that will receive information about splice junctions
Returns:
    (total number of reads, number of used reads, total number of bases, number of used bases)
    
=back

=cut

sub countSpliceJunctions($$\%\%) {

	my ( $inFileHandleRef, $srasMin, $Position, $junctionStorage ) = @_;

	# Find the input directories
	my ($lineCount, $usedReads, $usedBases) = (0, 0, 0);

   	while (<$inFileHandleRef>) {
        chomp;
		++$lineCount;
		my $line = $_;
		my @read = split /\t/, $line;
		
		next unless ( $read[ $exportFields{SingleScore} ] >= $srasMin );
		
		my $featureName = $read[ $exportFields{Contig} ];
		my $readLength  = length( $read[ $exportFields{Seq} ] );
		errorExit "ERROR: Contig is empty in: $line" if ( !defined $featureName );
		my $spliceJunctionFeatureRef = $junctionStorage->{feature}->{$featureName};
		my $new = 0;
		if ( !$spliceJunctionFeatureRef ) {
			my %spliceJunctionFeature = ();
			$new = 1;
			$spliceJunctionFeatureRef = \%spliceJunctionFeature;
		}    # if
		parseSpliceName( $featureName, %{$spliceJunctionFeatureRef} );
		my $start     = $spliceJunctionFeatureRef->{start};
		my $chromName = $spliceJunctionFeatureRef->{scaffold};

		$junctionStorage->{feature}->{$featureName} = $spliceJunctionFeatureRef;
		$junctionStorage->{readLength} = $readLength;
		if ( $new ) {
			push @{ $junctionStorage->{positions}->{$start} }, $spliceJunctionFeatureRef;
		}    # if
		++$junctionStorage->{count};

        ++$spliceJunctionFeatureRef->{coverage};
            
        ++$Position->{$start};
        ++$usedReads;
        $usedBases += $readLength;

	}

    return ($lineCount, $usedReads, $usedBases, $usedBases);
}

=pod

=head1 Extract list of genes from list of exons

=over 4

=item getGenesFromExons($exonfeatureStorage, $genesfeatureStorage)

Parameters:
    $exonfeatureStorage    - HASH MAP REF to exon features read with readFeatures procedure 
    $genesfeatureStorage   - HASH MAP REF to genes features read with readFeatures procedure 
Returns:
    nothing
    
=back

=cut

sub getGenesFromExons(\%\%) {
	my ( $exonfeatureStorage, $genesfeatureStorage ) = @_;
	my @featuresIds = keys %{ $exonfeatureStorage->{feature} };
	my @featuresSorted = sort {
		$exonfeatureStorage->{feature}->{$a}->{start} <=> $exonfeatureStorage->{feature}->{$b}->{start}
	} @featuresIds;

	foreach my $featureId (@featuresSorted) {
		my @genes = @{ $exonfeatureStorage->{feature}->{$featureId}->{genes} };
		foreach my $geneId (@genes) {

			$genesfeatureStorage->{feature}->{$geneId}->{size} +=
			  $exonfeatureStorage->{feature}->{$featureId}->{end} -
			  $exonfeatureStorage->{feature}->{$featureId}->{start} + 1;

            $genesfeatureStorage->{feature}->{$geneId}->{coverage} +=
              $exonfeatureStorage->{feature}->{$featureId}->{coverage} 
              if defined $exonfeatureStorage->{feature}->{$featureId}->{coverage};
			
			$genesfeatureStorage->{feature}->{$geneId}->{name} = $geneId;

			push @{ $genesfeatureStorage->{feature}->{$geneId}->{exons} }, $featureId;

			if ( !defined $genesfeatureStorage->{feature}->{$geneId}->{start} )
			{
				$genesfeatureStorage->{feature}->{$geneId}->{start} =
				  $exonfeatureStorage->{feature}->{$featureId}->{start};
			}    # if
			else {
				if ( $genesfeatureStorage->{feature}->{$geneId}->{start} >
					$exonfeatureStorage->{feature}->{$featureId}->{start} )
				{
					$genesfeatureStorage->{feature}->{$geneId}->{start} =
					  $exonfeatureStorage->{feature}->{$featureId}->{start};
				}    # if
			}    # else

			if ( !defined $genesfeatureStorage->{feature}->{$geneId}->{end} ) {
				$genesfeatureStorage->{feature}->{$geneId}->{end} =
				  $exonfeatureStorage->{feature}->{$featureId}->{end};
			}    # if
			else {
				if ( $genesfeatureStorage->{feature}->{$geneId}->{end} <
					$exonfeatureStorage->{feature}->{$featureId}->{end} )
				{
					$genesfeatureStorage->{feature}->{$geneId}->{end} =
					  $exonfeatureStorage->{feature}->{$featureId}->{end};
				}    # if
			}
		}
	}
}

1;           # says use was ok
__END__

