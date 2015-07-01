package Casava::PostAlignment::Sequencing::GenomicIOLib;

# PROJECT: CASAVA
# MODULE:  $RCSfile: GenomicIOLib.pm,v $
# AUTHOR:  Lukasz Szajkowski, Richard Carter
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#
# The library contains procedures and variables usefull in operating on Illumina
# genomic data e.g.: export files. The library is trying to abstract access to
# files.

=pod

=head1 NAME

Casava::Sequencing::GenomicIOLib.pm - The library contains procedures and variables 
usefull in operating on Illumina genomic data.

=head2 SYNOPSIS

use Casava::Sequencing::GenomicIOLib.pm qw();  

=head2 AUTHORSHIP

Copyright (c) 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

Created by Lukasz Szajkowski <lszajkowski@illumina.com>

=head1 DESCRIPTION

=head2 Overview

The library contains procedures and variables usefull in operating on Illumina 
genomic data e.g.: export files. The library is trying to abstract access to 
files.

=head2 Exports:
    writeRead(\%;\@)
	configureReadStorage(\%)
	flushReadFiles(\%;$)
	getNextRead(\%)
	rewindTo(\%;$;$)
	getExperimentName(\@;\%)
	readFeaturesENST($;\%)
	readFeaturesENST2($;$;\%)	

	parseSpliceName($\%)
		    
Global variables:
	%sortFields
	%exportFields
	%positionFields

=head1 SOURCE CODE

The most current release of the Perl source code 
for this module is available through CVS at genie01 a.k.a. 10.44.0.81 
cvs co BullFrog

=cut

BEGIN {
    use Exporter();
    @ISA       = qw(Exporter);
    @EXPORT    = qw();
    @EXPORT_OK = qw(&writeRead2Scaffold &writeRead &configureReadStorage &flushReadFiles &readFasta
      &rewindTo &getNextRead &readFeaturesENST &readFeaturesUCSCB &readFeaturesENST2
      &printFeatureCounts &saveFeatureCounts &readFeatureCounts %positionFields
      &parseSpliceName &sortFields %exportFields %doubleExportFields %featureCounts &outputSEread
      &outputPEread %predictedMateStrand $readsIdxHeader %readsIdxFields @nmnmTags $VERSION);
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);

use IO::File;
use Carp;

use Casava::Common::Log;
use Casava::Common::IOLib qw(bufferedPrint);

use Casava::PostAlignment::Sequencing::Config qw(%chrEnds %CONF_PROJ);

sub writeRead2Scaffold(\%;$;\@);
sub writeRead(\%;\@);
sub configureReadStorage(\%);
sub flushReadFiles(\%;$);
sub getPathForReadType(\%;$;$;$);
sub getNextRead(\%);
sub rewindTo(\%;$;$);
sub getExperimentName(\@;\%);
sub readFeaturesENST($;\%);
sub readFasta($);
sub printFeatureCounts($$\%$$);
sub saveFeatureCounts($$\%$);
sub readFeatureCounts($;$;\%);
sub readFeaturesUCSCB($;$;\%);
sub readFeaturesENST2($;$;\%);
sub getGenesFromExons(\%;\%);
sub parseSpliceName($\%);
sub outputSEread(\%$$\@$);
sub outputPEread(\%$$$\@\@$);

## expot.txt format fields
our %exportFields = (
    MachineName   => 0,
    RunNumber     => 1,
    Lane          => 2,
    Tile          => 3,
    X             => 4,
    Y             => 5,
    Index         => 6,
    ReadNum       => 7,
    Seq           => 8,
    QualityString => 9,
    Chr           => 10,
    Contig        => 11,
    Posn          => 12,
    Strand        => 13,
    Descriptor    => 14,
    SingleScore   => 15,
    PairedScore   => 16,
    PartnerChr    => 17,
    PartnerContig => 18,
    PartnerOffset => 19,
    PartnerStrand => 20,
    Filter        => 21
);
## double expot.txt format fields
our %doubleExportFields = (
    MachineName        => 0,
    RunNumber          => 1,
    Lane               => 2,
    Tile               => 3,
    X                  => 4,
    Y                  => 5,
    Index              => 6,
    ReadNum            => 7,
    Seq                => 8,
    QualityString      => 9,
    Chr                => 10,
    Contig             => 11,
    Posn               => 12,
    Strand             => 13,
    Descriptor         => 14,
    SingleScore        => 15,
    PairedScore        => 16,
    PartnerChr         => 17,
    PartnerContig      => 18,
    PartnerOffset      => 19,
    PartnerStrand      => 20,
    Filter             => 21,
    read2MachineName   => 22,
    read2RunNumber     => 23,
    read2Lane          => 24,
    read2Tile          => 25,
    read2X             => 26,
    read2Y             => 27,
    read2Index         => 28,
    read2ReadNum       => 29,
    read2Seq           => 30,
    read2QualityString => 31,
    read2Chr           => 32,
    read2Contig        => 33,
    read2Posn          => 34,
    read2Strand        => 35,
    read2Descriptor    => 36,
    read2SingleScore   => 37,
    read2PairedScore   => 38,
    read2PartnerChr    => 39,
    read2PartnerContig => 40,
    read2PartnerOffset => 41,
    read2PartnerStrand => 42,
    read2Filter        => 43
);
## sort.txt format fields
our %sortFields = (
    ExpName            => 0,
    Seq                => 1,
    QualityString      => 2,
    Chr                => 3,
    Contig             => 4,
    Posn               => 5,
    Strand             => 6,
    Descriptor         => 7,
    SingleScore        => 8,
    PairedScore        => 9,
    PartnerChr         => 10,
    PartnerContig      => 11,
    PartnerOffset      => 12,
    PartnerStrand      => 13,
    read2ExpName       => 15,
    read2Seq           => 16,
    read2QualityString => 17,
    read2Chr           => 18,
    read2Contig        => 19,
    read2Posn          => 20,
    read2Strand        => 21,
    read2Descriptor    => 22,
    read2SingleScore   => 23,
    read2PairedScore   => 24,
    read2PartnerChr    => 25,
    read2PartnerContig => 26,
    read2PartnerOffset => 27,
    read2PartnerStrand => 28
);
## sort.count format fields
our %positionFields = (
    position => 0,
    nA       => 1,
    nC       => 2,
    nG       => 3,
    nT       => 4,
    num_read => 6,
    num_call => 7,
    score    => 8
);
## codingExons_unig format fields
our %featureFields = (
    startPos => 0,
    endPos   => 1,
    name     => 2
);
## refFlat exons format fields
our %refFlatExonsFields = (
    chrom => 0,
    start => 1,
    end   => 2,
    names => 3
);
our %defaultConfiguration = (
    isAppend        => 0,
    binSize         => 10000000,
    currentBinStart => -1,
    currentBinEnd   => -1,
    readType        => 'export',
    baseDir         => '',
    fields          => \%exportFields,
);
## codingExons  nonRedundantSet_Ensembl_49 format fields
our %featureFieldsENS2 = (
    gene_stable_id => 0,
    struct_gene_id => 1,
    exon_stable_id => 2,
    str_chrom_name => 3,
    start          => 4,
    end            => 5,
    strand         => 6,
    struct_biotype => 7
);
## _count.txt fields
our %featureCounts = (
    scaffold => 0,
    start    => 1,
    end      => 2,
    id       => 3,
    count    => 4,
    absCount => 5,
    score    => 6,
);

our %predictedMateStrand = (
    'F' => 'R',
    'R' => 'F'
);

our %readsIdxFields = (
    totalReads=>0,
    totalBases=>1,
    usedReads=>2,
    usedBases=>3,
    goodReads=>4,
    goodReadBases=>5,
    failedTileFilt=>6,
    failedFilt=>7,
    failedQC=>8,
    singleExclude=>9,
    nonUniqueAlign=>10,
    goodPairs=>11,
    mixedPairs=>12,
    riboCount=>13,
    mitoCount=>14,
    splicedReads=>15,
    exportFileDir=>16,
    exportSet=>17,
    date=>18,
    instrumentName=>19,
    runID=>20,
    laneNumber=>21,
    barcode=>22,
    highSDFragmentLength=>23,
    lowSDFragmentLength=>24,
    maxFragmentLength=>25,
    medianFragmentLength=>26,
    minFragmentLength=>27,
    nominalOrientation=>28,
    read1Length=>29,
    read2Length=>30
);

our $readsIdxHeader = '#\$ COLUMNS ' . join ("\t",
    sort {$readsIdxFields{$a} <=> $readsIdxFields{$b}} keys(%readsIdxFields));

our @nmnmTags = qw(TAG_noMatch TAG_qcFail TAG_nonUnique TAG_rm TAG_mixed TAG_control);

=pod

=head1 FUNCTIONS

=over 4

=cut

=pod

=item configureReadStorage($configRef)

The procedure configures the read storage.

B<Parameters:>

	$configRef         - configuration HASH MAP REF

B<Returns:> 

	HASH MAP REf to read storage	

=cut

sub configureReadStorage(\%) {
    croak "ERROR: configureReadStorage wrong parameters \n"
      unless ( 1 == @_ );
    my ( $configRef ) = @_;
    my $storageRef;
    foreach my $key ( keys %defaultConfiguration ) {
        if ( defined $configRef->{$key} ) {
            $storageRef->{config}->{$key} = $configRef->{$key};
        }    # if
        else {
            $storageRef->{config}->{$key} = $defaultConfiguration{$key};
        }    # else
    }    # foreach

    #my @fields = keys %{ $storageRef->{config}->{fields} };

    #my $strFi = join (',', @fields);
    #print "Init configureReadStorage $strFi\n";
    if ( !defined $storageRef->{config}->{fields}->{Posn} ) {
        errorExit "ERROR configureReadStorage Fields->Posn not defined\n";
    }    # if
    if ( !defined $storageRef->{config}->{fields}->{Chr} ) {
        errorExit "ERROR configureReadStorage Fields->Chr not defined\n";
    }    # if
    return $storageRef;
}

=pod

=item rewindTo($storageRef, $conting, $position)

The procedure moves the read pointer to chromosome/position

B<Parameters:>

    $storageRef  - storage HASH MAP REF
    $conting     - chromosome/scaffold or contig name
    $position    - position    

B<Returns:>
    read ARRAY  

=cut

sub rewindTo(\%;$;$) {
    croak "ERROR: rewindTo wrong parameters \n"
      unless ( @_ == 3 );
    my ( $storageRef, $conting, $position ) = @_;
    my $binId = sprintf "%04d", $position / $storageRef->{config}->{binSize};
    my $readType = $storageRef->{config}->{readType};
    if ( defined $storageRef->{currentFile} ) {
        my $file = $storageRef->{currentFile};
        close($file);
    }
    my $filePath =
      getPathForReadType( %{$storageRef}, $readType, $conting, $binId );
    $storageRef->{currentFile} = IO::File->new( '<' . $filePath )
      || errorExit "ERROR rewindTo() Couldn't "
      . "open file handle for $filePath $!\n";
    $storageRef->{currentPos}    = $position;
    $storageRef->{currentBin}    = $binId;
    $storageRef->{currentContig} = $conting;

    #my $binId = sprintf "%04d", $position / $storageRef->{config}->{binSize};
}

=pod

=item getNextRead($readRef)

The procedure takes next read from a file.

B<Parameters:>

    $storageRef  - storage HASH MAP REF

B<Returns:>
    read ARRAY  

=cut

sub getNextRead(\%) {
    croak "ERROR: getNextRead wrong parameters \n"
      unless ( @_ == 1 );
    my ($storageRef) = @_;
    my $file         = $storageRef->{currentFile};
    my $line         = <$file>;
    my @read;
    if ( !defined $line || $line eq '' ) {
        close($file);
        my $conting  = $storageRef->{currentContig};
        my $position = $storageRef->{currentPos};
        my $readType = $storageRef->{config}->{readType};
        my $binId    = sprintf "%04d", $storageRef->{currentBin} + 1;
        my $filePath =
          getPathForReadType( %{$storageRef}, $readType, $conting, $binId );
        ## skip empty bins
        while ( -z $filePath ) {
            my $currentBinStart =
              ( $binId + 1 ) * $storageRef->{config}->{binSize};
            $binId = sprintf "%04d",
              $currentBinStart / $storageRef->{config}->{binSize};
            $filePath =
              getPathForReadType( %{$storageRef}, $readType, $conting, $binId );
        }
        if ( -e $filePath ) {
            $storageRef->{currentFile} = IO::File->new( '<' . $filePath )
              || errorExit "ERROR getNextRead() Couldn't "
              . "open file handle for $filePath $!\n";
            $file                     = $storageRef->{currentFile};
            $line                     = <$file>;
            $storageRef->{currentBin} = $binId;
        }    # if
        else {
            ## Scaffold end - there is no more bins
            return 0;
        }    # else
    }    # if
    if ( !defined $line ) {
        errorExit "ERROR: getNextRead() inknown problem\n";
    }
    @read = split "\t", $line;
    $storageRef->{currentPos}++;
    return \@read;
}

=pod

=item writeRead2Scaffold($storageRef, $scaffoldFile, $readRef)

The procedure writes read to a buffred file using the 
supplied scaffold file name instead of the chromosome name
in the $readRef

B<Parameters:>

    $storageRef    - storage HASH MAP REF
    $scaffoldFile  - scaffold name to use
    $readRef       - read ARRAY REF

B<Returns:>
    status -0 ok; -1 error  

=cut

sub writeRead2Scaffold(\%;$;\@) {
    croak "ERROR: writeRead2Scaffold wrong parameters \n"
      unless ( @_ == 3 );
    my ( $storageRef, $scaffoldFile, $readRef ) = @_;

    my $position     = $readRef->[ $storageRef->{config}->{fields}->{Posn} ];
    my $scaffold     = $scaffoldFile;

    my $readType = $storageRef->{config}->{readType};
    my $binId    = sprintf "%04d", $position / $storageRef->{config}->{binSize};

    my $fileHashKey = $scaffold . $binId . $readType;

    my $file = $storageRef->{files}->{$fileHashKey};
    if ( !defined $file )
    {
        my $filePath =
          getPathForReadType( %{$storageRef}, $readType, $scaffold, $binId );

        my $operation = $storageRef->{config}->{isAppend} ? '>>' : '>';

        $file = IO::File->new( $operation . $filePath )
              || errorExit "$0::ERROR writeRead2Scaffold Couldn't create/open file handle for $filePath $!\n";
        $storageRef->{files}->{$fileHashKey} = $file;
    }    # if

    print $file ( join "\t", @{$readRef} ) . "\n";
}

=pod

=item writeRead($storageRef, $readRef)

The procedure writes read to a buffered file.

B<Parameters:>

    $storageRef  - storage HASH MAP REF
    $readRef     - read ARRAY REF

B<Returns:>
    status -0 ok; -1 error  

=cut

sub writeRead(\%;\@) {
    croak "ERROR: writeRead wrong parameters \n"
      unless ( @_ == 2 );
    my ( $storageRef, $readRef ) = @_;

    my $chrom = $readRef->[ $storageRef->{config}->{fields}->{Chr} ];
    return writeRead2Scaffold(%$storageRef, $chrom, @$readRef);
}

=pod

=item flushReadFiles($storageRef, $isClose)

The procedure writes flushes open storage files.

B<Parameters:>

    $storageRef  - storage HASH MAP REF
    $isClose     - 1 to close flushed files.

B<Returns:>
    nothing

=cut

sub flushReadFiles(\%;$) {
    croak "ERROR: flushReadFiles wrong parameters \n"
      unless ( @_ == 2 );
    my ( $storageRef, $isClose ) = @_;
    foreach my $fileKey (keys %{$storageRef->{files}})
    {
        my $file = $storageRef->{files}->{$fileKey};
        $file->flush();
        if ($isClose) {
            close($file);
            delete $storageRef->{files}->{$fileKey};
        }
    }
}

=pod

=item getPathForReadType(\%;$;$;$)

The procedure maps file type to file path.

B<Parameters:>

    $storageRef  - storage HASH MAP REF
    $readType    - read type

B<Returns:>
    Path to file

=cut

sub getPathForReadType(\%;$;$;$) {
    croak "ERROR: getPathForReadType wrong parameters \n"
      unless ( @_ == 4 );
    my ( $storageRef, $readType, $chrom, $binId ) = @_;
    my $fileName;
    my $baseDir = $storageRef->{config}->{baseDir};
    my $path = File::Spec->catdir( $baseDir, $chrom, $binId );
    if ( $readType eq 'sort' ) {
        $fileName    = 'export.txt';
    }
    if ( $readType eq 'sorted' ) {
        $fileName    = 'sorted.txt';
    }
    elsif ( $readType eq 'rmDup' ) {
        $fileName    = 'sort_export.txt';
    }
    elsif ( $readType eq 'rmDup_orph' ) {
        $fileName    = 'sort_orph_export.txt';
    }
    elsif ( $readType eq 'rmDup_anom' ) {
        $fileName    = 'sort_anom_export.txt';
    }
    elsif ( $readType eq 'corect_sort' ) {
        $fileName    = 'corect_sort_' . $binId . '_export.txt';
        $path        = $baseDir;
    }
    elsif ( $readType eq 'corect_sort_anom' ) {
        $fileName    = 'corect_sort_anom_' . $binId . '_export.txt';
        $path        = $baseDir;
    }
    elsif ( $readType eq 'corect_sort_orph' ) {
        $fileName    = 'corect_sort_orph_' . $binId . '_export.txt';
        $path        = $baseDir;
    }
    elsif ( $readType eq 'unsort.wrong' ) {
        $fileName    = 'wrong_export.txt';
    }
    elsif ( $readType eq 'anom.unsort.wrong' ) {
        $fileName    = 'anom_wrong_export.txt';
    }
    elsif ( $readType eq 'sort.wrong' ) {
        $fileName    = 'sort_wrong_export.txt';
    }
    elsif ( $readType eq 'unsort_single' ) {
        $fileName    = 'single_export.txt';
    }
    else {
        errorExit "ERROR: getPathForReadType() unsuported read type [$readType].\n";
    }
    my $filePath = File::Spec->catfile( $path, $fileName );
    return $filePath;
}

=pod

=item getExperimentName($readRef)

The procedure generates experiment name (unique read id).

B<Parameters:>

    $readRef   - read ARRAY REF
    $fieldsRef - HASH MAP REF with fields

B<Returns:>
    experiment name (read id)

=cut

sub getExperimentName(\@;\%) {
    croak "ERROR: getExperimentName wrong parameters \n"
      unless ( @_ == 2 );
    my ( $readRef, $fieldsRef ) = @_;
    my $expName = @{$readRef}[ $fieldsRef->{MachineName} ];
    $expName .= "_" . @{$readRef}[ $fieldsRef->{RunNumber} ];
    $expName .= ":" . @{$readRef}[ $fieldsRef->{Lane} ];
    $expName .= ":" . @{$readRef}[ $fieldsRef->{Tile} ];
    $expName .= ":" . @{$readRef}[ $fieldsRef->{X} ];
    $expName .= ":" . @{$readRef}[ $fieldsRef->{Y} ];
    $expName .= ":" . @{$readRef}[ $fieldsRef->{Index} ];
    $expName .= ":-" . @{$readRef}[ $fieldsRef->{ReadNum} ];
    return $expName;
}

=pod

=item readFeatureCoverage($featureStorage)

The procedure reads features and their coverage from a file printed by readFeatureCoverage

B<Parameters:>

    $featurefileName  - feature count file (e.g.: c22_genes_count.txt)
    $fileName         - feature count file
    $chrom            - chromosome name
    $featureStorage   - HASH MAP REF to features read with readFeatures procedure

B<Returns:>
    nothing

=cut

sub readFeatureCounts($;$;\%) {
    croak "ERROR: readFeatureCounts wrong parameters \n"
      unless ( @_ == 3 );
    my ( $featurefileName, $chrom, $featureStorage ) = @_;
    my $featurefile = IO::File->new($featurefileName)
      || errorExit
"ERROR: readFeatureCounts() Couldn't create/open file handle for $featurefileName $!\n";
    while (<$featurefile>) {
        chomp;
        next if (/^#/);
        my $featureRef;
        my @row       = split "\t", $_;
        my $featureId = $row[ $featureCounts{id} ];

        #if ( !defined $featureStorage->{feature}->{$featureId} )
        #{
        my %feature = ();

        #}
        $featureStorage->{feature}->{$featureId} = \%feature;
        $featureStorage->{feature}->{$featureId}->{start} =
          $row[ $featureCounts{start} ];
        $featureStorage->{feature}->{$featureId}->{end} =
          $row[ $featureCounts{end} ];
        $featureStorage->{feature}->{$featureId}->{name} =
          $row[ $featureCounts{id} ];
        $featureStorage->{feature}->{$featureId}->{id} =
          $row[ $featureCounts{id} ];
        push @{ $featureStorage->{positions}
              ->{ $featureStorage->{feature}->{$featureId}->{start} } },
          \%feature;
        $featureStorage->{count}++;
        $featureStorage->{feature}->{$featureId}->{normCoverage} =
          $row[ $featureCounts{count} ];
        $featureStorage->{feature}->{$featureId}->{coverage} =
          $row[ $featureCounts{absCount} ];

        if ( !defined $row[ $featureCounts{score} ] ) {
            $featureStorage->{feature}->{$featureId}->{score} =
              $row[ $featureCounts{score} ];
        }
    }    # while
    close($featurefile);
}

=pod

=item printFeatureCounts

The procedure prints feature counts (number of bases or reads starts in a feature) into an open handle
e.g.:
c22 15977189    15982257    ENSG00000183307 0.52    2648
c22 15998411    16026177    ENSG00000069998 0.78    21791

B<Parameters:>

    $outputHandleRef  - feature count file (e.g.: c22_genes_count.txt)
    $fileName         - feature count file
    $chrom            - chromosome name
    $featureStorage   - HASH MAP REF to features read with readFeatures procedure
    $rpkmFactor       - If undef, the legacy Illumina normalized value is printed instead of
                        RPKM. For genes and exons rpkmFactor is total number of mapped bases 
                        in the experiment. For splice junctions, the rpkmFactor is total number
                        of mapped reads in the experiment 

B<Returns:>
    nothing

=cut

sub printFeatureCounts($$\%$$) {
    my ( $outputHandleRef, $chrom, $featureStorage, $printScore, $rpkmFactor ) = @_;

    my @featuresIds = keys %{ $featureStorage->{feature} };

    ## assertion
    foreach my $featureId (@featuresIds) {
        if ( !defined $featureStorage->{feature}->{$featureId}->{start} ) {
            use Data::Dumper;
            print Dumper( $featureStorage->{feature}->{$featureId} );
            print "$featureId\n";
            errorExit "ERROR: saveFeatureCounts() critical problem\n";
        }
    }

    # end

    my @featuresSorted = sort {
        $featureStorage->{feature}->{$a}->{start} <=> $featureStorage->{feature}
          ->{$b}->{start}
    } @featuresIds;
    my %genes   = ();
    my $content = "";

    foreach my $featureId (@featuresSorted) {
        my $coverage    = 0;
        my $coverageAvg = 0;
        my $size        = 0;
        if ( defined $featureStorage->{feature}->{$featureId}->{coverage} ) {
            $coverage = $featureStorage->{feature}->{$featureId}->{coverage};
            if ( defined $featureStorage->{feature}->{$featureId}->{size}
                && $featureStorage->{feature}->{$featureId}->{size} > 0 )
            {

                $size        = $featureStorage->{feature}->{$featureId}->{size};
                $coverageAvg = $coverage / $size;
                # Normalize for RPKM if total number of mapped bases is given
                if (defined $rpkmFactor)
                {
                    $coverageAvg *= (10**9 / $rpkmFactor);
                }
            }
        }    # if

        if ( defined $printScore && $printScore > 0 ) {
            my $score = 0;
            $score = defined $featureStorage->{feature}->{$featureId}->{score};
            $content = sprintf(
                "%s\t%d\t%d\t%s\t%.5f\t%d\t%d\n",
                $chrom,
                $featureStorage->{feature}->{$featureId}->{start},
                $featureStorage->{feature}->{$featureId}->{end},
                $featureStorage->{feature}->{$featureId}->{name},
                $coverageAvg,
                $coverage,
                $score
            );
        }
        else {

            $content = sprintf(
                "%s\t%d\t%d\t%s\t%.5f\t%d\n",
                $chrom,
                $featureStorage->{feature}->{$featureId}->{start},
                $featureStorage->{feature}->{$featureId}->{end},
                $featureStorage->{feature}->{$featureId}->{name},
                $coverageAvg,
                $coverage
            );
        }

        print $outputHandleRef $content;
    }
}

=pod

=item saveFeatureCounts

The procedure saves feature counts (number of bases or reads starts in a feature) into a file
e.g.:
c22	15977189	15982257	ENSG00000183307	0.52	2648
c22	15998411	16026177	ENSG00000069998	0.78	21791

B<Parameters:>

    $featurefileName  - feature count file (e.g.: c22_genes_count.txt)
    $fileName         - feature count file
    $chrom            - chromosome name
    $featureStorage   - HASH MAP REF to features read with readFeatures procedure

B<Returns:>
    nothing

=cut

sub saveFeatureCounts($$\%$) {
    my ( $featurefileName, $chrom, $featureStorage, $printScore ) = @_;

    my $featurefile = IO::File->new( ">" . $featurefileName )
      || errorExit "ERROR: saveFeatureCounts() Couldn't create/open file handle for $featurefileName $!\n";

    printFeatureCounts($featurefile, $chrom, %$featureStorage, $printScore, undef);

    close $featurefile;
}

=pod

=item readFeaturesENST($featuFileName, $featureStorage)

The procedure reads features from a file in ENST format.

B<Parameters:>

    $featuFileName    -
    $featureStorage   -

B<Returns:>
    nothing

=cut

sub readFeaturesENST($;\%) {
    croak "ERROR: readFeaturesENST \n" unless ( @_ == 2 );
    my ( $featuFileName, $featureStorage ) = @_;
    $featureStorage->{count} = 0;
    open( FEATURE, "<$featuFileName" )
      || errorExit "ERROR: readFeaturesENST() Could not open $featuFileName: $!\n";
    while (<FEATURE>) {
        next if (/^#/);
        my %feature = ();
        if (/(\d+)\s+(\d+)\s+(\S+)/) {
            $feature{start} = $1;
            $feature{end}   = $2;
            $feature{name}  = $3;
            $feature{id}    = $feature{name};
            $feature{size}  = $feature{end} - $feature{start} + 1;
            $feature{name} =~ /^(S+)\./;
            my @ganes = ();
            push @ganes, $1;
            $feature{genes} = \@ganes;
        }    # if
        $featureStorage->{feature}->{ $feature{id} } = \%feature;
        push @{ $featureStorage->{positions}->{ $feature{start} } }, \%feature;
        $featureStorage->{count}++;
    }    # while
    close(FEATURE);
}    # readFeaturesENST

=pod

=item readFeaturesENST2( $chrom, $featuFileName, $featureStorage)

The procedure reads features from a file in Santa Cruz Browser format.
e.g.:
gene_stable_id	struct_external_gene_id	exon_stable_id	str_chrom_name	exon_chrom_start	exon_chrom_end	transcript_chrom_strand	struct_biotype
ENSG00000146556	Q7Z3N4_HUMAN	ENSE00001367146	1	4274	4365	-1	protein_coding
ENSG00000146556	Q7Z3N4_HUMAN	ENSE00001383334	1	4863	4901	-1	protein_coding
ENSG00000146556	Q7Z3N4_HUMAN	ENSE00001388009	1	5659	5764	-1	protein_coding
ENSG00000146556	Q7Z3N4_HUMAN	ENSE00001375216	1	5767	5810	-1	protein_coding

%featureFieldsENS2
	gene_stable_id   => 0,
	struct_gene_id   => 1,    
	exon_stable_id   => 2,
	str_chrom_name   => 3,
	exon_chrom_start => 4,
	exon_chrom_end   => 5,
	strand           => 6,
	struct_biotype   => 7,

B<Parameters:>

    $chrom            - if not "" then read only $chrom chromosome
    $featuFileName    -
    $featureStorage   -

B<Returns:>
    nothing

=cut

sub readFeaturesENST2($;$;\%) {
    croak "ERROR: readFeaturesENST2 \n" unless ( 3 == @_);
    my ( $chrom, $featuFileName, $featureStorage ) = @_;
    $featureStorage->{count} = 0;
    open( FEATURE, "<$featuFileName" )
      || errorExit "ERROR: readFeaturesENST2() Could not open $featuFileName: $!\n";
    while (<FEATURE>) {
        next if (/^#/);
        chomp;
        my @row = split "\t", $_;
        my $chromTmp = $row[ $featureFieldsENS2{str_chrom_name} ];
        if ( $chrom ne "" && $chrom eq $chromTmp ) {

            my %feature = ();
            $feature{chrom} = $chromTmp;
            my $start = 0;
            my $end   = 0;

            $start = $row[ $featureFieldsENS2{start} ];
            $end   = $row[ $featureFieldsENS2{end} ];

            $feature{start} = $start;
            $feature{end}   = $end;
            $feature{name}  = $row[ $featureFieldsENS2{exon_stable_id} ];
            $feature{size}  = $feature{end} - $feature{start} + 1;
            my @ganes = ();
            push @ganes, $row[ $featureFieldsENS2{gene_stable_id} ];
            $feature{genes} = \@ganes;

            $feature{id} = sprintf( "%s:%d-%d",
                $feature{name}, $feature{start}, $feature{end} );
            $featureStorage->{feature}->{ $feature{id} } = \%feature;
            push @{ $featureStorage->{positions}->{ $feature{start} } },
              \%feature;
            $featureStorage->{count}++;
        }    # if
    }    # while
    close(FEATURE);
}    # readFeatures

=pod

=item readFeaturesUCSCB($chrom, $featuFileName, $featureStorage)

The procedure reads features from a file in Santa Cruz Browser format.
e.g.:
chr1	58953	59871	OR4F5
chr1	357521	358460	OR4F16,OR4F29,OR4F3

B<Parameters:>

	$chrom            - if not "" then read only $chrom chromosome
    $featuFileName    -
    $featureStorage   -

B<Returns:>
    nothing

=cut

sub readFeaturesUCSCB($;$;\%) {
    croak "ERROR: readFeaturesUCSCB \n" unless ( 3 == @_ );
    my ( $chrom, $featuFileName, $featureStorage ) = @_;
    $featureStorage->{count} = 0;
    open( FEATURE, "<$featuFileName" )
      || errorExit "ERROR: readFeaturesUCSCB() Could not open $featuFileName: $!\n";
    while (<FEATURE>) {
        next if (/^#/);
        chomp;
        my @row = split "\t", $_;

        #my @names = split ",",  $row[ $refFlatExonsFields{names} ];
        #foreach my $geneName (@names) {
        my $chromTmp = $row[ $refFlatExonsFields{chrom} ];
        if ( $chrom ne "" && $chrom eq $chromTmp ) {

            #print "$chrom $_";
            my %feature = ();
            $feature{chrom} = $chromTmp;
            $feature{start} = $row[ $refFlatExonsFields{start} ];
            $feature{end}   = $row[ $refFlatExonsFields{end} ];
            $feature{name}  = $row[ $refFlatExonsFields{names} ];
            $feature{size}  = $feature{end} - $feature{start} + 1;
            my @ganes = ();
            push @ganes, split ",", $row[ $refFlatExonsFields{names} ];
            $feature{genes} = \@ganes;

            #$feature{name}  = $geneName;
            $feature{id} = sprintf( "%s:%d-%d",
                $feature{name}, $feature{start}, $feature{end} );
            $featureStorage->{feature}->{ $feature{id} } = \%feature;
            push @{ $featureStorage->{positions}->{ $feature{start} } },
              \%feature;
            $featureStorage->{count}++;

            #print "$feature{id} \n";
            #}    # foreach
            #return;
        }    # if
    }    # while
    close(FEATURE);
}    # readFeatures

=pod

=item readFasta($genomeSeqFile)

The procedure reads the sequence from the file (fasta format).

B<Parameters:>

    $genomeSeqFile    - sequence file name

B<Returns:>
    genome sequence

=cut

sub readFasta($) {
    my ($genomeSeqFile) = @_;
    my $genomeSeq = "";

    #get the genome reference seq
    open( GEN, "<$genomeSeqFile" )
      or errorExit "ERROR: readFasta() Can't open $genomeSeqFile $!\n";
    while (<GEN>) {
        chomp;
        next if ( $_ =~ m/>/ || length ($_) == 0);
        my $line = uc($_);
        if ( $line !~ /[ACGTN]+/ ) {
            errorExit "ERROR: Wrong fasta file(accepted are only ACGTN):\n[$line]\n ";
        }    # if
        $genomeSeq .= $_;
    }    # while
    close(GEN);
    return $genomeSeq;
}    #readFasta

=pod

=item parseSpliceName($spliceJunctioName, $spliceJunctioDesc)

The procedure pareses splice jucntion name and returns information 
about scaffold name start stop and splice site size.


B<Parameters:>

    $spliceJunctioName    - splice Junction Name
    $spliceJunctioDesc    - HASH MAP REF returns the description of splice junction 

B<Returns:>
    genome sequence
        


=cut

sub parseSpliceName($\%) {
    my ( $spliceJunctioName, $spliceJunctioDesc ) = @_;
    if ( $spliceJunctioName =~
        /^(\S+)_([0-9]+)\-([0-9]+)_([0-9]+)\-([0-9]+)_(\S+)_([0-9]+)\-([0-9]+)$/
      )
    {
        $spliceJunctioDesc->{length1}  = $2;
        $spliceJunctioDesc->{length2}  = $3;
        $spliceJunctioDesc->{scaffold} = $6;
        $spliceJunctioDesc->{start}    = $7;
        $spliceJunctioDesc->{end}      = $8;
    }
    elsif ( $spliceJunctioName =~
        /^(.+)_([0-9]+)\_([0-9]+)_(\S+)_([0-9]+)\_([0-9]+)$/ )
    {
        my $scaffoldFile = $4;
        
        $spliceJunctioDesc->{length1} = $2;
        $spliceJunctioDesc->{length2} = $3;
        $spliceJunctioDesc->{start}   = $5;
        $spliceJunctioDesc->{end}     = $6;

        $spliceJunctioDesc->{scaffold} = $scaffoldFile;
        
    }
    else {
        errorExit "ERROR: parseSpliceName() [$spliceJunctioName]\n";
    }    # else
    $spliceJunctioDesc->{size} =
      $spliceJunctioDesc->{length1} + $spliceJunctioDesc->{length2};
    $spliceJunctioDesc->{name} = $spliceJunctioName;
    $spliceJunctioDesc->{id}   = $spliceJunctioName;
}

=pod

=item outputSEread($ref2Output, $refDirs, $start, 
        $chr, $ref2read2, $fileType)

The procedure prints one reads into bin files.

B<Parameters:>

    $exportStorageRef - HASH MAP REF with file handles 
    $start            - start position (of read1)
    $chr              - chromosome folder name
    $ref2read1        - ARRAY REF with read
    $fileType         - type of read (each type goes to diffrent file)

B<Returns:>
    nothing

=cut

## TODO marge with writeRead
sub outputSEread(\%$$\@$) {

    my ( $exportStorageRef, $start, $chr, $ref2read1, $fileType ) = @_;

    if ( $fileType eq 'single' ) {
        my $bin    = sprintf "%04d", $start / $exportStorageRef->{binSize};
        my $binDir = File::Spec->catdir($exportStorageRef->{resultsDir}, $chr, $bin);
        my $file   = File::Spec->catfile($binDir,"single_export.txt");

        my $ref2Output = $exportStorageRef->{ref2Output};
        my $thisFile = $ref2Output->{files}->{$file};
        unless (defined $thisFile)
        {
            if (!($thisFile = IO::File->new(">$file")))
            {
                if( !exists $chrEnds{$chr} )
                {
                    errorExit "Chromosome name $chr in input data is unexpected "
                             ."while processing the following read:"
                             .join( "\t", @{$ref2read1});
                }
                errorExit "ERROR: Couldn't create file handle for $file $!\n";
            }

            $ref2Output->{files}->{$file} = $thisFile;
        }
        print $thisFile join( "\t", @{$ref2read1} ), "\n";
    }
    else {
        errorExit "ERROR outputSEread unsuported fileType [$fileType]\n";
    }    # else
}    # outputSEread

=pod

=item outputPEread($exportStorageRef, $start, $end, $chr, $ref2read1, $ref2read2, $fileType)

The procedure prints two reads into bin files.

B<Parameters:>

    $exportStorageRef - HASH MAP REF with file handles and buffers 
    $start            - start position (of read1 or read2)
    $end              - end position (of read1 or read2)
    $chr              - chromosome folder name
    $ref2read1        - ARRAY REF with read1
    $ref2read2        - ARRAY REF with read2
    $fileType         - type of read (each type goes to diffrent file)

B<Returns:>
    nothing

=cut

sub outputPEread(\%$$$\@\@$) {

    my ( $exportStorageRef, $start, $end, $chr, $ref2read1, $ref2read2,
        $fileType )
      = @_;

    my $bin    = sprintf "%04d", $start / $exportStorageRef->{binSize};
    my $fileHashKey = $chr . $bin . $fileType;
    my $ref2Output = $exportStorageRef->{ref2Output};

    if($exportStorageRef->{isCompressPair}) {
        # the file split is done by compressXPair here, but perl has
        # more context for error messages so we maintain per-file
        # checking code in this option:
        #
        unless(exists $ref2Output->{files}->{$fileHashKey}){
            my $binDir = File::Spec->catdir($exportStorageRef->{resultsDir}, $chr, $bin);
            unless(-d $binDir) {
                unless(exists $chrEnds{$chr}) {
                    errorExit "ERROR: Chromosome name $chr in input data is unexpected "
                        ."while processing the following read pair:"
                        .join( "\t", $start, $end, @{$ref2read1}, @{$ref2read2});
                }
                errorExit "ERROR: Can't find expected directory '$binDir'\n";
            }
            $ref2Output->{files}->{$fileHashKey} = undef;
        }

        print ${$exportStorageRef->{compFH}} ( join( "\t", "$chr/$bin" , $fileType, $start, $end, @{$ref2read1}, @{$ref2read2} )  . "\n");
    } else {
        my $thisFile =  $ref2Output->{files}->{$fileHashKey};
        unless (defined $thisFile)
        {
            my $binDir = File::Spec->catdir($exportStorageRef->{resultsDir}, $chr, $bin);
            my $file   = '';
            if ( $fileType eq 'norm' ) {
                $file = File::Spec->catfile($binDir,"unsort.txt");
            }
            elsif ( $fileType eq 'orphan' ) {
                $file = File::Spec->catfile($binDir,"unsort_orph.txt");
            }
            elsif ( $fileType eq 'anom' ) {
                $file = File::Spec->catfile($binDir,"unsort_anom.txt");
            }
            else {
                errorExit "Invalid read file type '$fileType'";
            }    # else

            my $cmd = "> $file";
            unless( open($thisFile,$cmd) )
            {
                unless( exists $chrEnds{$chr} )
                {
                    errorExit "ERROR: Chromosome name $chr in input data is unexpected "
                        ."while processing the following read pair:"
                        .join( "\t", $start, $end, @{$ref2read1}, @{$ref2read2});
                }
                errorExit "ERROR: Couldn't open output file '$file' $!\n";
            }

            $ref2Output->{files}->{$fileHashKey} = $thisFile;
        }
        print $thisFile  ( join( "\t", $start, $end, @{$ref2read1}, @{$ref2read2} )  . "\n");
    }
}     # outputPERead
1;    # says use was ok
__END__

