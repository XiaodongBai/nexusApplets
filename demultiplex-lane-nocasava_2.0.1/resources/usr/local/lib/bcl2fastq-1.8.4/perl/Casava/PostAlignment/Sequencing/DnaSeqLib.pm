package Casava::PostAlignment::Sequencing::DnaSeqLib;

# PROJECT: CASAVA
# MODULE:  $RCSfile: DnaSeqLib.pm,v $
# AUTHOR:  Lukasz Szajkowski, Richard Carter
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#
# The library contains procedures and variables specific to DNA-Seq analysis
# and some additional general purpose functions.

=pod

=head1 NAME

Casava::PostAlignment::Sequencing::DnaSeqLib.pm - Perl utility library for running BullForg analysis.

=head1 SYNOPSIS

The library contains procedures and variables specific to CASAVA
application. 
 
use Casava::PostAlignment::Sequencing::DnaSeqLib.pm qw();  

=head1 DESCRIPTION

Exports: 
	positionCoverage(\%;\%;\@);
    
# Global variable
    %predictedMateStrand;



=head1 AUTHORSHIP

Copyright (c) 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

Created by Lukasz Szajkowski <lszajkowski@illumina.com>

=cut

BEGIN {
    use Exporter();
    @ISA       = qw(Exporter);
    @EXPORT    = qw($VERSION);
    @EXPORT_OK =
      qw(&positionCoverage
      &summaryHtmlToShortXml
      &summaryXmlToShortXml
       &splitSortNormInit &splitSortNormRecord &splitSortNormClose
       &splitSortOrphInit &splitSortOrphRecord &splitSortOrphClose
       &splitSortAnomInit &splitSortAnomRecord &splitSortAnomClose);
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);

use Carp;

use Casava::Common::Log;
use Casava::PostAlignment::Sequencing::GenomicIOLib qw(writeRead2Scaffold writeRead configureReadStorage
  flushReadFiles %positionFields parseSpliceName  %exportFields);
use Casava::Common::IOLib qw(bufferedPrint executeCmd formatValueWithDeviation);
use Casava::PostAlignment::Sequencing::Config qw(%CONF_APP isSpliceJunctionChrom);

sub splitSortFilesNorm($$$$\%);
sub splitSortFilesOrph($$$$\%);
sub splitSortFilesAnom($$$$$\%);
sub summaryHtmlToShortXml(\%;\%;);
sub summaryXmlToShortXml($;\%;);
sub positionCoverage($;\%;\%;\@);

our %predictedMateStrand = (
    'F' => 'R',
    'R' => 'F'
);

our @strands = ( 'F', 'R' );
our @bases   = ( 'A', 'C', 'G', 'T' );

our %htmlTablesNames = (
    'Lane Parameter Summary'         => 'lane_parameter',
    'Lane Results Summary : Read 1'  => 'lane_results_r1',
    'Lane Results Summary : Read 2'  => 'lane_results_r2',
    'Expanded Lane Summary : Read 1' => 'expanded_lane_r1',
    'Expanded Lane Summary : Read 2' => 'expanded_lane_r2',
    'Lane Results Summary'           => 'lane_results',
    'Expanded Lane Summary'          => 'expanded_lane'
);
our %htmlFieleds = (
    'Equiv Perfect Clusters (PF)'  => 'equiv_perfect_clusters_pf',
    'Sample ID'                    => 'sample_id',
    'Clusters (tile mean) (raw)'   => 'clusters_tile_mean_raw',
    '% Prephasing'                 => 'prephasing_pr',
    'Clusters (PF)'                => 'clusters_pf',
    '% Error Rate (raw)'           => 'error_rate_raw_pr',
    '% Error Rate (PF)'            => 'error_rate_pf_pr',
    'Equiv Perfect Clusters (raw)' => 'equiv_perfect_clusters_raw',
    'Lane Yield (kbases)'          => 'lane_yield_kb',
    'Yield (kbases)'               => 'yield_kb',
    'Num Tiles'                    => 'numTiles',
    'Cycle 10-20 Av % Loss (PF)'   => 'cycle_10_20_av_loss_pf_pr',
    'Cycle 2-4 Av Int (PF)'        => 'cycle_2_4_av_int_pf',
    '% Align (PF)'                 => 'align_pf_pr',
    'Tiles'                        => 'tiles',
    'Cycle 2-10 Av % Loss (PF)'    => 'cycle_2_10_av_loss_pf_pr',

    #    'Clusters (PF)'                => 'clusters_PF',
    'Sample Target'        => 'sample_target',
    'Alignment Score (PF)' => 'alignment_score_pf',
    'Filter'               => 'filter',
    'Clusters'             => 'clusters',
    '% retained'           => 'retained_pr',
    '1st Cycle Int (PF)'   => 'firts_st_cycle_int_pf',
    'Lane'                 => 'lane',

    #	'% Error Rate (PF)'            => 'error_rate_pf',
    #    'Lane'                         => 'lane',
    'Sample Type'                      => 'sample_type',
    'Length'                           => 'length',
    '% PF Clusters'                    => 'pf_clusters_pr',
    '% Phasing'                        => 'phasing_pr',
    'Chast. Thresh.'                   => 'chast_thresh',
    'Chastity Threshold'               => 'chast_thresh',
    '% intensity after 20 cycles (PF)' => 'intensity_after_20_cycles_pf',
    'Clusters (raw)'                   => 'clusters_raw'
);



sub splitSortNormInit($$$\%){
    my ( $binSize, $resultsDir, $fieldsRef, $stateRef ) = @_;

    my %correctReadConf = (
        binSize  => $binSize,
        readType => 'corect_sort',
        baseDir  => $resultsDir,
        fields   => $fieldsRef
    );
    my %wrongReadConf = (
        binSize  => $binSize,
        readType => 'unsort.wrong',
        baseDir  => $resultsDir,
        fields   => $fieldsRef
    );

    $stateRef-> {correctStorage} =
      configureReadStorage( %correctReadConf );
    $stateRef-> {wrongStorage} =
      configureReadStorage( %wrongReadConf );
}


my $ssfield=$exportFields{SingleScore};

sub splitSortNormRecord($$;$) {
    my ( $tempRef, $stateRef, $extraVal ) = @_;

    my @read1 = @{$tempRef}[ 0 .. 21 ];
    my @read2 = @{$tempRef}[ 22 .. 43 ];

    if( defined $extraVal ) {
        push @read1, $extraVal;
        push @read2, $extraVal;
    }

    # Additional optional field check in Norm only:
    if( ($read1[$ssfield] == 0) &&
        ($read2[$ssfield] == 0) ) {
        push @read1, "AM:i:0";
        push @read2, "AM:i:0";
    }

    my $correct = $stateRef-> {correctStorage};
    my $wrong = $stateRef-> {wrongStorage};

    if ( $read1[ $correct->{config}->{fields}->{Strand} ] eq 'F' ) {
        writeRead( %{$correct}, @read1 );
        writeRead( %{$wrong},  @read2 );
    }
    else {
        writeRead( %{$correct}, @read2 );
        writeRead( %{$wrong},  @read1 );
    }
}



sub splitSortNormClose(\%) {
    my ( $stateRef ) = @_;

    flushReadFiles( %{$stateRef->{correctStorage}}, 1 );
    flushReadFiles( %{$stateRef->{wrongStorage}},  1 );
}



=pod

=item splitSortFilesNorm

The procedure splits the sort.txt file from rmDupFromSort into sorted and
unsorted reads. Unsorted reads are split into chromosomes.

Note this procedure is no longer used directly, but is left here to
document the basic usage of the Init,Record and Close functions to
split the reads in a file.

B<Parameters:>

    $binSize    - build bin size
    $sortPref   - sort file prefix (eg.: sort)
    $sourceDir  - source directory
    $resultsDir - result directory
    $fieldsRef  - Hash map with double export format

B<Returns:> 
    nothing

=cut

sub splitSortFilesNorm($$$$\%) {
    my ( $binSize, $sortPref, $sourceDir, $resultsDir, $fieldsRef ) = @_;

    my $sortedPairFilePath = File::Spec->catfile( $sourceDir , $sortPref . ".txt" );
    open( SRT,  "<$sortedPairFilePath" )
      or errorExit
      "ERROR: splitSortFilesNorm() Couldn't open $sortedPairFilePath $!\n";

    my %state;
    splitSortNormInit($binSize, $resultsDir, $fieldsRef, %state);

    while (<SRT>) {
        chomp;
        my @temp  = split /\t/, $_;
        splitSortNormRecord(\@temp,\%state);
    }    # while
    close(SRT);
    splitSortNormClose(%state);
}




sub splitSortOrphInit($$$\%) {
    my ( $binSize, $resultsDir, $fieldsRef, $stateRef ) = @_;

    my %correctReadConf = (
        binSize  => $binSize,
        readType => 'corect_sort_orph',
        baseDir  => $resultsDir,
        fields   => $fieldsRef
    );

    $stateRef-> {correctStorage} =
      configureReadStorage( %correctReadConf );
}



sub splitSortOrphRecord($$;$) {
    my ( $tempRef, $stateRef, $extraVal ) = @_;

    my @read1 = @{$tempRef}[ 0 .. 21 ];
    my @read2 = @{$tempRef}[ 22 .. 43 ];

    if(defined $extraVal) {
        push @read1, $extraVal;
        push @read2, $extraVal;
    }

    my $correctRef = $stateRef->{correctStorage};
    my $fieldsRef = $correctRef->{config}->{fields};

    my $position1 = $read1[ $fieldsRef->{Posn} ];

    my $binId = sprintf "%04d", $position1 / $correctRef->{config}->{binSize};
    if ( (not defined $stateRef->{prevBinId}) or ($stateRef->{prevBinId} != $binId) ) {
        if (defined $stateRef->{prevBinId}) {
            flushReadFiles( %{$correctRef}, 1 );
        }
        $stateRef->{prevBinId} = $binId;
    }

    $read2[ $fieldsRef->{Posn} ] = $position1;
    $read2[ $fieldsRef->{Strand} ] =
      $predictedMateStrand{ $read1[ $fieldsRef->{Strand} ] };
    $read2[ $fieldsRef->{Descriptor} ] = '-';
    $read2[ $fieldsRef->{SingleScore} ] = -1;
    $read2[ $fieldsRef->{PairedScore} ] = 0;
    $read2[ $fieldsRef->{PartnerStrand} ] = 'N';
    $read2[ $fieldsRef->{PartnerOffset} ] = 0;

    writeRead( %{$correctRef}, @read1 );
    writeRead2Scaffold( %{$correctRef}, $read1[ $fieldsRef->{Chr} ], @read2 );
}



sub splitSortOrphClose(\%) {
    my ( $stateRef ) = @_;

    flushReadFiles( %{$stateRef->{correctStorage}}, 1 );
}



=pod

=head1 The splitSortFilesOrph splits the sort_orph.txt.



=item splitSortFilesOrph

The procedure splits the sort_orph.txt file from rmDupFromSort into sorted and
unsorted reads. Unsorted reads are split into chromosomes.

Note this procedure is no longer used directly, but is left here to
document the basic usage of the Init,Record and Close functions to
split the reads in a file.

B<Parameters:>
    $binSize    - build bin size
    $sortPref   - sort file prefix (eg.: sort)
    $sourceDir  - source directory
    $resultsDir - result directory
    $fieldsRef  - Hash map with double export format
B<Returns:>
    nothing


=cut

sub splitSortFilesOrph($$$$\%) {
    my ( $binSize, $sortPref, $sourceDir, $resultsDir, $fieldsRef ) = @_;

    my $sortedPairFilePath = File::Spec->catfile($sourceDir, "$sortPref.txt");
    open( SRT, "<$sortedPairFilePath" )
      or errorExit "$0:ERROR: Couldn't open $sortedPairFilePath $!\n";

    my %state;
    splitSortOrphInit($binSize, $resultsDir, $fieldsRef, %state);

    while (<SRT>) {
        chomp;
        my @temp = split("\t");
        splitSortOrphRecord(\@temp,\%state);
    }
    close(SRT);
    splitSortOrphClose(%state);
}



sub splitSortAnomInit($$$$\%){
    my ( $currentChrom, $binSize, $resultsDir, $fieldsRef, $stateRef ) = @_;

    $stateRef->{currentChrom} = $currentChrom;

    my %correctReadConf = (
        binSize  => $binSize,
        readType => 'corect_sort_anom',
        baseDir  => $resultsDir,
        fields   => $fieldsRef
    );
    my %wrongReadConf = (
        binSize  => $binSize,
        readType => 'anom.unsort.wrong',
        baseDir  => $resultsDir,
        fields   => $fieldsRef
    );

    $stateRef->{correctStorage} =
      configureReadStorage( %correctReadConf );
    $stateRef->{wrongStorage} =
      configureReadStorage( %wrongReadConf );
}



sub splitSortAnomRecord($$;$) {
    my ( $tempRef, $stateRef, $extraVal ) = @_;

    my @read1 = @{$tempRef}[ 0 .. 21 ];
    my @read2 = @{$tempRef}[ 22 .. 43 ];

    if(defined $extraVal) {
        push @read1, $extraVal;
        push @read2, $extraVal;
    }

    my $correctRef = $stateRef->{correctStorage};
    my $wrongRef = $stateRef->{wrongStorage};

    my $fieldsRef = $correctRef->{config}->{fields};

    my $chr1 = $read1[ $fieldsRef->{Chr} ];
    my $chr2 = $read2[ $fieldsRef->{Chr} ];

    my $read = 0;
    if ( $read1[ $fieldsRef->{Strand} ] ne $read2[ $fieldsRef->{Strand} ] ) {
        if ( $chr1 ne $chr2 ) {
            # cross chromosome match
            # put it in the anom file to be sorted
            # store everything by its first read chromosome
            $read = 1 if ( $stateRef->{currentChrom} eq $chr1 );
        } else {
            $read = ( ($read1[$fieldsRef->{Strand}] eq 'F') ? 1 : 2 );
        }

        if ( $read != 0 ) {
            my $binReadRef = ( ($read == 1) ? \@read1 : \@read2 );
            my $binId = sprintf "%04d", $binReadRef->[ $fieldsRef->{Posn} ] / $correctRef->{config}->{binSize};
            if ( (not defined $stateRef->{prevBinId}) or ($stateRef->{prevBinId} != $binId) ) {
                if (defined $stateRef->{prevBinId}) {
                    flushReadFiles( %{$correctRef}, 1 );
                }
                $stateRef->{prevBinId} = $binId;
            }
        }
    }

    my $read1Ref = ( ($read==1) ? $correctRef : $wrongRef );
    my $read2Ref = ( ($read==2) ? $correctRef : $wrongRef );
    writeRead( %{$read1Ref}, @read1 );
    writeRead( %{$read2Ref}, @read2 );
}



sub splitSortAnomClose(\%) {
    my ( $stateRef ) = @_;

    flushReadFiles( %{$stateRef->{correctStorage}}, 1 );
    flushReadFiles( %{$stateRef->{wrongStorage}},  1 );
}



=pod

=head1 The procedure splits the sort_anom.txt.



=item splitSortFilesAnom

The procedure splits the sort_anom.txt file from rmDupFromSort into sorted and
unsorted reads. Unsorted reads are split into chromosomes.

Note this procedure is no longer used directly, but is left here to
document the basic usage of the Init,Record and Close functions to
split the reads in a file.

B<Parameters:>
    $currentChrom - currently proccess chromosome name
    $binSize      - build bin size
    $sortPref     - sort file prefix (eg.: sort)
    $sourceDir    - source directory
    $ResultsDir   - result directory
    $fieldsRef    - Hash map with double export format

B<Returns:> 
    nothing
        


=cut

sub splitSortFilesAnom($$$$$\%) {
    my (
        $currentChrom, $binSize, $sortPref,
        $sourceDir, $resultsDir,
        $fieldsRef
      )
      = @_;

    my $sortedPairFilePath = File::Spec->catfile($sourceDir, "$sortPref.txt");
    open( SRT, "<$sortedPairFilePath" )
      or errorExit "$0:ERROR: Couldn't open $sortedPairFilePath $!\n";

    my %state;
    splitSortAnomInit($currentChrom, $binSize, $resultsDir, $fieldsRef, %state);

    while (<SRT>) {
        chomp;
        my @temp = split /\t/, $_;
        splitSortAnomRecord(@temp,%state);
    }
    close(SRT);
    splitSortAnomClose(%state);
}




=pod

=head1 The procedure calculates the mean coverage over all the features 
	in the input HASH MAP over one position.



=item calculatesMeanCove($inputFile, $chromosomeEnd, $dirBuildParsed, 
	$chrom, $binSize, $minStep, $allele_sort_f)

The procedure calculates the mean coverage over all the features 
	in the input HASH MAP over one position.

B<Parameters:>
	$threshold			  - ignore bases with score below threshold
    $featureStorage       - HASH MAP REF with features read with readFeatures procedure
    $tempFeatureStorage   - temporary storage HASH MAP REF (can be empty)
    $positionArray        - split row from sort.count file in $positionFields format
B<Returns:> 
    nothing
        


=cut

our @posFeaturesSorted     = ();
our $posFeaturesSortedSize = 0;
our $isInit                = 0;
our $currentCountPosIndex  = 0;
our $highestEnd            = 0;

sub positionCoverage($;\%;\%;\@) {
    croak "ERROR: positionCoverage \n" unless ( @_ == 4 );
    my ( $threshold, $featureStorage, $tempFeatureStorage, $positionArray ) =
      @_;
    my $position = $positionArray->[ $positionFields{position} ];

    #return;
    #my $position;
    my $featureId;
    my $isFound;
    my $stop = 0;
    if ( $isInit == 0 ) {

        #my @features          = keys %{ $tempFeatureStorage->{feature} };
        my @tmp = keys %{ $featureStorage->{positions} };
        @posFeaturesSorted     = sort { $a <=> $b } @tmp;
        $posFeaturesSortedSize = scalar(@posFeaturesSorted);
        $isInit                = 1;
        $currentCountPosIndex  = 0;
    }

    #	print "$currentCountPosIndex $posFeaturesSortedSize $stop\n ";
    my $loopsize = 0;

    for (
        my $i = $currentCountPosIndex ;
        ( $i < $posFeaturesSortedSize && $stop == 0 ) ;
        $i++
      )
    {
        my $featurePos = $posFeaturesSorted[$i];
        if ( $position < $featurePos ) {
            $stop = 1;
            next;
        }
        my $featuresStartingAt = $featureStorage->{positions}->{$featurePos};
        foreach my $featureRef (@$featuresStartingAt) {
            $featureId = $featureRef->{id};
            $isFound   = 0;

            if (
                $featureStorage->{feature}->{$featureId}->{start} <= $position )
            {
                if ( $featureStorage->{feature}->{$featureId}->{end} >=
                    $position )
                {

                    if ( $positionArray->[ $positionFields{num_call} ] >
                        $threshold )
                    {
                        if ( $isFound == 0 ) {
                            $isFound = 1;
                        }
                        $featureStorage->{feature}->{$featureId}->{coverage} +=
                          $positionArray->[ $positionFields{num_call} ];
                        $featureStorage->{feature}->{$featureId}->{size}++;
                    }
                    else {
                        $featureStorage->{feature}->{$featureId}->{coverage} +=
                          0;
                    }

                    #print "loop $loopsize\n";
                }
                else {

                    #if ($featureStorage->{feature}->{$featureId}->{end})
                    $currentCountPosIndex = $i;
                }
            }
            else {

            }
        }
        $loopsize++;
    }
    return;
}

=pod

=item summaryHtmlToShortXml( $htmlMapRef, $xmlMapRef )

The The procedure parses Summary.html in HASH MAP and returns shorter version 
	(only lanes) as xml HASH MAP.

B<Parameters:>
    $htmlMapRef - in HASH MAP REF with Summary.html data
    $xmlMapRef  - out HASH MAP REF with short Summary.xml data           

B<Returns:> 
    nothing   

=cut

sub summaryHtmlToShortXml(\%;\%;) {
    croak "ERROR: summaryHtmlToShortXml " unless ( @_ == 2 );
    my ( $htmlMapRef, $xmlMapRef ) = @_;
    foreach my $bodyKey ( sort keys %{$htmlMapRef} ) {
        my @topElements    = sort keys %{ $htmlMapRef->{$bodyKey} };
        my $tableTitle     = '';
        my @tableList      = ();
        my %tableMap       = ();
        my $disableSection = 1;
        foreach my $topElement ( sort keys %{ $htmlMapRef->{$bodyKey} } ) {
            if ( $topElement eq 'table' ) {
                my @bigTable = @{ ${ $htmlMapRef->{$bodyKey} }{$topElement} };
                my $tableId  = 0;
                foreach my $eleRef (@bigTable) {
                    my ( $tableKey, $tableRef );
                    foreach $tableKey ( keys %{$eleRef} ) {
                        if ( ref( ${$eleRef}{$tableKey} ) eq 'ARRAY'
                            && defined $tableList[$tableId] )
                        {
                            my @rows       = @{ ${$eleRef}{$tableKey} };
                            my $tableTitle = $tableList[$tableId];

                            #print "$tableId $tableTitle\n";
                            my $tableName = $tableTitle;
                            if ( defined $htmlTablesNames{$tableTitle} ) {
                                $tableName      = $htmlTablesNames{$tableTitle};
                                $disableSection = 0;
                            }
                            else {
                                $disableSection = 1;
                            }
                            my %table = ();
                            if ( $tableTitle eq 'Chip Summary' ) {
                                $xmlMapRef->{chip}{machine} =
                                  ${$eleRef}{$tableKey}[0]{td}[1];
                                $xmlMapRef->{chip}{run_folder} =
                                  ${$eleRef}{$tableKey}[1]{td}[1];
                                $xmlMapRef->{chip}{chip_id} =
                                  ${$eleRef}{$tableKey}[2]{td}[1];
                            }
                            elsif ( $tableTitle eq 'Chip Results Summary' ) {
                                $xmlMapRef->{chip_results}{clusters} =
                                  ${$eleRef}{$tableKey}[1]{td}[0];
                                $xmlMapRef->{chip_results}{clusters_pf} =
                                  ${$eleRef}{$tableKey}[1]{td}[1];
                                $xmlMapRef->{chip_results}{yield_kb} =
                                  ${$eleRef}{$tableKey}[1]{td}[2];
                            }
                            else {
                                my $rowId      = 0;
                                my @fields     = ();
                                my $headerSize = 2;
                                if ( $disableSection == 1 ) {
                                    ## Skip unknown sections
                                    next;
                                }

                                $table{index} = $tableId - 2;
                                $table{title} = $tableTitle;
                                $table{name}  = $tableName;
                                $tableName =~ s/\s/_/g;
                                if (   $tableTitle eq 'Lane Parameter Summary'
                                    || $tableTitle =~ /Expanded Lane Summary/ )
                                {
                                    $headerSize = 1;
                                }

                 #                                if ( $tableTitle =~ /Lane/ ) {
                 #print "$tableTitle \n";
                 #                                    $headerSize = 1;
                 #next;
                 #                                }
                                if ( $tableTitle =~ /Pair Summary/ ) {

                                    #print "Pair Summary\n";
                                    next;
                                }
                                foreach my $rowRef (@rows) {
                                    my $cellId = 0;
                                    my $laneId = -1;
                                    my %newRow = ();
                                    foreach my $rowRef2 ( values %{$rowRef} ) {
                                        if ( ref($rowRef2) eq 'ARRAY' ) {
                                            if ( $rowId < $headerSize - 1 ) {
                                                next;
                                            }
                                            foreach my $cell ( @{$rowRef2} ) {
                                                if ( $headerSize - 1 == $rowId )
                                                {

                                                    if ( ref($cell) eq 'HASH' )
                                                    {
                                                        $headerSize = 2;
                                                        last;
                                                    }
                                                    $cell =~ s/^\s+//g;
                                                    $cell =~ s/\s+$//g;

                                                    if (
                                                        defined
                                                        $htmlFieleds{$cell} )
                                                    {
                                                        if ( $htmlFieleds{$cell}
                                                            eq
                                                            "error_rate_pf_pr"
                                                            && $tableName =~
                                                            "expanded_lane" )
                                                        {
                                                            $cell =
                                                              "error_rate_pf";
                                                            push @fields,
                                                              "error_rate_pf";
                                                        }
                                                        else {
                                                            push @fields,
                                                              $htmlFieleds{
                                                                $cell };
                                                        }
                                                    }
                                                }
                                                else {
                                                    if ( !defined
                                                        $fields[$cellId] )
                                                    {

                         #print "Warning: Someting is wrong with Summary.htm\n";
                                                        next;

                       #print join (@fields, ",") . "\n";
                       #print "$cellId $rowId $tableTitle $tableName $cellId\n";
                                                    }
                                                    if ( $fields[$cellId] eq
                                                        'tiles' )
                                                    {
                                                        next;
                                                    }
                                                    if ( $cellId == 0 ) {
                                                        $laneId = $cell;
                                                        if ( $laneId !~ /\d+/ )
                                                        {
                                                            next;
                                                        }
                                                        $newRow{ $fields
                                                              [ $cellId ] } =
                                                          $cell;
                                                    }
                                                    else {
                                                        $newRow{ $fields
                                                              [ $cellId ] } =
                                                          $cell;
                                                    }
                                                }
                                                $cellId++;
                                            }
                                        }
                                        if ( scalar keys(%newRow) > 0
                                            && $newRow{lane} < 9 )
                                        {
                                            push @{ $table{lane} }, \%newRow;
                                        }
                                    }
                                    $rowId++;
                                }

                                #print "[$tableName]\n";
                                $xmlMapRef->{$tableName} = \%table;

                                #push @{ $xmlMapRef->{table} }, \%table;
                            }
                        }
                    }
                    $tableId++;
                }
            }
            if ( $topElement eq 'h1' ) {
                $tableTitle = ${ $htmlMapRef->{$bodyKey} }{$topElement};
            }
            if ( $topElement eq 'h2' ) {
                my @titleTable = @{ ${ $htmlMapRef->{$bodyKey} }{$topElement} };
                for ( my $i = 0 ; $i < scalar(@titleTable) ; $i++ ) {
                    my $title;
                    if (
                        ref(
                            ${ ${ $htmlMapRef->{$bodyKey} }{$topElement} }[$i]
                        ) eq 'HASH'
                      )
                    {
                        $title =
                          ${ ${ $htmlMapRef->{$bodyKey} }{$topElement} }[$i]
                          ->{content};
                    }
                    else {
                        $title =
                          ${ ${ $htmlMapRef->{$bodyKey} }{$topElement} }[$i];
                    }
                    $tableMap{$title} = $i;
                    push @tableList, $title;
                }
            }
            if ( $topElement eq 'h3' ) {
                $tableTitle = ${ $htmlMapRef->{$bodyKey} }{$topElement};
            }
        }
    }
    return;
}

=pod

=item summaryXmlToShortXml( $summaryXMl, $xmlMapRef )

The The procedure parses Summary.xml in HASH MAP and returns shorter version 
	(only lanes) as xml HASH MAP.

B<Parameters:>
    $summaryXMl - in path to Summary.xml
    $xmlMapRef  - out HASH MAP REF with short runs_summary.xml data           

B<Returns:> 
    nothing   

=cut

sub summaryXmlToShortXml($;\%;) {
    croak "ERROR: summaryXmlToShortXml " unless ( @_ == 2 );
    my ( $summaryXMl, $xmlMapRef ) = @_;
    
    my $xs = new XML::Simple(forcearray => [ qw(Read Lane) ]);
    my $summaryXmlMapRef = $xs->XMLin($summaryXMl);
    
    my $readsNumber = scalar(@{ $summaryXmlMapRef->{LaneResultsSummary}{Read} });
    
    my $index = 0;
    $xmlMapRef->{chip}{machine} = $summaryXmlMapRef->{ChipSummary}{Machine};
    $xmlMapRef->{chip}{run_folder} = $summaryXmlMapRef->{ChipSummary}{RunFolder};
    $xmlMapRef->{chip}{chip_id} = $summaryXmlMapRef->{ChipSummary}{ChipID};

    $xmlMapRef->{chip_results}{clusters} = $summaryXmlMapRef->{ChipResultsSummary}{clusterCountRaw};
    $xmlMapRef->{chip_results}{clusters_pf} = $summaryXmlMapRef->{ChipResultsSummary}{clusterCountPF};
    $xmlMapRef->{chip_results}{yield_kb} = sprintf "%.0f", ( $summaryXmlMapRef->{ChipResultsSummary}{yield} / 1000 );

    $xmlMapRef->{lane_parameter}{name} = "lane_parameter"; 
    $xmlMapRef->{lane_parameter}{title} = "Lane Parameter Summary"; 
    $xmlMapRef->{lane_parameter}{index} = $index++; 
    foreach my $lane ( @{ $summaryXmlMapRef->{LaneParameterSummary}{Lane} } ) {
        my $laneNumber = $lane->{laneNumber} - 1;
        $xmlMapRef->{lane_parameter}{lane}[$laneNumber]{chast_thresh} = $lane->{chastityThreshold};
        $xmlMapRef->{lane_parameter}{lane}[$laneNumber]{filter} = $lane->{purity};
        $xmlMapRef->{lane_parameter}{lane}[$laneNumber]{lane} = $laneNumber + 1;
        $xmlMapRef->{lane_parameter}{lane}[$laneNumber]{length} = ("unknown" eq $lane->{lengthsList}) ? $lane->{originalReadLength} : $lane->{lengthsList};
        $xmlMapRef->{lane_parameter}{lane}[$laneNumber]{numTiles} = $lane->{tileCountRaw};
        $xmlMapRef->{lane_parameter}{lane}[$laneNumber]{sample_id} = $lane->{sample};
        $xmlMapRef->{lane_parameter}{lane}[$laneNumber]{sample_target} = $lane->{template};
        $xmlMapRef->{lane_parameter}{lane}[$laneNumber]{sample_type} = $lane->{type};
    }
    
    foreach my $read ( @{ $summaryXmlMapRef->{LaneResultsSummary}{Read} } ) {
        my $readNumber = $read->{readNumber};
        my $readElementName = (1 == $readsNumber) ? "lane_results" : "lane_results_r${readNumber}";
        $xmlMapRef->{$readElementName}{name} = $readElementName; 
        $xmlMapRef->{$readElementName}{title} = (1 == $readsNumber) ? 
            "Lane Results Summary" : "Lane Results Summary : Read ${readNumber}"; 
        $xmlMapRef->{$readElementName}{index} = $index++; 
        foreach my $lane ( @{ $read->{Lane} } ) {
            my $laneNumber = $lane->{laneNumber} - 1;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{align_pf_pr} = formatValueWithDeviation($lane->{percentUniquelyAlignedPF});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{alignment_score_pf} = formatValueWithDeviation($lane->{averageAlignScorePF});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{clusters_pf} = formatValueWithDeviation($lane->{clusterCountPF});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{clusters_raw} = formatValueWithDeviation($lane->{clusterCountRaw});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{error_rate_pf_pr} = formatValueWithDeviation($lane->{errorPF});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{firts_st_cycle_int_pf} = formatValueWithDeviation($lane->{oneSig});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{intensity_after_20_cycles_pf} = formatValueWithDeviation($lane->{signal20AsPctOf1});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{lane} = $laneNumber + 1;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{lane_yield_kb} = $lane->{laneYield};
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{pf_clusters_pr} = formatValueWithDeviation($lane->{percentClustersPF});
        }
    }

    foreach my $read ( @{ $summaryXmlMapRef->{ExpandedLaneSummary}{Read} } ) {
        my $readNumber = $read->{readNumber};
        my $readElementName = (1 == $readsNumber) ? "expanded_lane" : "expanded_lane_r${readNumber}";
        $xmlMapRef->{$readElementName}{name} = $readElementName; 
        $xmlMapRef->{$readElementName}{title} = (1 == $readsNumber) ? 
            "Expanded Lane Summary" : "Expanded Lane Summary : Read ${readNumber}"; 
        $xmlMapRef->{$readElementName}{index} = $index++; 
        foreach my $lane ( @{ $read->{Lane} } ) {
            my $laneNumber = $lane->{laneNumber} - 1;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{align_pf_pr} = (defined $lane->{percentUniquelyAlignedPF}{mean}) ?
                $lane->{percentUniquelyAlignedPF}{mean} : 0;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{clusters_tile_mean_raw} = sprintf "%.0f", (defined $lane->{clusterCountRaw}{mean}) ?
                $lane->{clusterCountRaw}{mean} : 0;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{cycle_10_20_av_loss_pf_pr} = formatValueWithDeviation($lane->{signalLoss10to20});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{cycle_2_10_av_loss_pf_pr} = formatValueWithDeviation($lane->{signalLoss2to10});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{cycle_2_4_av_int_pf} = formatValueWithDeviation($lane->{signalAverage2to4});
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{equiv_perfect_clusters_pf} = (defined $lane->{infoContentPF}{mean}) ? 
                $lane->{infoContentPF}{mean} : 0;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{equiv_perfect_clusters_raw} = (defined $lane->{infoContentRaw}{mean}) ?
                $lane->{infoContentRaw}{mean} : 0;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{error_rate_pf} = (defined $lane->{errorPF}{mean}) ? $lane->{errorPF}{mean} : 0;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{error_rate_raw_pr} = (defined $lane->{errorRaw}{mean}) ? $lane->{errorRaw}{mean} : 0;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{lane} = $laneNumber + 1;
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{phasing_pr} = $lane->{phasingApplied};
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{prephasing_pr} = $lane->{prephasingApplied};
            $xmlMapRef->{$readElementName}{lane}[$laneNumber]{retained_pr} = $lane->{percentClustersPF}{mean};
        }
    }

    return;
}



1;                                                   # says use was ok
__END__

