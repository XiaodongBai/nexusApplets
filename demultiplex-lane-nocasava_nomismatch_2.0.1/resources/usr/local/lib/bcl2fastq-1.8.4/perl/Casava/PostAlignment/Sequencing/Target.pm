package Casava::PostAlignment::Sequencing::Target;

# PROJECT: CASAVA
# MODULE:  $RCSfile: Target.pm,v $
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

Casava::PostAlignment::Sequencing::Target.pm - Perl library with all targets (steps of CASAVA)

=head1 SYNOPSIS

 
use Casava::PostAlignment::Sequencing::Target.pm qw();  

=head1 DESCRIPTION

Exports: 
    checkTarget($;\%;\%;\%;$)
	    
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

Created by Lukasz Szajkowski <lszajkowski@illumina.com>

=head1 SOURCE CODE

The most current release of the Perl source code 
for this module is available through CVS at genie01 a.k.a. 10.44.0.81 
cvs co BullFrog

=cut

BEGIN {
    use Exporter();
    @ISA       = qw(Exporter);
    @EXPORT    = qw();
    @EXPORT_OK =
      qw(&checkTarget &splitPEexport &splitSEexport &commitExportSplit $VERSION);
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);

use File::Basename qw(dirname);
use IO::File;
use Carp;

use Casava::Common::Log;
use Casava::Common::IOLib qw(executeCmd executeCmdEx);
use Casava::PostAlignment::Sequencing::GenomicIOLib
  qw(%exportFields %doubleExportFields $readsIdxHeader outputSEread
  outputPEread flushReadFiles parseSpliceName %predictedMateStrand @nmnmTags);
use Casava::Common::Utils qw(apply_mismatch convertToAlignmentDescriptor expandMd compressMd getChromosomeSizes);
use Casava::PostAlignment::Sequencing::Config qw(isSpliceJunctionChrom %chrEnds readProjectParameters);
use Casava::PostAlignment::Sequencing::RnaSeqLib qw($tmpSpliceCountFileName);
use Casava::PostAlignment::Sequencing::SamLib qw(matchDescToCigar);
use Casava::PostAlignment::QC::Build_Spec;
use Casava::PostAlignment::QC::Export_Dir_Spec;
use Casava::PostAlignment::QC::Lane_Spec;

#use Casava::Common::IOLib qw(bufferedPrint);

sub checkTarget($\%\%\%);
sub splitPEexport($$$$$\%\%\%);
sub splitSEexport($$$$$\%\%\%);
sub commitExportSplit($\%\%\%);
sub calculateCrc32($);

# save hash lookup time:
my $exportMachineNameIndex = $exportFields{'MachineName'};
my $exportRunNumberIndex   = $exportFields{'RunNumber'};
my $exportLaneIndex        = $exportFields{'Lane'};
my $exportTileIndex        = $exportFields{'Tile'};
my $exportXIndex           = $exportFields{'X'};
my $exportYIndex           = $exportFields{'Y'};
my $exportBarcodeIndex     = $exportFields{'Index'};
my $exportFilterIndex      = $exportFields{'Filter'};
my $exportChrIndex         = $exportFields{'Chr'};
my $exportStrandIndex      = $exportFields{'Strand'};
my $exportSingleScoreIndex = $exportFields{'SingleScore'};
my $exportPairedScoreIndex = $exportFields{'PairedScore'};
my $exportPosnIndex        = $exportFields{'Posn'};
my $exportContigIndex      = $exportFields{'Contig'};
my $exportSequenceIndex    = $exportFields{'Seq'};
my $exportQualityIxdex     = $exportFields{'QualityString'};
my $exportDescriptorIndex  = $exportFields{'Descriptor'};
my $exportPartnerChrIndex  = $exportFields{'PartnerChr'};
my $exportPartnerContigIndex  = $exportFields{'PartnerContig'};
my $exportPartnerOffsetIndex  = $exportFields{'PartnerOffset'};
my $exportPartnerStrandIndex  = $exportFields{'PartnerStrand'};

=pod

=head1 The procedure calculates a quick crc32 checksum for a given file

=over 4

=item configureTargets($confStatus, $targetRef)

The procedure calculates a checksum for a given file. In essence, this is just
used to detect if several files are different or not. A real crc32 or md5 hash
would be optimal, but this simple function does the job without adding more
module dependencies.

Parameters:
    $filename - input filename

Returns:
    nothing
=back

=cut

sub calculateCrc32($) {
    croak "ERROR: calculateCrc32\n" unless (@_ == 1);
    my ($filename) = @_;

    open(IN, $filename) or die "ERROR: Cannot open $filename: $!";
    binmode(IN);

    # initialize the crc32 table
    my $CRC_POLY = 0x04c11db;
    my $INT_MASK = 0xffffffff;

    my @crc32_table = ();
    for(my $i = 0; $i < 256; $i++) {
        my $c = $i << 24;
        for(my $j = 8; $j > 0; $j--) {
            $c = ($c & 0x80000000 ? ($c << 1) ^ $CRC_POLY : ($c << 1));
        }
        push @crc32_table, $c & $INT_MASK;
    }

    # calculate the crc32 value
    my $crc = $INT_MASK;
    while(my $numBytesRead = read(IN, my $buffer, 4096)) {
        foreach my $byte (split(//, $buffer)) {
            $crc = ($crc << 8) ^ $crc32_table[(($crc >> 24) & 0xff) ^ ord($byte)];
        }
    }

    return (~$crc & $INT_MASK);
}

=pod

=head1 The procedure check if all required by CASAVA files exists and are non zero size

=over 4

=item checkTarget

The procedure check if all required by CASAVA files exists and are non zero size.


Parameters:
    $projectDir    - CASAVA project directory 
    $CONF_APP      - HASH MAP REF with application specific configuration
    $CONF_PROJ     - HASH MAP REF with project specific configuration
    $runsConfig    - HASH MAP REF with runs specific configuration
	
Returns:
    genome sequence
        
=back

=cut

sub checkTarget($\%\%\%) {
    my ( $projectDir, $CONF_APP_REF, $CONF_PROJ_REF, $runsConfigRef ) = @_;
    my $reads_indx_f = $CONF_APP_REF->{f_reads_indx};
    
    errorExit "ERROR: No input directories were specified. Please use the --inputSampleDir option."
        unless exists $runsConfigRef->{inputDirectories};

    # check for the sample-specific metadata files
    my @genomeSizeFiles = ();

    foreach my $alignedSampleDir (@{$runsConfigRef->{inputDirectories}->{sampleDirectory}}) {

        # grab the filenames in this directory
        opendir(DIR, $alignedSampleDir);
        my @dirFiles = map { "$alignedSampleDir/$_" } readdir(DIR);
        closedir(DIR);

        # add the genome size files to our main array so that we can check
	# consistency between all genome size files.
        my @dirGenomeSizeFiles = grep(/_genomesize.xml$/, @dirFiles);
        my $numGenomeSizeFiles = scalar(@dirGenomeSizeFiles);
        push (@genomeSizeFiles, @dirGenomeSizeFiles);

        errorExit("ERROR: Could not find any genome size files (*_genomesize.xml) in $alignedSampleDir.") if ($numGenomeSizeFiles == 0);

        # check the genome size file
        my $genomesizeXml = $dirGenomeSizeFiles[0];
        if(!$CONF_PROJ_REF->{ignoreMissingMetadata}) {
            my $chromSizesRef = getChromosomeSizes($genomesizeXml,
                {contigName=>'contigName', fileName=>'fileName'}->{$CONF_PROJ_REF->{chromNameSource}});

            for my $chrom (keys %{$chromSizesRef}) {
 
                errorExit "ERROR: Chromosome $chrom found in $genomesizeXml is not present in $CONF_PROJ_REF->{dirRefSeq}"
                    unless (exists $chrEnds{$chrom});

                errorExit "ERROR: Chromosome $chrom size found in $genomesizeXml is different than expected. "
                         ."Expected: $chrEnds{$chrom}->{length}. Actual: $chromSizesRef->{$chrom}. "
                         ."Reference: $CONF_PROJ_REF->{dirRefSeq}"
                    if ( $chromSizesRef->{$chrom} != $chrEnds{$chrom}->{length} );
            }
        }
    }

    # check the consistency of all genome size files
    my @genomeSizeChecksums = ();
    foreach my $genomeSizeFilename (@genomeSizeFiles) {
	my $crc = calculateCrc32($genomeSizeFilename);
	push (@genomeSizeChecksums, $crc);
    }

    my %genomeSizeChecksumSet = map { $_ => 1 } @genomeSizeChecksums;
    my $numDistinctChecksums = scalar keys(%genomeSizeChecksumSet);

    errorExit("ERROR: Expected to find identical genome size files (*_genomesize.xml) in all aligned sample directories, but found $numDistinctChecksums different genome size files.") unless ($numDistinctChecksums == 1);

    # check all the export files
    foreach my $exportSet (@{$runsConfigRef->{exportFiles}}) {

        # check the read 1 export files
        if($exportSet->{Read1}) {
            errorExit "ERROR: Could not find the read 1 export file (" . $exportSet->{Read1} . ")."
                unless (-s $exportSet->{Read1});
        }
        
        # check the read 2 export files
        if($exportSet->{Read2}) {
            errorExit "ERROR: Could not find the read 2 export file (" . $exportSet->{Read2} . ")."
                unless (-s $exportSet->{Read2});
        }
        
        # check for paired.xml if we have both read 1 and read 2 export files
        if('paired' eq $CONF_PROJ_REF->{readMode}) {
            errorExit "ERROR: Could not find the pair.xml file associated with " . $exportSet->{Read1} . " and " . $exportSet->{Read2} . "."
                unless (-s $exportSet->{PairXml});
        }

    }
}

=pod

=head1 Produces hash of filtered tiles from build file.

=over 4

=item getFilteredTileHash($buildFilePath)

Parameters:
    $buildFilePath    - build file path
    $exportFileDir    - export directory (GERALD directory)
    $laneNum          - lane number
Returns:
    Hash mapping filtered tiles to 1; mappings for other tiles do not exist.

=back

=cut

sub getFilteredTileHash($;$;$)
{
    my ($buildFilePath, $exportFileDir, $laneNum) = @_;

    if (!$buildFilePath) {
	return ();
    }

    my $buildSpec = Casava::PostAlignment::QC::Build_Spec->new($buildFilePath);

    my $exportDirSpec = $buildSpec->get_export_dir_spec($exportFileDir);
    my $laneSpec = $exportDirSpec->get_lane_spec($laneNum);
    my @badTiles = $laneSpec->get_bad_tiles();

    # The `int' is needed to strip off leading zeros for comparison with
    # the tile column in an export file.
    my %isFilteredTile = map {int($_) => 1} @badTiles;
    return %isFilteredTile;
}



# part of the annoation is shared by both reads:
#
sub markShadowReadCommon(\@) {
    my ( $ref ) = @_;
    $ref->[$exportPairedScoreIndex]   = 0;
    $ref->[$exportPartnerChrIndex]    = '';
    $ref->[$exportPartnerContigIndex] = '';
    $ref->[$exportPartnerOffsetIndex] = 0;
    $ref->[$exportPartnerStrandIndex] = 'N';
}



#
# makes field changes required to indicate a mapped/shadow read pair
#
sub markShadowReadPair(\@;\@) {
    my ( $mappedReadRef , $shadowReadRef ) = @_;

    $shadowReadRef->[$exportChrIndex]           = "NM";
    $shadowReadRef->[$exportPosnIndex]          = $mappedReadRef->[$exportPosnIndex];
    $shadowReadRef->[$exportStrandIndex]        = $predictedMateStrand{$mappedReadRef->[$exportStrandIndex]};
    $shadowReadRef->[$exportDescriptorIndex]    = '-';
    $shadowReadRef->[$exportSingleScoreIndex]   = -1;
    markShadowReadCommon(@$shadowReadRef);
    markShadowReadCommon(@$mappedReadRef);
}



#
sub unMatchRead(\@) {
    my ( $ref ) = @_;
    $ref->[$exportChrIndex] = "NM";
    for ($exportContigIndex..$exportPartnerStrandIndex) { $ref->[$_] = ''; }
}



#
# consolidate repeated index open code:
#
sub getIDXFH($) {
    my ( $readsIndexFile ) = @_;

    # relax restriction against rerunning the 'sort' target:
    #
    # if ( -e $readsIndexFile ) {
    #     errorExit "ERROR: sortExport $readsIndexFile "
    #       . "already added to this build. "
    #       . "Each lane can be processed only once.";
    # }

    my $IDXFH;
    open($IDXFH, ">$readsIndexFile" )
      or errorExit "ERROR: Could not find the index file "
        . "$readsIndexFile for writing $!\n";
    return $IDXFH;
}



sub getUnmapFH($;$) {
    my ( $CONF_APP_REF, $NMNM_DIR ) = @_;

    my %unmapFH;
    for my $key (@nmnmTags) {
        my $label = $CONF_APP_REF->{$key};
        my $path = File::Spec->catfile($NMNM_DIR,  $label .".txt");
        $unmapFH{$key} = IO::File->new(">$path")
          or errorExit "ERROR: Couldn't create file handle for $path $!\n";
    }
    return \%unmapFH;
}



=pod

=head1 The procedure splits pair-ended PE export files into chromosomes and bins

=over 4

=item splitPEexport

The procedure splits pair-ended export files into chromosomes and bins


Parameters:
    $projectDir          - CASAVA project directory 
    $buildFilePath       - build file path
    $read1ExportFilename - read 1 export filename
    $read2ExportFilename - read 2 export filename
    $exportSet           - the export set ID
    $CONF_PROJ_REF       - HASH MAP REF with application specific configuration
    $CONF_APP_REF        - HASH MAP REF with project specific configuration
    $runsConfigRef       - HASH MAP REF with runs specific configuration	
Returns:
    genome sequence
        
=back

=cut

sub splitPEexport($$$$$\%\%\%) {
    my ($projectDir, $buildFilePath, $read1ExportFilename,
        $read2ExportFilename, $exportSet, $CONF_APP_REF, $CONF_PROJ_REF,
        $runsConfigRef) = @_;

    my $totalReads     = 0;
    my $totalBases     = 0;
    my $usedReads      = 0;
    my $usedBases      = 0;
    my $goodReads      = 0;
    my $goodReadBases  = 0;
    my $originalReads  = 0;
    my $includedReads  = 0;
    my $failedTileFilt = 0;
    my $failedFilt     = 0;
    my $failedQC       = 0;
    my $unanchored     = 0;
    my $nonUniquePairs = 0;
    my $nmnmPairs      = 0;
    my $riboCount      = 0;
    my $goodPairs      = 0;
    my $mixedPairs     = 0;
    my $controlPairs   = 0;
    my $bufferSize     = 0;
    my %OUTPUTFILES;
    my $NMNM_TAG        = $CONF_APP_REF->{TAG_NMNM};
    my $timeStampFormat = $CONF_APP_REF->{formatTimeStamp};
    my $reads_indx_f    = $CONF_APP_REF->{f_reads_indx};
    my $binSize         = $CONF_PROJ_REF->{binSizeProject};
    my $toNMScore       = (defined $CONF_PROJ_REF->{toNMScore}) ? $CONF_PROJ_REF->{toNMScore} : -1;
    my $readSampleRateInput = $CONF_PROJ_REF->{readSampleRateInput};
    my $resultsDir = File::Spec->catdir( $CONF_PROJ_REF->{dirBuildExport}, 'sets', $exportSet );
    my $isKeepAllReads = ((defined $CONF_PROJ_REF->{sortKeepAllReads}) and $CONF_PROJ_REF->{sortKeepAllReads});
    my $isCompressPair = (not ((defined $CONF_PROJ_REF->{sortNoCompressPair}) and ($CONF_PROJ_REF->{sortNoCompressPair})));
    my $isKeepFilteredReads = $isKeepAllReads;
    my $isKeepUnmappedReads = $isKeepAllReads;

    my $NMNM_DIR  = File::Spec->catdir( $resultsDir, $NMNM_TAG );
    my %Fields    = %exportFields;
    my $fieldSize = scalar( keys %Fields );

    # TODO: Make it possible to get the filtered tile
    #my %isFilteredTile
    #= getFilteredTileHash($buildFilePath, $exportFileDir, $lane);

    my $IDXFH = getIDXFH( File::Spec->catfile( $resultsDir, $reads_indx_f ) );

    my $confDir = $CONF_APP_REF->{dirConf};
    my $dirsPath = File::Spec->catfile( $projectDir , $confDir , "project.dirs.tar" );
    my $cmd = "tar -C $resultsDir/ -xf $dirsPath";
    executeCmd( $cmd, 5 );

    my $unmap;
    $unmap = getUnmapFH($CONF_APP_REF,$NMNM_DIR) if($isKeepUnmappedReads);

    my $compFH;
    if($isCompressPair) {
        my $libexecDir = '/usr/local/libexec/bcl2fastq-1.8.4';
        my $cxp = File::Spec->catfile($libexecDir,'compressXPair');
        my $cmd = "| $cxp -s '$resultsDir'";
        open($compFH,$cmd) or errorExit("ERROR: failed to open process '$cmd' $!\n");
    }
    
    # grab the fragment length statistics from the pair.xml file
    my $pairXmlFilename = $read1ExportFilename;
    $pairXmlFilename =~ s/.gz$//;
    $pairXmlFilename =~ s/_export.txt/_pair.xml/;
    $pairXmlFilename =~ s/_R1_/_/;

    errorExit("ERROR: Unable to find the pair.xml file (${pairXmlFilename}).\n") unless (-s $pairXmlFilename);
            
    # open the XML file
    my $xs = new XML::Simple(
        KeyAttr    => [],
    );
            
    my $insertSizeRef = $xs->XMLin($pairXmlFilename);
        
    # grab the fragment length data 
    my $highSDFragmentLength = $insertSizeRef->{InsertSize}->{HighSD};
    my $lowSDFragmentLength  = $insertSizeRef->{InsertSize}->{LowSD};
    my $maxFragmentLength    = $insertSizeRef->{InsertSize}->{Max};
    my $medianFragmentLength = $insertSizeRef->{InsertSize}->{Median};
    my $minFragmentLength    = $insertSizeRef->{InsertSize}->{Min};
    my $nominalOrientation   = $insertSizeRef->{Orientation}->{Nominal};
    
    #print STDERR "\n";
    printLog( "Looking at $read1ExportFilename & $read2ExportFilename\n", 5 );
    my $exportFile1 = $read1ExportFilename;
    my $exportFile2 = $read2ExportFilename;
    open( EXPT_1, "gunzip -fc $exportFile1 |" )
      || errorExit "ERROR: Could not open export file $exportFile1 $!\n";
    open( EXPT_2, "gunzip -fc $exportFile2 |" )
      || errorExit "ERROR: Could not open export file $exportFile2 $!\n";
    my %exportStorage = ();
    $exportStorage{ref2Output}     = \%OUTPUTFILES;
    $exportStorage{resultsDir}     = $resultsDir;
    $exportStorage{binSize}        = $binSize;
    $exportStorage{isCompressPair} = $isCompressPair;
    $exportStorage{compFH}         = $compFH;
    my $needRunIdFix = 0;
    
    my $needRunConfigData = 1;
    my $instrumentName;
    my $runID;
    my $laneNumber;
    my $barcodeSequence;
        
    #my $machineName  = $runsConfigRef->{run}{$runId}{machine};
    #my $exportRunId  = $runsConfigRef->{run}{$runId}{exportRunId};

    #if ( $runsConfigRef->{run}{$runId}{fixId} ne 'none' ) {
    #    $needRunIdFix = 1;
    #}
    
    my ($read1SeqLength, $read2SeqLength) = (0,0);
    while (<EXPT_1>) {
        chomp;
        my @read1 = split /\t/, $_;    # read1 array
        my $read2Str = <EXPT_2>;

	# retrieve the instrument name, run ID, lane, and barcode
	if($needRunConfigData) {
	    $instrumentName  = $read1[$exportMachineNameIndex];
	    $runID           = $read1[$exportRunNumberIndex];
	    $laneNumber      = $read1[$exportLaneIndex];
	    $barcodeSequence = $read1[$exportBarcodeIndex];
	    $needRunConfigData = 0;
	}

        ++$originalReads;
        next if($readSampleRateInput && $originalReads * $readSampleRateInput <= $includedReads);
        ++$includedReads;

        chomp($read2Str);
        my @read2 = split /\t/, $read2Str;    # read2 array
        if ( scalar(@read1) != $fieldSize ) {
            errorExit "ERROR: parsing $exportFile1 file at $originalReads (wrong number of fields "
              . scalar(@read1) . " != " . $fieldSize . ")\n";
        }

        if ( scalar(@read2) != $fieldSize ) {
            errorExit "ERROR: parsing $exportFile2 file at $originalReads (wrong number of fields "
              . scalar(@read2) . " != " . $fieldSize . ")\n";
        }

        $totalReads += 2;
        $read1SeqLength = length($read1[$exportSequenceIndex]);
        $read2SeqLength = length($read2[$exportSequenceIndex]);
        $totalBases += $read1SeqLength + $read2SeqLength;


        # Ignore tiles specified as filtered out.
        # TODO: Make it possible to get the filtered tile
        #if (exists($isFilteredTile{$read1[$exportTileIndex]})) {
        #    ++$failedTileFilt;
        #    next;
        #}

        #if ( $needRunIdFix == 1 ) {
        #    ## those reads need to be fixed (machine name and run id)
        #    $read1[ $exportMachineNameIndex ] = $machineName;
        #    $read2[ $exportMachineNameIndex ] = $machineName;
        #    $read1[ $exportRunNumberIndex ]   = $exportRunId;
        #    $read2[ $exportRunNumberIndex ]   = $exportRunId;
        #}

        # #
        # # assert read-pair clusters, offset and PE score:
        # #
        # if( ($read1[ $exportMachineNameIndex ] ne $read2[ $exportMachineNameIndex ]) or
        #     ($read1[ $exportRunNumberIndex ] ne $read2[ $exportRunNumberIndex ]) or
        #     ($read1[ $exportLaneIndex ] ne $read2[ $exportLaneIndex ]) or
        #     ($read1[ $exportTileIndex ] ne $read2[ $exportTileIndex ]) or
        #     ($read1[ $exportXIndex ] ne $read2[ $exportXIndex ]) or
        #     ($read1[ $exportYIndex ] ne $read2[ $exportYIndex ])) {
        #     errorExit "ERROR: Mate pair reads with different clusters. Line Number: $originalReads in pair export files:\n".
        #               "(1) $exportFile1\n".
        #               "(2) $exportFile2\n";
        # }

        # my $PEVal1 = $read1[ $exportPairedScoreIndex ];
        # my $PEVal2 = $read2[ $exportPairedScoreIndex ];
        # my $chrVal1    = $read1[$exportChrIndex];
        # my $chrVal2    = $read2[$exportChrIndex];
        # my $POffsetVal1 = $read1[ $exportPartnerOffsetIndex ];
        # my $POffsetVal2 = $read2[ $exportPartnerOffsetIndex ];

        # if(($PEVal1 ne "") and ($PEVal2 ne "") and ($PEVal1 ne $PEVal2)) { 
        #     errorExit "ERROR: Mate pair reads with different paired-end score. Line Number: $originalReads in pair export files:\n".
        #               "(1) $exportFile1\n".
        #               "(2) $exportFile2\n";
        # }

        # if(($POffsetVal1 ne "") and ($POffsetVal2 ne "") and ($chrVal1 eq $chrVal2) and ($POffsetVal1 != -$POffsetVal2)) { 
        #     errorExit "ERROR: Mate pair reads with incompatible paired-offset score. Line Number: $originalReads in pair export files:\n".
        #               "(1) $exportFile1\n".
        #               "(2) $exportFile2\n";
        # }


        # make adjustment for SAM:
        if( ($read1[ $exportPosnIndex ] ne '') and
            ($read1[ $exportPosnIndex ] < 1))
        {
            if( ($read2[ $exportPosnIndex] eq '') or
                ($read2[ $exportPosnIndex ] < 1))
            {   # not a shadow -- both reads are bad:
                unMatchRead(@read1);
                unMatchRead(@read2) if($read2[ $exportPosnIndex ] ne '');
            }
            else
            {  # read1 shadow case:
                markShadowReadPair(@read2,@read1);
            }
        }
        elsif(($read2[ $exportPosnIndex ] ne '') and
              ($read2[ $exportPosnIndex ] < 1))
        { 
            if( ($read1[ $exportPosnIndex] eq '') or
                ($read1[ $exportPosnIndex ] < 1))
            {   # not a shadow -- both reads are bad:
                unMatchRead(@read2);
                unMatchRead(@read1) if($read1[ $exportPosnIndex ] ne '');
            }
            else
            {  # read2 shadow case:
                markShadowReadPair(@read1,@read2);
            }
        }

        my $filterVal1 = $read1[$exportFilterIndex] eq "N";
        my $filterVal2 = $read2[$exportFilterIndex] eq "N";

        errorExit ("ERROR: Ends of a pair have different filter values ("
            . join (" ", @read1, "...", join " ", @read2)
            . ")\n Run aborted \n") unless ($filterVal1 eq $filterVal2);

        my $read1Sras = $read1[ $exportSingleScoreIndex ];
        $read1Sras = 0 if ('' eq $read1Sras);
        my $read2Sras = $read2[ $exportSingleScoreIndex ];
        $read2Sras = 0 if ('' eq $read2Sras);
        my $pras = $read1[ $exportPairedScoreIndex ];
        $pras = 0 if ('' eq $pras);
        my $pras2 = $read2[ $exportPairedScoreIndex ];
        $pras2 = 0 if ('' eq $pras2);

        errorExit ("ERROR: Ends of a pair have different pair alignment scores ("
            . join (" ", @read1, "...", join " ", @read2)
            . ")\n Run aborted \n") unless ($pras eq $pras2);

        ## Re-classify anomalous reads
        if (   ($read1[ $exportPosnIndex ] ne "")
            && ($read2[ $exportPosnIndex ] ne "")
            && (0 == $pras) )
        {
	    	# read 1 shadow, read 2 singleton ?
            if (   ( $read1Sras < $toNMScore )
                && ( $read2Sras >= $toNMScore ) )
            {
                markShadowReadPair(@read2,@read1);
            }    # if

	    	# read 2 shadow, read 1 singleton ?
            elsif (   ( $read2Sras < $toNMScore )
                   && ( $read1Sras >= $toNMScore ) )
            {
                markShadowReadPair(@read1,@read2);
            }    # if
        }

        my $strandVal1 = $read1[$exportStrandIndex];
        my $strandVal2 = $read2[$exportStrandIndex];
        my $strandPosn1 = $read1[$exportPosnIndex];
        my $strandPosn2 = $read2[$exportPosnIndex];
        my $chrVal1    = $read1[$exportChrIndex];
        my $chrVal2    = $read2[$exportChrIndex];

        if ( $filterVal1 )
        {   # enumerate any pairs that fail the filters
            $failedFilt++;
            next unless($isKeepFilteredReads);
        }

        my $unmapTag;

        if ( ( $chrVal1 eq "QC" ) && ( $chrVal2 eq "QC" ) )
        {   # enumerate any reads where neither reads passes QC
            $failedQC++;
            $unmapTag="TAG_qcFail";
        }
        elsif ( ( $chrVal1 =~ /\d+:\d+:\d+/ )
         	&& ( $chrVal2 =~ /\d+:\d+:\d+/ ) )
        {   # enumerate any reads where neither reads aligns uniquely
            $nonUniquePairs++;
            $unmapTag="TAG_nonUnique";
        }
        elsif ( ( $chrVal1 eq "NM" ) && ( $chrVal2 eq "NM" ) )
        {
            # grab the NM-NM matches and put them to one side
            # in case they are needed for de novo
            $nmnmPairs++;
            $unmapTag="TAG_noMatch";
        }
        elsif ( ($chrVal1 eq "RM") and ($chrVal2 eq "RM") )
        {
            $riboCount++;
            $unmapTag="TAG_rm";
        }
        elsif ( ($chrVal1 eq "CONTROL") and ($chrVal2 eq "CONTROL") )
        {
            $controlPairs++;
            $unmapTag="TAG_control";
        }
        elsif (($chrVal1 =~ /^(NM|QC|RM|CONTROL|\d+:\d+:\d+)$/) and
               ($chrVal2 =~ /^(NM|QC|RM|CONTROL|\d+:\d+:\d+)$/))
        {
            $mixedPairs++;
            $unmapTag="TAG_mixed";
        }

        if(defined $unmapTag)
        {
            if($isKeepUnmappedReads) {
                my $unmapFH = $unmap->{$unmapTag};
                print $unmapFH join("\t", @read1),"\n", join("\t", @read2),"\n";
            }
            next;
        }

        if ( $read1Sras == 0 && $read2Sras == 0 )
        {   # all alignment scores are zero
            $unanchored++;
        }

        ## look at the normal pairs ##
        if ( ( $pras =~ /^\d+$/ ) && ( $pras > 0 ) )
        {
            # if we have a paired read score
            if ( $chrVal1 eq $chrVal2 )
            {
                # both reads on the same chromosome
                # ignores contigs currently
                # would have to check that both reads are on same contig
                my $chrom = $chrVal1;
                if ( ( $strandVal1 eq "F" ) && ( $strandVal2 eq "R" ) )
                {
                    my $start = $strandPosn1;
                    my $end = $strandPosn2 + $read2SeqLength;
                    outputPEread( %exportStorage, $start, $end, $chrom, @read1, @read2, "norm" );
                }
                elsif ( ( $strandVal1 eq "R" ) && ( $strandVal2 eq "F" ) )
                {
                    my $start = $strandPosn2;
                    my $end = $strandPosn1 + $read1SeqLength;
                    outputPEread( %exportStorage, $start, $end, $chrom, @read1, @read2, "norm" );
                }
                else {
                    errorExit ("ERROR: Non standard entry (". join ("\t", @read1) . ")\n Run aborted \n");
                }

                ++$goodPairs if (($pras > $CONF_PROJ_REF->{QVCutoff}) and (not $filterVal1));
            }
            else {
                errorExit "ERROR: splitPEexport: PairedScore > 0 when chr_read1 != chr_read2 \n";
            }
        }
        ## look at the anomalous matches
        else {
            # if only one read maps
            # put it in the orphan file
            if ( ( $chrVal1 =~ /^(NM|QC|RM|CONTROL|\d+:\d+:\d+)$/ ) )
            {   # if read 1 is bad
                my $start = $strandPosn2;
                my $end   = $strandPosn2 + $read2SeqLength;
                my $chrom = $chrVal2;
                outputPEread( %exportStorage, $start, $end, $chrom, @read2, @read1, "orphan" );
            }
            elsif ( ( $chrVal2 =~ /^(NM|QC|RM|CONTROL|\d+:\d+:\d+)$/ ) )
            {   # if read2 is bad
                my $start = $strandPosn1;
                my $end   = $strandPosn1 + $read1SeqLength;
                my $chrom = $chrVal1;
                outputPEread( %exportStorage, $start, $end, $chrom, @read1, @read2, "orphan" );
            }
            elsif (($chrVal2) && 1 == ( $chrVal1 cmp $chrVal2 ) )
            {        # cross chromosome match
                     # put it in the anom file to be sorted
                     # store everything by chrVal1
                my $start = $strandPosn1;
                errorExit ("Invalid position:\n" . join(' ', @read1) . '\n' .  join(' ', @read2)) 
                    unless ( $strandPosn1 =~ /^\d+$|^-\d+$/ ) && ( $strandPosn2 =~ /^\d+$|^-\d+$/ );
                my $end = ( $strandVal1 eq "F" )
                  ? $strandPosn2 + $read2SeqLength
                  : $strandPosn2;
                $end = "$chrVal2:$end";
                outputPEread( %exportStorage, $start, $end, $chrVal1, @read1, @read2, "anom" );
            }
            elsif (($chrVal2) && 1 == ( $chrVal2 cmp $chrVal1 ) )
            {        # cross chromosome match
                     # put it in the anom file to be sorted
                     # store everything by chrVal2
                my $start = $strandPosn2;
                errorExit ("Invalid position:\n" . join(' ', @read1) . '\n' .  join(' ', @read2)) 
                    unless ( $strandPosn1 =~ /^\d+$|^-\d+$/ ) && ( $strandPosn2 =~ /^\d+$|^-\d+$/ );
                my $end = ( $strandVal2 eq "F" )
                  ? $strandPosn1 + $read1SeqLength
                  : $strandPosn1;
                $end = "$chrVal1:$end";
                outputPEread( %exportStorage, $start, $end, $chrVal2, @read2, @read1, "anom" );
            }
            elsif ( ( $strandVal1 eq "F" ) && ( $strandVal2 eq "R" ) )
            {
                # oversized pairs to go into anomaly file
                my $start = $strandPosn1;
                errorExit ("Invalid position:\n" . join(' ', @read1) . '\n' .  join(' ', @read2)) 
                    unless ( $start =~ /^-?\d+$/ );
                my $end = $strandPosn2 + $read2SeqLength;
                outputPEread( %exportStorage, $start, $end, $chrVal1, @read1, @read2, "anom" );
            }
            elsif ( ( $strandVal1 eq "R" ) && ( $strandVal2 eq "F" ) )
            {            # oversized
                         # to go into anomaly file
                my $start = $strandPosn2;
                errorExit ("Invalid position:\n" . join(' ', @read1) . '\n' .  join(' ', @read2)) 
                    unless ( $start =~ /^-?\d+$/ );
                my $end   = $strandPosn1 + $read1SeqLength;
                outputPEread( %exportStorage, $start, $end, $chrVal1, @read1, @read2, "anom" );
            }
            elsif ( ($strandVal1 eq $strandVal2) &&
                    ('R' eq $strandVal1 || 'F' eq $strandVal1) )
            {  # Fwd-Fwd OR Rev-Rev anomalous match to go into anomaly file

                my ($start, $end) = ( $strandPosn1 < $strandPosn2 )
                  ? ($strandPosn1, $strandPosn2 + $read2SeqLength )
                  : ($strandPosn2, $strandPosn1 + $read1SeqLength );

                outputPEread( %exportStorage, $start, $end, $chrVal1, @read1, @read2, "anom" );
            }
            else {
                errorExit "ERROR: Unexpected read record content:\n" . join (" ", @read1) . "\n" . join(" ", @read2);
            }
        }
        if(not $filterVal1){
            ++$usedReads;
            $usedBases += $read1SeqLength + $read2SeqLength;
            if ($read1Sras > $CONF_PROJ_REF->{QVCutoffSingle})
            {
                ++$goodReads;
                $goodReadBases += $read1SeqLength;
            }
            if ($read2Sras > $CONF_PROJ_REF->{QVCutoffSingle})
            {
                ++$goodReads;
                $goodReadBases += $read2SeqLength;
            }
        }
    }
    close(EXPT_1);
    close(EXPT_2);
    if($isKeepUnmappedReads) {
        for(values %$unmap) { close }
    }
    if($isCompressPair) {
        close($compFH) or errorExit("ERROR: failed to close process: $!\n");
    } else {
        flushReadFiles(%OUTPUTFILES, 1);
    }

    my $date = strftime $timeStampFormat, localtime;
    print $IDXFH "$readsIdxHeader\n";
    print $IDXFH "$totalReads"
      . "\t$totalBases"
      . "\t$usedReads"
      . "\t$usedBases"
      . "\t$goodReads"
      . "\t$goodReadBases"
      . "\t$failedTileFilt"
      . "\t$failedFilt"
      . "\t$failedQC"
      . "\t$unanchored"
      . "\t$nonUniquePairs"
      . "\t$goodPairs"
      . "\t$mixedPairs"
      . "\t$riboCount"
      . "\t0" #$mitoCount"
      . "\t0" #$splicedReads"
      . "\t" . ($CONF_PROJ_REF->{skipVariableMetadata} ? '' :  dirname($exportFile1))
      . "\t$exportSet" #$lane"
      . "\t$date"
      . "\t$instrumentName"
      . "\t$runID"
      . "\t$laneNumber"
      . "\t$barcodeSequence"
      . "\t" . (defined $medianFragmentLength ? $highSDFragmentLength : '')
      . "\t" . (defined $medianFragmentLength ? $lowSDFragmentLength : '')
      . "\t" . (defined $medianFragmentLength ? $maxFragmentLength : '')
      . "\t" . (defined $medianFragmentLength ? $medianFragmentLength : '')
      . "\t" . (defined $medianFragmentLength ? $minFragmentLength : '')
      . "\t" . (defined $medianFragmentLength ? $nominalOrientation : '')
      . "\t$read1SeqLength"
      . "\t$read2SeqLength"
      . "\n" unless $needRunConfigData; #don't print the line if we did not parse a single read

    close($IDXFH);
}



#
# convert export MD to alignment
#
sub getAlignment($) {
    my ($read) = @_;

    my $cigar = matchDescToCigar( $read->[$exportDescriptorIndex] );
    if($cigar !~ /^(\d+[MID])*$/) {
        errorExit("ERROR: unexected cigar string '$cigar', generated for read: ".join("\t",@$read)."\n");
    }

    my @align;
    push @align, [$1,$2] while($cigar =~ /(\d+)([MID])/g);

    # eliminate any edge deletions:
    while(scalar(@align) and ($align[0][1] eq "D")) {shift @align;}
    while(scalar(@align) and ($align[-1][1] eq "D")) {pop @align;}

    @align = reverse(@align) if($read->[$exportStrandIndex] eq 'R');
    return @align;
}



#
# trim any part of alignment below position 1 and return new start
# position
#
sub trimAlignment($;$) {
    my ($alignRef,$startPos) = @_;

    my @align2;
    my $pos = $startPos;
    my ($i,$readPos) = (0,0);
    for (@$alignRef){
        my ($length,$type) = @$_;
        if(($type eq 'M') and
           ($pos+$length >= 1)) {
            my $left=$readPos+(1-$pos);
            push @align2, [$left,'S'] if($left);
            push @align2, [($length-$left),$type] if($length-$left);
            $startPos=1;
            last;
        } elsif(($type eq 'D') and
                (($pos+$length) >= 1)) {
            push @align2, [$readPos,'S'] if($readPos);
            $startPos=$pos+$length;
            last;
        }
        $pos+=$length if($type =~ /[MD]/);
        $readPos+=$length if($type =~ /[MI]/);
        $i++;
    }
    splice(@$alignRef,0,$i+1,@align2) if(@align2);
    return $startPos;
}



=pod

=head1 outputSpliceJunctionReadGap($exportStorage, $read1, $inputFileLine)

The procedure outputs SE splice junction-aligned read as a single
genome-aligned read with an invalid MD string, '-', and a BAM cigar
string in export field 23. In case of an off-junction read or
exotic/non-colinear splice site, the function returns a nonzero error
value without writing the read.

=over 4

=item #

Bases that are overhanging at either end of the splice site are marked
as soft clipped in BAM cigar

=item #

If read does not cover the splice junction by at least one base, a
warning is emitted and no read is output by this function. An error
code is returned in this case.

Rationale: this means a bug in pickBestAlignmentRNA and it should be
fixed there. For archival bin/sort, it is expected that these error
reads will be written as unmapped.

=item #

If an indel occurs adjacent to or across the splice junction no read
is output and an error code is returned.

Rationale: Allowing alternate splice sites would be a form a splice
site discovery which would make the current splice and exon statistics
inaccurate.

=back

Parameters:
    $exportStorage  - output storage
    $read           - Refernce to array containing components of a read line 
    $inputLineNumber  Line location in input file. Used for diagnostics output

Returns:
    zero if the read covers junction point

=cut

sub outputSpliceJunctionReadGap(\%\@$\%$)
{
    use constant {
        SUCCESS => 0,
        ERROR_NON_COLINEAR => 1,
        ERROR_OFF_JUNCTION => 2,
        ERROR_JUNCTION_INDEL => 3,
    };

    my ( $exportStorage, $read, $inputLineNumber, $spliceJunctionCount, $isCountable) = @_;
    my %spliceJunctionFeature = ();
    my $contig = $read->[$exportContigIndex];
    parseSpliceName( $contig, %spliceJunctionFeature );

    # sleft = number of bases in the splice read to the left of the splice site
    # sright = number of bases in the splice read to the right of the splice site
    # eleft = 1-indexed location of the last base of the left exon
    # eright = 1-indexed location of the first base of the right exon
    #
    my ( $spliceLeftSize, $spliceRightSize, $chrom, $exonLeftEnd, $exonRightStart ) = 
        ($spliceJunctionFeature{length1}, $spliceJunctionFeature{length2}, 
         $spliceJunctionFeature{scaffold}, 
         $spliceJunctionFeature{start}, $spliceJunctionFeature{end});

    unless($exonLeftEnd < $exonRightStart) {
        logWarning "Skipping non-colinear splice site: $exportContigIndex";
        return ERROR_NON_COLINEAR;
    }

    my @align = getAlignment($read);
    my $startPos = $read->[$exportPosnIndex];
    my $endPos= $startPos-1;
    for(@align) { $endPos += $_->[0] if($_->[1] =~ /[MD]/) };

    if(($startPos > $spliceLeftSize) or ($endPos < $spliceLeftSize)) {
        logWarning "Splice read does not cross splice junction: " . join ("\t", @$read);
        return ERROR_OFF_JUNCTION;
    }

    # soft clip segments of read hanging off either end of the splice read:
    #
    if($startPos<1){
        $startPos = trimAlignment(\@align,$startPos);
    }

    my $spliceSize=($spliceLeftSize+$spliceRightSize);
    if($endPos>$spliceSize){
        @align = reverse(@align);
        trimAlignment(\@align,(1+$spliceSize-$endPos));
        @align = reverse(@align);
#        my $str; $str .= $_->[0] . $_->[1] for(@align);
#        print STDERR "END SPLICE: $spliceSize $endPos $str\n";
    }

    # find genomic start coordinate and spliced CIGAR.
    #
    my ($startPosGenome, $splicedCigar);
    {
        my $splicePos = $startPos;
        my (@align2, $index2);
        for my $i (0..$#align){
            my ($length,$type) = @{$align[$i]};
            if(($type eq 'M') and
               ($splicePos <= $spliceLeftSize) and
               (($splicePos+$length) > $spliceLeftSize)) {
                my $left = 1+$spliceLeftSize-$splicePos;
                $index2=$i;
                push @align2, [ $left, $type ] if($left);
                push @align2, [ ($exonRightStart-$exonLeftEnd)-1 , 'N' ];
                push @align2, [ ($length-$left), $type ] if($length-$left);
                $startPosGenome=1+$exonLeftEnd+($startPos-($splicePos+$left));
            }
            if((($type eq 'D') and
                ($splicePos <= ($spliceLeftSize+1)) and
                (($splicePos+$length) >= $spliceLeftSize)) or
               (($type eq 'I') and
                ($splicePos == ($spliceLeftSize+1)))) {
                # This is not an error condition so much as an
                # expression of ELAND-RNA's restriction against splice
                # site discovery, so leave this warning out of
                # production version:
                #
                # logWarning "Splice read contains indel at junction point. " . join ("\t", @$read);
                return ERROR_JUNCTION_INDEL;
            }
            $splicePos += $length if($type =~ /[MD]/);
        }
        splice(@align,$index2,1,@align2);
        $splicedCigar .= $_->[0] . $_->[1] for(@align);
    }

    if((not $startPosGenome) or (not $splicedCigar)){
        errorExit("ERROR: failed to splice read: ".join("\t",@$read)."\n");
    }

    $read->[$exportChrIndex] = $chrom;
    $read->[$exportContigIndex] = '';
    $read->[$exportDescriptorIndex]  = '-';
    $read->[$exportPosnIndex]  = $startPosGenome;
    push @$read, "CIGAR:".$splicedCigar;

    outputSEread( %$exportStorage, $startPosGenome, $chrom, @$read, "single" );
    $spliceJunctionCount->{$chrom}{$contig}++ if($isCountable);

    return SUCCESS;
}



=pod

=head1 The procedure splits single-ended SE export files into chromosomes and bins

=over 4

=item splitSEexport

The procedure splits single-ended SE export files into chromosomes and bins


Parameters:
    $projectDir          - CASAVA project directory 
    $buildFilePath       - build file path
    $read1ExportFilename - read 1 export filename
    $exportSet           - export set ID
    $CONF_PROJ_REF       - HASH MAP REF with application specific configuration
    $CONF_APP_REF        - HASH MAP REF with project specific configuration
    $runsConfigRef       - HASH MAP REF with runs specific configuration	
	
Returns:
    genome sequence
        
=back

=cut

sub splitSEexport($$$$$\%\%\%) {

    my ($projectDir, $buildFilePath, $read1ExportFilename, $read2ExportFilename, $exportSet,
        $CONF_APP_REF, $CONF_PROJ_REF, $runsConfigRef) = @_;

    my $totalReads       = 0;
    my $totalBases       = 0;
    my $usedReads        = 0;
    my $usedBases        = 0;
    my $goodReads        = 0;
    my $goodReadBases    = 0;
    my $originalReads       = 0;
    my $includedReads       = 0;
    my $failedTileFilt      = 0;
    my $failedFilt          = 0;
    my $mitoCount           = 0;
    my $riboCount           = 0;
    my $controlCount        = 0;
    my $splicedReads = 0;
    my $failedQC            = 0;
    my $singleExclude       = 0;
    my $nonUniqueAlign      = 0;
    my $nonMappedAlign      = 0;
    my %OUTPUTFILES;
    my $NMNM_TAG        = $CONF_APP_REF->{TAG_NMNM};
    my $timeStampFormat = $CONF_APP_REF->{formatTimeStamp};
    my $reads_indx_f    = $CONF_APP_REF->{f_reads_indx};
    my $resultsDir = File::Spec->catdir( $CONF_PROJ_REF->{dirBuildExport}, 'sets', $exportSet );
    my $isKeepAllReads = ((defined $CONF_PROJ_REF->{sortKeepAllReads}) and $CONF_PROJ_REF->{sortKeepAllReads});
    my $isKeepFilteredReads = $isKeepAllReads;
    my $isKeepUnmappedReads = $isKeepAllReads;

    my $regExpMitochondria = $CONF_PROJ_REF->{regExpMitochondria};
    my $readSampleRateInput = $CONF_PROJ_REF->{readSampleRateInput};
    my $minSEScore          = $CONF_PROJ_REF->{QVCutoffSingle};

    my $NMNM_DIR  = File::Spec->catdir( $resultsDir, $NMNM_TAG );
    my %Fields    = %exportFields;
    my $fieldSize = scalar( keys %Fields );

    my $binSize          = $CONF_PROJ_REF->{binSizeBuild};

     # TODO: Make it possible to get the filtered tile
    #my %isFilteredTile = getFilteredTileHash($buildFilePath, $exportFileDir, $lane);

    my $IDXFH = getIDXFH( File::Spec->catfile( $resultsDir, $reads_indx_f ) );

    my $confDir = $CONF_APP_REF->{dirConf};
    my $dirsPath = File::Spec->catfile( $projectDir , $confDir , "project.dirs.tar" );
    my $cmd = "tar -C $resultsDir/ -xf $dirsPath";
    executeCmd( $cmd, 5 );

    my $unmap;
    $unmap = getUnmapFH($CONF_APP_REF,$NMNM_DIR) if($isKeepUnmappedReads);

    # my $rmHandle = undef;
    # if ( 'RNA' eq $CONF_PROJ_REF->{applicationType} )
    # {
    #     $rmHandle = IO::File->new(">$resultsDir/RM/rm.txt")
    #       || errorExit "ERROR: Couldn't create file handle for $resultsDir/RM/rm.txt $!\n";
    # }

    my %exportStorage = ();
    $exportStorage{ref2Output}         = \%OUTPUTFILES;
    $exportStorage{resultsDir}         = $resultsDir;
    $exportStorage{binSize}            = $binSize;
    my $needRunIdFix = 0;
    
    my $needRunConfigData = 1;
    my $instrumentName;
    my $runID;
    my $laneNumber;
    my $barcodeSequence;
    
    #my $machineName  = $runsConfigRef->{run}{$runId}{machine};
    #my $exportRunId  = $runsConfigRef->{run}{$runId}{exportRunId};

    #if ( $runsConfigRef->{run}{$runId}{fixId} ne 'none' ) {
    #    $needRunIdFix = 1;
    #}

    # try to replicate splice junction counting logic for each lane
    # during spliced read splitting:
    #
    # this structure stores: {chrom}{splice_site_id}{count}
    #
    my ($read1SeqLength, $read2SeqLength) = ('','');
    my %spliceJunctionCount  = ();
    my $exportFile;
    foreach my $exportFile1 (grep {$_} ($read1ExportFilename, $read2ExportFilename))
    {
        $exportFile = $exportFile1;
        printLog( "Looking at $exportFile1\n", 4 );
        open( EXPT_1, "gunzip -fc $exportFile1 |" ) || errorExit "ERROR: Could not open export file $exportFile1 $!\n";
        my $readSeqLength='';
        while (<EXPT_1>) {
    
            $originalReads++;
            next if($readSampleRateInput && $originalReads * $readSampleRateInput <= $includedReads);
            ++$includedReads;
    
            chomp;
            my $line = $_;
            my @read1 = split /\t/, $line;    # read1 array
            if ( scalar(@read1) != $fieldSize ) {
                errorExit "ERROR: parsing $exportFile1 file at $originalReads (wrong number of fields "
                  . scalar(@read1) . " != " . $fieldSize . ")\n";
            }
    	
    	# retrieve the instrument name, run ID, lane, and barcode
    	if($needRunConfigData) {
    	    $instrumentName  = $read1[$exportMachineNameIndex];
    	    $runID           = $read1[$exportRunNumberIndex];
    	    $laneNumber      = $read1[$exportLaneIndex];
    	    $barcodeSequence = $read1[$exportBarcodeIndex];
    	    $needRunConfigData = 0;
    	}
    	
            $totalReads++;
            $readSeqLength = length($read1[$exportSequenceIndex]);
            $totalBases += $readSeqLength;
    
            # Ignore tiles specified as filtered out.
            # TODO: Make it possible to get the filtered tile
            #if (exists($isFilteredTile{$read1[$exportTileIndex]})) {
            #    ++$failedTileFilt;
            #    next;
            #}
    
            #if ( $needRunIdFix == 1 ) {
            #    ## those reads need to be fixed (machine name and run id)
            #    $read1[ $Fields{'MachineName'} ] = $machineName;
            #    $read1[ $Fields{'RunNumber'} ]   = $exportRunId;
            #}
    
            # make adjustment for SAM:
            if(($read1[ $exportPosnIndex ] ne '') and
               ($read1[ $exportPosnIndex ] < 1 ) and
               (not isSpliceJunctionChrom( $read1[$exportChrIndex], %$CONF_PROJ_REF)))
            {
                unMatchRead(@read1);
            }
    
            my $filterVal1  = $read1[$exportFilterIndex] eq 'N';
            my $chrVal1     = $read1[$exportChrIndex];
            my $strandVal1  = $read1[$exportStrandIndex];
            my $strandPosn1 = $read1[$exportPosnIndex];
    
            if ( $filterVal1 ) {    # ignore any reads which fail the filters
                $failedFilt++;
                next unless($isKeepFilteredReads);
            }    # if
    
            my $unmapTag;
    
            if ( $chrVal1 eq "QC" )
            {
                $failedQC++;
                $unmapTag="TAG_qcFail";
            }
            elsif ( $chrVal1 =~ /\d+:\d+:\d+/ )
            {    # ignore any read which doesn't aligns uniquely
                $nonUniqueAlign++;
                $unmapTag="TAG_nonUnique";
            }
            elsif ( $chrVal1 eq "NM" )  # grab the NM matches and put them to one side
            {                        # in case it are needed for de novo
                $nonMappedAlign++;
                $unmapTag="TAG_noMatch";
            }
            elsif ( $chrVal1 eq "RM" )
            { # grab the RM matches and put them to one side
                $riboCount++;
                $unmapTag="TAG_rm";
            }
            elsif ( $chrVal1 eq "CONTROL" )
            {
                $controlCount++;
                $unmapTag="TAG_control";
            }
    
            if(defined $unmapTag)
            {
                if($isKeepUnmappedReads) {
                    my $unmapFH = $unmap->{$unmapTag};
                    print $unmapFH join("\t", @read1),"\n";
                }
                next;
            }
    
            if (   $chrVal1 =~ /$regExpMitochondria/)
            {    # count mitochondria (there is no standard for Mitochondria name)
                $mitoCount++;
            }    # if
    
            my $readSras = $read1[ $exportSingleScoreIndex ];
            if ( !$readSras )
            {   # alignment score is zero:
                #$singleExclude++;
            }
    
            if ( !isSpliceJunctionChrom($chrVal1, %$CONF_PROJ_REF) )
            {
                outputSEread( %exportStorage, $strandPosn1, $chrVal1, @read1, "single" );
            }
            else
            {
                $splicedReads++;
                my $isSpliceFail = 1;
                if ( 'readBases' eq $CONF_PROJ_REF->{rnaCountMethod} )
                {
                    my $isCountable = ((not $filterVal1) and ($readSras >= $minSEScore));
                    $isSpliceFail = outputSpliceJunctionReadGap (%exportStorage, @read1, $totalReads,%spliceJunctionCount,$isCountable);
                }
                if($isSpliceFail)
                {
                    if($isKeepUnmappedReads) {
                        unMatchRead(@read1);
                        my $unmapFH=$unmap->{TAG_noMatch};
                        print $unmapFH join("\t", @read1),"\n";
                    }
                    next;
                }
            }
    
            if(not $filterVal1)
            {
                ++$usedReads; #count only the orignal read (not the generated parts)
                $usedBases += $readSeqLength;
                if ($readSras >= $CONF_PROJ_REF->{QVCutoffSingle})
                {
                    ++$goodReads;
                    $goodReadBases += $readSeqLength;
                }
            }
        }
        close(EXPT_1);
        $read1SeqLength = $readSeqLength if ($exportFile1 eq $read1ExportFilename);
        $read2SeqLength = $readSeqLength if ($exportFile1 eq $read2ExportFilename);
    }
    
    
    if($isKeepUnmappedReads) {
        for(values %$unmap) { close }
    }
    flushReadFiles(%OUTPUTFILES, 1);

    my $date = strftime $timeStampFormat, localtime;
    print $IDXFH "$readsIdxHeader\n";
    print $IDXFH "$totalReads"
      . "\t$totalBases"
      . "\t$usedReads"
      . "\t$usedBases"
      . "\t$goodReads"
      . "\t$goodReadBases"
      . "\t$failedTileFilt"
      . "\t$failedFilt"
      . "\t$failedQC"
      . "\t$singleExclude"
      . "\t$nonUniqueAlign"
      . "\t0" #$goodPairs
      . "\t0" #$mixedPairs
      . "\t$riboCount"
      . "\t$mitoCount"
      . "\t$splicedReads"
      . "\t" . ($CONF_PROJ_REF->{skipVariableMetadata} ? '' :  dirname($exportFile))
      . "\t$exportSet" #$lane"
      . "\t$date"
      . "\t$instrumentName"
      . "\t$runID"
      . "\t$laneNumber"
      . "\t$barcodeSequence"
      . "\t" #. (defined $medianFragmentLength ? $highSDFragmentLength : '')
      . "\t" #. (defined $medianFragmentLength ? $lowSDFragmentLength : '')
      . "\t" #. (defined $medianFragmentLength ? $maxFragmentLength : '')
      . "\t" #. (defined $medianFragmentLength ? $medianFragmentLength : '')
      . "\t" #. (defined $medianFragmentLength ? $minFragmentLength : '')
      . "\t" #. (defined $medianFragmentLength ? $nominalOrientation : '')
      . "\t$read1SeqLength"
      . "\t$read2SeqLength"
      . "\n" unless $needRunConfigData; #don't print the line if we did not parse a single read
      
    close($IDXFH);

    return if ( 'RNA' ne $CONF_PROJ_REF->{applicationType} );

    # if RNA mode, write out splice junction counts:
    #
    for my $chrom (keys %spliceJunctionCount) {
        my $chromDir = File::Spec->catdir($resultsDir, $chrom);
        errorExit "ERROR: can't find input lane chromosome directory: $chromDir\n" if(not -d $chromDir);
        my $countPath = File::Spec->catfile($chromDir, $tmpSpliceCountFileName);
        open(my $countFH,">$countPath") or
          errorExit "ERROR: can't write to file: $countPath\n";
        my $cRef = $spliceJunctionCount{$chrom};
        print $countFH $_ ."\t". $cRef->{$_} ."\n" for (keys %$cRef);
        close($countFH);
    }
}



=pod

=head1 The procedure commits all changes after export files has been split

=over 4

=item commitExportSplit

The procedure commits all changes after export files has been split. At this point
we can start marging bins in each lane.

Parameters:
    $projectDir       - CASAVA project directory 
    $CONF_PROJ_REF    - HASH MAP REF with application specific configuration
    $CONF_APP_REF     - HASH MAP REF with project specific configuration
    $runsConfig       - HASH MAP REF with runs specific configuration
	
Returns:
    genome sequence
        
=back

=cut

sub commitExportSplit($\%\%\%) {

    my ( $projectDir, $CONF_PROJ_REF, $CONF_APP_REF, $runsConfigRef ) = @_;
    my $exportSetsDir = File::Spec->catdir( $CONF_PROJ_REF->{dirBuildExport}, 'sets');
    
    # (1) merge all lane read index files into a single copy in 'stats'
    #
    my $reads_indx_f = $CONF_APP_REF->{f_reads_indx};
    my $dscIndex = File::Spec->catfile( $projectDir, 'stats', $reads_indx_f );
    executeCmd( "echo '$readsIdxHeader' > $dscIndex", 5 );
    
    for my $i (1 .. scalar(@{$runsConfigRef->{exportFiles}})) {
        my $setDir =  File::Spec->catdir( $exportSetsDir, $i );
        errorExit "ERROR: Could not find the set directory $setDir $!\n" unless ( -d $setDir );
        my $srcIndex = File::Spec->catfile( $setDir, $reads_indx_f );
        executeCmdEx( "cat $srcIndex | grep -vE '^#' >> $dscIndex", 0, 5 );
    }
    
    errorExit("ERROR: The provided samples had absolutely no data.") 
        if executeCmdEx( "cat $dscIndex | grep -vE '^#' 1>&2", 0, 5 );

    return if ( 'RNA' ne $CONF_PROJ_REF->{applicationType} );

    # (2) for RNA mode, merge all splice junction counts for each
    # chrom into a single copy:
    #
    my $currentBuildDir = $CONF_PROJ_REF->{dirBuildParsed};
    my %buildChromsBinSizes = ();
    readProjectParameters( %buildChromsBinSizes, "BUILD_BIN_SIZES", $projectDir);
    my @chroms = keys %buildChromsBinSizes;

    for my $chrom (@chroms) {
        next if (isSpliceJunctionChrom($chrom, %$CONF_PROJ_REF));
        my %mergedCounts = ();
        
        for my $set (1 .. scalar(@{$runsConfigRef->{exportFiles}})) {
            my $exportSetDir = File::Spec->catdir($exportSetsDir, $set);
            errorExit "ERROR: Could not find the export set directory ($exportSetDir) $!\n" unless ( -d $exportSetDir );
            my $inCountPath = File::Spec->catfile( $exportSetDir, $chrom, $tmpSpliceCountFileName );
            
            next unless(-e $inCountPath);
            
            open(my $CFH,"<$inCountPath") or
            errorExit "ERROR: can't open file: $inCountPath\n";
            
            while(<$CFH>) {
                chomp;
                my @s = split(/\t/);
                $mergedCounts{$s[0]} += $s[1];
            }
            
            close($CFH);
        }
        
        my $outCountPath = File::Spec->catfile($currentBuildDir,$chrom,$tmpSpliceCountFileName);
        open(my $OFH,">$outCountPath") or
          errorExit "ERROR: can't write to file: $outCountPath\n";
        print $OFH $_ . "\t" . $mergedCounts{$_} . "\n" for (keys %mergedCounts);
        close($OFH);
    }
}

1;    # says use was ok
__END__

