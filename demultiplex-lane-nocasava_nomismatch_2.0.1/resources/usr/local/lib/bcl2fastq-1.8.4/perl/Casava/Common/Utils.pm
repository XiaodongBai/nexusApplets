package Casava::Common::Utils;
use Carp;
use Cwd qw(realpath);
use File::Spec;
use Exporter;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Common::IOLib qw(executeCmd);


@ISA=("Exporter");
@EXPORT =
    qw(
       &adjustScoreForIncorrectBase
       &apply_mismatch
       &bestAlignFromList
       &conditionally_init_sigpipe_handler
       &conditionally_report_sigpipes
       &convertProbToQual
       &convertProbCorrectToPhred
       &convertQualToProbCorrect
       &convertQualToProbWrong
       &convertPhredToProbCorrect
       &convertPhredToProbWrong
       &convertToAlignmentDescriptor
       &compareVersions
       &determine_missing_tile_nums
       &expandLegacyUseBasesString
       &expandSymbolicQualityString
       &expandUseBasesString
       &expandUseBases
       &explodeChrom
       &explodePartnerOffset
       &get_expected_tile_nums
       &get_paths
       &get_read_use_bases_str
       &get_tiles_by_lane
       &getAlignmentScore
       &getBustardRunInfo
       &getChromosomeSizes
       &getConfigurationFilePath
       &getCycleMap
       &getEOL
       &getExactScoreForRead
       &getFlowCell
       &getFlowCellId
       getInstrumentName
       getRunNumber
       &getMatchDetails
       &getMeanBaseQuality
       &getPhasingFilename
       &getPhasingOptions
       &getQualFromChar 
       &getReadID
       &getReadLengths
       &getReadStartCycles
       &getRestOfGenomeCorrection
       &getRunFolderFromBustard
       &getRunFolderFromGerald
       &getRunFolderName
       &getRunFolderPath
       &getRunParameter
       &getShortRunFolderName
       &getSmtFilter
       &getTileList
       &getTotalGenomeSize
       &getWorstBases
       &gnuplotImage
       &lane_prefix
       &lane_num
       &pad_seq_with_masked_cycles
       &parseAnalysisDirName
       parseRunFolderName
       &parseUseBases
       &pickBestSingleReadAlignment
       &scoreAlignFromNhood
       &sigpipe_handler
       &updateCycleOffsets
       expandMd
       compressMd
       reallyRealPath
       );


#------------------------------------------------------------------------------

# PROJECT: GERALD
# MODULE:  Common.pm
# AUTHOR:  A. J. Cox
#
# Copyright (c) 2007 Solexa
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

# Functions common to the GERALD module

=pod

=head1 NAME

GERALD::Common.pm - Functions common to the GERALD module

=head2 SYNOPSIS

use GERALD::Common.pm qw();  

=head2 AUTHORSHIP

Copyright (c) 2007 Solexa, 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

=head1 DESCRIPTION

=cut

use warnings 'all';
use strict;
use POSIX;
use XML::Simple;
use constant DEBUG => 0; # set to 1 to get debug info
use Casava::Common::Log;

sub getShortRunFolderName(;$$);
sub getRunFolderPath($);
sub getRunFolderFromBustard($);
sub getRunFolderFromGerald($);
sub getTileList($$);
sub getFlowCell($$);
sub getRunFolderName($);
sub getFlowCellId($);
sub getInstrumentName($);
sub getRunNumber($);
sub gnuplotImage($$);
sub getRunParameter($$$);
sub getReadID($$$$$$$$;$);
sub expandUseBasesString($);
sub expandUseBases($$);
sub getChastityThreshold($$);
sub expandLegacyUseBasesString($$);
sub parseUseBases($;$$);
sub get_read_use_bases_str($$$$);
sub pad_seq_with_masked_cycles($\@$$);
sub getCycleMap($\@);
sub compareVersions($$);
sub updateCycleOffsets(\@$$);
sub parseAnalysisDirName($$);
sub getBustardRunInfo($$$\@\%);
sub getBustardRunInfo2($$\@\%);
sub getPhasingOptions($$\%);
sub getReadStartCycles($$\@);
sub getReadLengths($$\@);
sub getPhasingFilename($$;$);
sub convertProbToQual($);
sub convertProbCorrectToPhred($);
sub getQualFromChar($);
sub convertQualToProbWrong($);
sub convertPhredToProbWrong($);
sub convertProbCorrectToPhred($);
sub convertQualToProbCorrect($);
sub convertPhredToProbCorrect($);
sub expandSymbolicQualityString($);
sub getMatchDetails($$$);
sub getExactScoreForRead($);
sub apply_mismatch($$$$);
sub getWorstBases($$);
sub adjustScoreForIncorrectBase($);
sub getExactScoreForRead($);
sub pickBestSingleReadAlignment($$$$$$$);
sub bestAlignFromList($$$$$$$);
sub scoreAlignFromNhood($$$$$$$$);
sub getRestOfGenomeCorrection($$);
sub getMeanBaseQuality($$);
sub getAlignmentScore($$$);
sub getChromosomeSizes($$);
sub getTotalGenomeSize(\%);
sub getRestOfGenomeCorrection($$);

sub get_tiles_by_lane($$\@);
sub get_expected_tile_nums($$\@);
sub determine_missing_tile_nums($\%$\@);
sub get_paths($$$$\@);

sub getAutoFilename(\%$$$);
sub expandMd($$$);
sub compressMd($$);
sub explodeChrom($$);
sub explodePartnerOffset($$);

my $pathGnuplot = "gnuplot";
my $pathConvert = "convert";

=pod

=head1 FUNCTIONS

=head2 General Functions

=over 4

=cut


=pod

=item getShortRunFolderName($longRunFolderName, $readPairMode)

The procedure makes a short run folder name to use in the read ID
Strip off date and try to shorten other fields

B<Parameters:>

    $longRunFolderName   - long name of run folder
    $readPairMode        - [0,1 or 2] Optional readPairMode appends '-1' or '-2' to readID

B<Returns:> 

    short name of run folder

=cut
sub getShortRunFolderName(;$$)
{
    my ( $longRunFolderName, $readPairMode ) = @_;
    return "s" unless defined($longRunFolderName);
    $longRunFolderName=~s/^\d{6}_//; # strip date
    $longRunFolderName=~s/^slxa-//; # strip unnecessary machine name
    $longRunFolderName=~s/_0+/_/; # strip leading zeroes
# remove any _FCxxxx idiocy at end
# TC 18.6.7 - also remove any _R1 or _R2 suffixes. Issue is that both read 1
# and read 2 end up suffixed _R1 because everything gets named after the first
# run folder
    $longRunFolderName=~s/_[Ff][Cc]\d+(_[Rr][12])?$//; 

    if (defined($readPairMode))
    {
        die "Invalid value $readPairMode for readPairMode - must be 0,1 or 2"
            unless (($readPairMode eq '0')
                    ||($readPairMode eq '1')
                    ||($readPairMode eq '2'));
        $longRunFolderName.="-$readPairMode"
            if (($readPairMode eq '1')||($readPairMode eq '2'));
    } # if
    return $longRunFolderName;
} # sub getShortRunFolderName

=pod

=item getRunFolderPath($path)

The procedure extracts the name of the run-folder from a Base Calls directory.

B<Parameters:>

    $path                - base calls directory path

B<Returns:> 

    run folder path

=cut
sub getRunFolderPath($) {
    my ( $path ) = @_;

    my $paramsfile    = getConfigurationFilePath($path);
    my @path          = grep{ $_ } File::Spec->splitdir($path);
    my @emptyArray    = ();
    my %emptyHashMap  = ();
    return getBustardRunInfo($paramsfile, $path[-1], 'RunParameters/RunFolder', @emptyArray, %emptyHashMap);
} # sub getRunFolderPath

=pod

=item getRunFolderFromBustard($path)

The procedure extracts the name of the run-folder from a Base Calls directory.

B<Parameters:>

    $path                - base calls directory path

B<Returns:> 

    run folder name

=cut
sub getRunFolderFromBustard($) {
    my ( $path ) = @_;
    my $name = undef;

    my $runFolderPath = getRunFolderPath($path);
    return undef  unless defined $runFolderPath;
    my @dirs = grep{ $_ } File::Spec->splitdir($runFolderPath);
    my $runFolder = $dirs[-1];
    if ($runFolder =~ m/(\d{6}_[^_]+_\d+.*)/) {
        $name = $1;
    } # if
    if ($runFolder =~ m/(\d{8})/) {
        $name = $1;
    } # if
    return $name;
} # sub getRunFolderFromBustard

=pod

=item getRunFolderFromGerald($path)

The procedure extracts the name of the run-folder from an Alignment directory.

B<Parameters:>

    $path                - alignment directory path

B<Returns:> 

    run folder name

=cut
sub getRunFolderFromGerald($) {
    my ( $path )      = @_;
    my $geraldXmlPath = File::Spec->catfile( $path, 'config.xml' );

    my $geraldXmlRef = new XML::Simple(suppressempty => '', XMLDecl => 1);

    my $configHash = $geraldXmlRef->XMLin($geraldXmlPath);
    
    return undef unless exists $configHash->{ChipWideRunParameters}->{EXPT_DIR};
    return getRunFolderFromBustard($configHash->{ChipWideRunParameters}->{EXPT_DIR});
} # sub getRunFolderFromGerald

=pod

=item getTileList($folder, $suffix)

The procedure gets a list of filenames from a Gerald folder with a certain suffix
indexed by lanes and tile

B<Parameters:>

    $folder             - path to Gerald folder
    $suffix             - tile suffix

B<Returns:> 

    HASH MAP with list of tiles

=cut
sub getTileList($$) {
    my ( $folder, $suffix ) = @_;

    my $file;
    my %results;
    opendir (DIR, $folder) or die "can't open ${folder}.\n";
    while (defined($file = readdir (DIR))) {
        next unless ($file=~/s_([1-9])_([0-9]{4})${suffix}/);
        $results{$1}{$2} = $file;
    } # while

    return %results;
} # getTileList

=pod

=item gnuplotImage($gnuplotInputFile, $imageFile)

The procedure generates image from gnuplot script file.

B<Parameters:>

    $gnuplotInputFile   - path to gnuplot script
    $imageFile          - path to image

B<Returns:> 

    Nothing

=cut
sub gnuplotImage($$) {
    my ( $gnuplotInputFile, $imageFile ) = @_;

    executeCmd("${pathGnuplot} ${gnuplotInputFile} | ${pathConvert} ps:- $imageFile 2>&1");
}

=pod

=item getReadID($machineName, $runNumber, $flowCellId, $lane, $tile, $x, $y, $readPairMode, $indexBases)

The procedure generates read id string
Make a read ID of the form b2_537_FC123:2:8:222:333
For paired reads use b2_537_FC123:2:8:222:333-1 or b2_537_FC123:2:8:222:333-2
Input - machine name, run number, flow cell id, lane, tile, x, y, read pair mode, 
(optionally) index id
readPairMode appends '-1' or '-2' to readID

B<Parameters:>

    $machineName        - machine name
    $runNumber          - run number
    $flowCellId         - flow cell id
    $lane               - lane id
    $tile               - tile id
    $x                  - x coordinates
    $y                  - y coordinates
    $readPairMode       - [0,1 or 2] readPairMode appends '-1' or '-2' to readID
    $indexBases         - index id

B<Returns:> 

    string with readId

=cut
sub getReadID($$$$$$$$;$)
{
    my ($machineName, $runNumber, $flowCellId, $lane, $tile, $x, $y,
        $readPairMode, $indexBases)=@_;

    my $paddedRunNumber = sprintf("%04d", $runNumber);
    my $fullRunId = "${machineName}_${paddedRunNumber}";
    $fullRunId = "${fullRunId}_${flowCellId}"  if $flowCellId;
    my $readID = "${fullRunId}:${lane}:${tile}:${x}:${y}";
    if (defined($indexBases))
    {
        $indexBases=~s/\./N/g;
        $readID.='#'.$indexBases;
    } # if

    if (defined($readPairMode))
    {
        die "Invalid value $readPairMode for readPairMode - must be 0 or the read number"
            unless ($readPairMode =~ /^\d+$/);
        $readID.="/$readPairMode" if ($readPairMode ne '0');
    } # if
    return $readID;
} # sub getReadID

=pod

=item expandUseBasesString($useBasesString)

The procedure converts a USE_BASES parameter from short to long form
eg nY5N2y3 becomes nYYYYYNNyyy

B<Parameters:>

    $useBasesString     - short version of useBasesString 

B<Returns:> 

    string with used bases

=cut
sub expandUseBasesString($)
{
    my ( $useBasesString ) = @_;
    my $a;
    while($useBasesString=~/([YyNnIi])(\d+)/) 
    { 
        $a=$1 x $2; 
        $useBasesString=~s/$&/$a/; 
    } # while
    return $useBasesString;
} # sub expandUseBasesString

=pod

=item expandUseBases($mask, $readLength)

The procedure converts a USE_BASES parameter from short to long form
eg nY5N2y* becomes nYYYYYNNyyy on 11 bases long read

B<Parameters:>

    $mask       - short version of useBasesString with up to 1 '*'
    $readLength - non-zero read length to expand '*'

B<Returns:> 

    string with used bases

=cut
sub expandUseBases($$)
{
    my ($mask, $readLength) = @_;
    my $useBases = '';
    my $originalMask = $mask;
    my $filler;
    # Note: The order in (\*|\d*) DOES matter!
    while (0 < length($mask) and $mask =~ /^(y|Y|n|N|i|I)(\*|\d*)(.*)$/)
    {
        if ('*' eq $2)
        {
            errorExit("ERROR: more than one filler ('*') in USE_BASES mask: $originalMask") if $filler;
            $filler = {position => length($useBases), value => lc $1};
        }
        else
        {
            my $count = $2;
            $count = 1 unless $count;
            $useBases .= join('', map {lc $1} (1..$count));
        }
        $mask = $3;
    }
    errorExit("ERROR: Invalid USE_BASES maks: $originalMask") unless 0 == length($mask);
    if ($filler)
    {
        my $count = $readLength - length($useBases);
        my $s1 = substr $useBases, 0, $filler->{position};
        my $s2 = substr $useBases, $filler->{position};
        $useBases = $s1 . join('', map {$filler->{value}} (1..$count)) . $s2;
    }
    # pad with 'n' until the end of the read
    $useBases .= join('', map {'n'} (length($useBases)..($readLength - 1)));
    return $useBases;
}


=pod

=item getChastityThreshold($paramsfile, $bustarddir)

The procedure to extracts the value of the chastity threshold.

B<Parameters:>

    $paramsfile        - name of bustard config file 
    $bustarddir        - path to bustard directory 

B<Returns:> 

    chastity threshold value

=cut
sub getChastityThreshold($$)
{
    my ( $paramsfile, $bustarddir ) = @_;
    my @emptyArray  = ();
    my %emptyHashMap  = ();
    my $paramsRef    = getBustardRunInfo($paramsfile, $bustarddir, 'BaseCallParameters', 
        @emptyArray, %emptyHashMap);
    return $paramsRef->{ChastityThreshold};
} # sub getChastityThreshold

=pod

=item getFlowCell($paramsfile, $bustarddir)

The procedure to extracts the value of the flow cell type.

B<Parameters:>

    $paramsfile        - name of bustard config file
    $bustarddir        - path to bustard directory

B<Returns:>

    the string representation of the flow-cell type (should be "1.4mm" or "v4")

=cut
sub getFlowCell($$)
{
    my ( $paramsfile, $bustarddir ) = @_;
    my @emptyArray  = ();
    my %emptyHashMap  = ();
    my $paramsRef    = getBustardRunInfo($paramsfile, $bustarddir, 'RunParameters',
        @emptyArray, %emptyHashMap);
    return $paramsRef->{FlowCell} if exists $paramsRef->{FlowCell};
    # try the AnalysisInfo.xml file
    my ($volume,$directories,$file) = File::Spec->splitpath( $paramsfile );
    my @dirs = File::Spec->splitdir( $directories );
    my $analysisInfoPath = File::Spec->catfile(@dirs[0..$#dirs-4], "Config", "AnalysisInfo.xml");
    return undef unless -f $analysisInfoPath;
    my $analysisInfoRef = XMLin($analysisInfoPath);
    croak "$analysisInfoPath: missing QualityTable" unless exists $analysisInfoRef->{QualityTable};
    my $qualityTable = $analysisInfoRef->{QualityTable};
    $qualityTable = "qtable-v4.txt" if $qualityTable eq "qtable-1.6mm.txt";
    croak "$qualityTable: quality table should be either qtable-1.4mm.txt or qtable-v4.txt" unless $qualityTable =~ /^qtable-(1.4mm|v4).txt$/;
    return "$1";
} # sub getFlowCell

=pod

=item parseRunFolderName($name)

The procedure splits run folder name into the components.

B<Parameters:>

    $name        - RunFolder element contents from name of bustard config file

B<Returns:>

    undef or a referenoce to the following array: 
        (date 6 digits, instrument, run number, flow cell id)

=cut
sub parseRunFolderName($)
{
    my ( $name ) = @_;
    return undef unless $name and $name =~ /^(\d{6})_([^_]*)_(\d+)_([^_]*)/;
    my @ret = ($1, $2, $3, $4);
    return \@ret;
}

=pod

=item getRunFolderName($paramsfile)

The procedure to extracts the value of the flow cell id.

B<Parameters:>

    $paramsfile        - name of bustard config file

B<Returns:>

    undef or the string content of RunFolder element with (possible) directory path removed.

=cut
sub getRunFolderName($)
{
    my ( $paramsfile ) = @_;
    my @emptyArray  = ();
    my %emptyHashMap  = ();
    my $paramsRef    = getBustardRunInfo2($paramsfile, 'RunParameters',
        @emptyArray, %emptyHashMap);
    return '' unless exists $paramsRef->{RunFolder} and defined $paramsRef->{RunFolder};
    my @dirs = File::Spec->splitdir($paramsRef->{RunFolder});
    my $runFolder = $dirs[-1];
    return $runFolder;
}

=pod

=item getFlowCellId($paramsfile)

The procedure to extracts the value of the flow cell id.

B<Parameters:>

    $paramsfile        - name of bustard config file

B<Returns:>

    the string representation of the flow-cell id

=cut
sub getFlowCellId($)
{
    # First, try the RunParameters/RunFlowcellId element
    my @emptyArray = ();
    my %emptyHash = ();
    my $paramsRef = getBustardRunInfo2($_[0], 'RunParameters', @emptyArray, %emptyHash);
    my $flowCellId = undef;
    $flowCellId = $paramsRef->{RunFlowcellId} if exists $paramsRef->{RunFlowcellId};
    return $flowCellId if defined $flowCellId and 0 < length($flowCellId);
    # Fallback to parsing the run folder name
    my $runFolderName = getRunFolderName($_[0]);
    my $parsedName = parseRunFolderName($runFolderName);
    return '' unless $parsedName;
    return $parsedName->[3];
}

=pod

=item getInstrumentName($paramsfile)

The procedure to extracts the value of the flow cell id.

B<Parameters:>

    $paramsfile        - name of bustard config file

B<Returns:>

    the string containing the instrument name parsed from Xml RunFolder element

=cut
sub getInstrumentName($)
{
    my $runFolderName = getRunFolderName($_[0]);
    my $parsedName = parseRunFolderName($runFolderName);
    return '' unless $parsedName;
    return (@$parsedName)[1];
}

=pod

=item getRunNumber($paramsfile)

The procedure to extracts the value of the flow cell id.

B<Parameters:>

    $paramsfile        - name of bustard config file

B<Returns:>

    the string representation of the run number parsed from Xml RunFolder element

=cut
sub getRunNumber($)
{
    my $runFolderName = getRunFolderName($_[0]);
    my $parsedName = parseRunFolderName($runFolderName);
    return '' unless $parsedName;
    return (@$parsedName)[2];
}

=pod

=item getRunParameter($paramsfile, $bustarddir, $parameter)

The procedure to extracts the value of the specified run parameter.

B<Parameters:>

    $paramsfile        - name of bustard config file
    $bustarddir        - path to bustard directory
    $parameter         - the run parameter of interest

B<Returns:>

    the string representation of the run parameter, or undef if not found

=cut
sub getRunParameter($$$)
{
    my ( $paramsfile, $bustarddir, $parameter ) = @_;
    my @emptyArray  = ();
    my %emptyHashMap  = ();
    my $paramsRef    = getBustardRunInfo($paramsfile, $bustarddir, 'RunParameters',
        @emptyArray, %emptyHashMap);
    return $paramsRef->{$parameter} if exists $paramsRef->{$parameter};
    return undef;
}

=pod

=item expandLegacyUseBasesString($useBasesString)

The procedure to extracts the value of the chastity threshold.

B<Parameters:>

    $useBases          - [all, old]
    $readLength        - read length

B<Returns:> 

    use bases string

=cut
sub expandLegacyUseBasesString($$)
{
    my ( $useBases, $readLength ) = @_;

    if ($useBases eq "all")
    {
        die "$0: USE_BASES 'all' specified but invalid read length $readLength given!"
            if ($readLength<=0); 
        $useBases= "Y" x $readLength;
    } # if
    elsif ($useBases eq "odd")
    {
        die "$0: USE_BASES 'odd' specified but invalid read length $readLength given!"
            if ($readLength<=0); 
        $useBases= "YN" x ($readLength);
        $useBases=substr($useBases,0,length($useBases)-1);
    } # elsif
    else
    {
        die "Use bases string $useBases must be odd, all or YNYN..."
            unless ($useBases=~/^[YNyn]+$/);
    } # else

    return $useBases;
} # sub expandLegacyUseBasesString

=pod

=item parseUseBases($useBases, $readLength, $readPairMode)

The procedure Parses a long form USE_BASES parameter into an array of indices

Inputs:
1. USE_BASES string - either: 
'all' - use all bases (probably an anachronism)
'odd' - use odd numbered bases (definitely an anachronism)
Some combination of 'Y', 'y', 'N' and 'n'

2. read length - used in 'all' or 'odd' modes to generate a string
of the right length, also used to check there are enough 'Y/y's in 
USE_BASES - again probably an anachronism, set to 0 to ignore

3. (optional) Which read of a paired read - either
'1' - first read: treat all 'y's as 'n's
'2' - second read: treat all 'Y's as 'N's

B<Parameters:>

    $useBases          - [all, odd]
    $readLength        - read length
    $readPairMode      - [0,1 or 2] readPairMode appends '-1' or '-2' to readID

B<Returns:> 

    ARRAY with bases index - Read used cycle inds

=cut
sub parseUseBases($;$$)
{
    my ($useBases, $readLength, $readPairMode)=@_;
    
    $readLength   = 0  unless defined $readLength;
    $readPairMode = 0  unless defined $readPairMode;
    my @baseIndex;

    $useBases = expandLegacyUseBasesString($useBases, $readLength);

    if ($readPairMode eq "1")
    {
        $useBases=~tr/y/N/;
    } # if
    elsif ($readPairMode eq "2")
    {
        $useBases=~tr/Y/N/;
    } # elsif

    $useBases = uc ($useBases);

    my @ub=split(//,$useBases);
  #  my @baseIndex;
    for (my $i=0;$i<@ub;$i++)
    {
        push (@baseIndex, $i) if ($ub[$i] eq "Y"); 
    } # if

    # Only check read length with use bases if positive read length given
    if ($readLength>0)
    {
        die "Need $readLength bases of base information, only ", 
        scalar(@baseIndex), " specified in $useBases"
            unless (scalar(@baseIndex)>=$readLength);
        @baseIndex=splice(@baseIndex,0,$readLength);
    } # if
    return @baseIndex;
} # sub parseUseBases

=pod

=item get_read_use_bases_str($full_use_bases_str, $orig_read_lengths_str,
                             $read_num, \$read_use_bases_str)

Extracts the section of a USE_BASES string that corresponds to the specified
read number (1-offset), given the set of original read lengths.

B<Parameters:>

    $full_use_bases_str      - USE_BASES 
    $orig_read_lengths_str   - colon-separated original read lengths
    $read_num                - read number (1-offset)
    $read_use_bases_str_ref  - reference to $read_use_bases to be derived

B<Returns:>

    1 if OK, dies on error.

=cut

sub get_read_use_bases_str($$$$)
{
    my ($use_bases_str, $orig_read_lengths_str, $read_num,
        $read_use_bases_str_ref) = @_;

    my @orig_read_lengths = split(/:/, $orig_read_lengths_str);
    my $start_cycle_ind = 0;

    for (my $read_ind = 0; $read_ind < ($read_num - 1); ++$read_ind) {
        $start_cycle_ind += $orig_read_lengths[$read_ind];
    }

    my $read_length = $orig_read_lengths[$read_num - 1];

    if (($start_cycle_ind + $read_length) > length($use_bases_str)) {
        my $end_cycle_ind = $start_cycle_ind + $read_length - 1;
        die("Attempt to extract past end of USE_BASES (`$use_bases_str' "
             . "$start_cycle_ind:$end_cycle_ind)\n");
        # return 0;
    }

    $$read_use_bases_str_ref
        = substr($use_bases_str, $start_cycle_ind, $read_length);

    return 1;
}

=pod

=item pad_seq_with_masked_cycles($seq_str, $read_used_cycle_inds_ref, 
$read_orig_length, $pad_ch)

Pads a supplied sequence string with the specified character at all 
positions previously masked out by a USE_BASES mask.
The USE_BASES mask must be specified in the array-of-cycle-indices format
returned by parseUseBases in conjunction with the original read length
(the same as the length of the USE_BASES string).

[The expectation is that USE_BASES will be pre-parsed before calling this
#  function many times.]

B<Parameters:>

    $seq_str                    - [all, old]
    $read_used_cycle_inds_ref   - ARRAY REF  Read used cycle inds
    $read_orig_length           - read length
    $pad_ch                     - padd with $pad_ch character

B<Returns:> 

    padded string

=cut
sub pad_seq_with_masked_cycles($\@$$)
{
    my ($seq_str, $read_used_cycle_inds_ref, $read_orig_length, $pad_ch) = @_;

    if (length($seq_str) != scalar(@$read_used_cycle_inds_ref)) {
        my $seq_str_len = length($seq_str);
        my $expected_length = scalar(@$read_used_cycle_inds_ref);

        die("Sequence `$seq_str' is of length $seq_str_len "
            . "- expected length $expected_length");
    } # if

    my $padded_str = $pad_ch x $read_orig_length;
    my @padded_str_chs = split(//, $padded_str);
    my @seq_chs = split(//, $seq_str);
    @padded_str_chs[@$read_used_cycle_inds_ref] = @seq_chs;
    $padded_str = join('', @padded_str_chs);

    return $padded_str;
} # sub pad_seq_with_masked_cycles

=pod

=item getCycleMap($use_bases_str, $cycle_map_ref)

The procedure assumes fully expanded USE_BASES, i.e. only [NYny]* with `Y's strictly 
before `y's.
Note that, for convenience, @cycle_map is indexed by (1-offset) cycle number 
so the first element (cycle 0) is undefined and should be ignored.

B<Parameters:>

    $use_bases_str              - used bases string see L<expandUseBasesString($useBasesString)>
    $cycle_map_ref              - the results will be stored in $cycle_map_ref

B<Returns:> 

    max read number 

=cut
sub getCycleMap($\@)
{
    my ($use_bases_str, $cycle_map_ref) = @_;

    my $read_num = 0;
    my $max_read_num = 0;

    my @cycle_nums;

    my @use_bases_chars = split(//, $use_bases_str);
    my $base_ind = 1;

    foreach my $use_bases_char (@use_bases_chars) {
        $read_num = 0;

        if ($use_bases_char eq 'Y') {
            $read_num = 1;

            if ($max_read_num == 2) {
                warn "Found `Y' after `y' in use_bases $use_bases_str"
                    . " after base $base_ind\n";
                return 0;
            }
        } elsif ($use_bases_char eq 'y') {
            $read_num = 2;

            if ($max_read_num == 0) {
                warn "Found `y' before any `Y' in use_bases $use_bases_str"
                    . " after base $base_ind\n";
                return 0;
            }
        } elsif (!(($use_bases_char eq 'N') || ($use_bases_char eq 'n'))) {
            warn "Found unexpected char `$use_bases_char' "
                . "in use_bases $use_bases_str after base $base_ind\n";
            return 0;
        }

        if ($read_num > $max_read_num) {
            $max_read_num = $read_num;
        } # if

        if (!defined($cycle_nums[$read_num])) {
            $cycle_nums[$read_num] = 0;
        } # if

        $cycle_map_ref->[$base_ind]{read_num} = $read_num;
        $cycle_map_ref->[$base_ind]{cycle_num} = ++$cycle_nums[$read_num];

        ++$base_ind;
    }

    return $max_read_num;
} # sub getCycleMap


=pod

=item updateCycleOffsets($use_bases_str, $cycle_map_ref)

The procedure updates Cycle Offsets
@cycle_offsets is indexed by 1-offset $read_num

B<Parameters:>

    $cycle_offsets_ref          - cycle Map - see L<getCycleMap>
    $curr_read_num              - ???
    $num_cycles_in_curr_read    - ???

B<Returns:> 

    Nothing 

=cut
sub updateCycleOffsets(\@$$)
{
    my ($cycle_offsets_ref, $curr_read_num, $num_cycles_in_curr_read) = @_;
        
    if (defined($cycle_offsets_ref->[$curr_read_num + 1])) {
        return;
    } # if

    if ($curr_read_num == 1) {
        $cycle_offsets_ref->[$curr_read_num] = 0;
    } # if

    if (!defined($cycle_offsets_ref->[$curr_read_num])) {
        die("Cannot calculate cycle offset for read $curr_read_num "
            . "as offset for previous read is unknown\n");
    } # if

    $cycle_offsets_ref->[$curr_read_num + 1]
        = $cycle_offsets_ref->[$curr_read_num] + $num_cycles_in_curr_read;
} # sub updateCycleOffsets

=pod

=item compareVersions($ver1_ref, $ver2_ref)

The procedure takes refs to two version triples (major, minor, release).
Returns : (v2 > v1) : 1; (v2 < v1) : -1; (v2 == v1) : 0
           or undef if either specified array is not a triple.

B<Parameters:>

    $ver1_ref          - string with version 1
    $ver2_ref          - string with version 2

B<Returns:> 

    -1, 0 or 1 

=cut
sub compareVersions($$)
{
    my ($ver1_ref, $ver2_ref) = @_;
    
    my @ver1 = @$ver1_ref;
    my @ver2 = @$ver2_ref;

    my $num_ver_pieces = 3;

    if ((scalar(@ver1) != $num_ver_pieces) 
        || (scalar(@ver1) != $num_ver_pieces)) {
        return undef;
    } # if

    for (my $piece_ind = 0; $piece_ind < $num_ver_pieces; ++$piece_ind) {
        my $ver1_piece = $ver1[$piece_ind];
        my $ver2_piece = $ver2[$piece_ind];

        if ($ver2_piece < $ver1_piece) {
            return -1;
        } elsif ($ver2_piece > $ver1_piece) {
            return 1;
        }
    } #  for

    return 0;
} # sub compareVersions

=pod

=item parseAnalysisDirName($dir_name, $dir_info_ref)

The procedure parses the name of an analysis directory Firecrest / Bustard / Gerald,
returning information extracted as a hash.
Possible keys : software_name, major_ver, minor_ver, release, user, date,
                 first_cycle, last_cycle

B<Parameters:>

    $dir_name          - path to Gerald folder
    $dir_info_ref      - HASH MAP REF where parsed info will be stored

B<Returns:> 

    status 0 or 1 

=cut
sub parseAnalysisDirName($$)
{
    my ($dir_name, $dir_info_ref) = @_;

    my @dir_pieces = split('_', $dir_name);

    if (@dir_pieces < 3) {
        carp("Found less than 3 parts in the directory name $dir_name");
        return 0;
    }

    my $date = $dir_pieces[-2];
    my @date_pieces = split('-', $date);

    if (@date_pieces != 3) {
        carp("Invalid date: $date");
        return 0;
    }

    $dir_info_ref->{user} = $dir_pieces[-1];
    $dir_info_ref->{date} = $date;

    my $soft_name_ver = $dir_pieces[-3];

    if ($soft_name_ver =~ /([a-zA-Z]+)(\d+)\.(\d+)\.(\d+)/) {
        $dir_info_ref->{software_name} = $1;
        $dir_info_ref->{major_ver} = $2;
        $dir_info_ref->{minor_ver} = $3;
        $dir_info_ref->{release} = $4;
    } elsif ($soft_name_ver =~ /([a-zA-Z]+)(\d+)\.(\d+)([a-z0-9]+)/) {
        $dir_info_ref->{software_name} = $1;
        $dir_info_ref->{major_ver} = $2;
        $dir_info_ref->{minor_ver} = $3;
        $dir_info_ref->{release} = 0;
        $dir_info_ref->{variant} = $4;
    } else {
        $dir_info_ref->{software_name} = $soft_name_ver;
    }
    
    if (@dir_pieces > 3) {
        my $cycle_str = $dir_pieces[-4];

        if ($cycle_str =~ /C(\d+)-(\d+)/) {
            $dir_info_ref->{first_cycle} = $1;
            $dir_info_ref->{last_cycle} = $2;
        }
    }

    return 1;
} # sub parseAnalysisDirName

=pod

=item getConfigurationFilePath($dir)

The procedure gets the path to the configuration file for Firecrest or Bustard.
Starting from PL version 1.3, the configuration file is the config.xml file in
the given directory. In previous versions, it was the .params file in the
parent directory.

Since the method can be called on non-existing directories, some logic has been 
added to provide an appropriate error processing, depending on the conditions:

- if the directory does not exist, produce a warning and return the path for the
current version of the pipeline

- if the directory exists, produce an error if none of the configuration files
are found.

B<Parameters:>

    $dir            - Full path of directory

B<Returns:> 

    Full path to the configuration file.

=cut

sub getConfigurationFilePath($)
{
    my ($dir) = @_;
    croak("ERROR: Undefined directory") unless defined $dir;
    my $configuration_file = File::Spec->catfile($dir, "config.xml");
    if (! -d $dir)
    {
        carp("Warning: the directory $dir does not exist") unless -d $dir;
        return $configuration_file;
    }
    return $configuration_file if -f $configuration_file;
    my $legacy_configuration_file = File::Spec->catfile($dir, "..", ".params");
    croak("ERROR: Couldn't find the configuration file for $dir") unless -f $legacy_configuration_file;
    return $legacy_configuration_file;
}

=pod

=item getBustardRunInfo($xml_file_path, $bustard_dir, $rel_xpath_str, 
$force_array_keys_ref, $key_attr_ref)

The procedure parses a Bustard configuration file (i.e. config.xml in PL>=1.3
or the .params in the parent parent Firecrest folder in PL<1.3)
returning a reference to the Run substructure for the specified
Bustard directory.

B<Parameters:>

    $xml_file_path          - configuration file path
    $bustard_dir            - Name of bustard directory
    $rel_xpath_str          - xpath expresion to get items from config file
    $force_array_keys_ref   - ARRAY REF force_array to XML::Simple
    $key_attr_ref           - HASH MAP REF $key_attr to XML::Simple

B<Returns:> 

    HASH MAP REF with parsed configuration file 

=cut
sub getBustardRunInfo2($$\@\%) {
    my ($xml_file_path, $rel_xpath_str, $force_array_keys_ref,
        $key_attr_ref) = @_;

    return undef  unless (-f $xml_file_path);

    my $xml_ref = new XML::Simple(suppressempty => '',
                                  XMLDecl => 1);

    my @force_array_keys = (@$force_array_keys_ref);
    my %key_attr_map = %$key_attr_ref;

    my $cfg_hash_ref = $xml_ref->XMLin($xml_file_path,
                                       forcearray => \@force_array_keys,
                                       keyattr => \%key_attr_map);
    return undef  unless defined $cfg_hash_ref;

    my @rel_xpath_dirs = split('/', $rel_xpath_str);
    my @layer_names = (@rel_xpath_dirs);
    carp "$xml_file_path does not have a 'Run' element" unless exists $cfg_hash_ref->{'Run'};
    my $curr_hash_ref = $cfg_hash_ref->{'Run'};
    my $next_hash_ref;

    foreach my $layer_name (@layer_names) {
        print "getBustardRunInfo: going though layer '$layer_name'...\n"  if DEBUG;
        my $next_hash_ref = $curr_hash_ref->{$layer_name};
        return undef  unless defined $next_hash_ref;

        $curr_hash_ref = $next_hash_ref;
    }

    return $curr_hash_ref;
}

sub getBustardRunInfo($$$\@\%) {
    my ($xml_file_path, $bustard_dir, $rel_xpath_str, $force_array_keys_ref,
        $key_attr_ref) = @_;

    return getBustardRunInfo2($xml_file_path, $rel_xpath_str, @$force_array_keys_ref,
        %$key_attr_ref);
}

=pod

=item getPhasingOptions($xml_file_path, $bustard_dir, $options_ref)

The procedure a Bustard configuration file (i.e. config.xml for PL>=1.3
or .params in the parent Firecrest folder in PL<1.3)
to extract a hash of phasing options indexed by read and then by option.

B<Parameters:>

    $xml_file_path          - Bustard configuration file path
    $bustard_dir            - Name of bustard directory
    $options_ref            - HASH MAP REF with phasing options indexed by read and then by option

B<Returns:> 

    status 0 or 1 

=cut
sub getPhasingOptions($$\%)
{
    my ($xml_file_path, $bustard_dir, $options_ref) = @_;

    my $rel_xpath_str = 'BaseCallParameters/Phasing';
    my @force_array_keys = ('Phasing');
    my %key_attr_map = ('Phasing' => 'Read');

    my $curr_hash_ref = getBustardRunInfo($xml_file_path, $bustard_dir,
                                          $rel_xpath_str,
                                          @force_array_keys, %key_attr_map);

    if (!defined($curr_hash_ref)) {
        return 0;
    }

    %$options_ref = %$curr_hash_ref;

    return 1;
} # sub getPhasingOptions

=pod

=item getReadStartCycles($xml_file_path, $bustard_dir, $read_start_cycles_ref)

The procedure a Bustard configuration file (i.e. config.xml in PL>=1.3 or
.params in the parent Firecrest folder in PL<1.3)
to extract an array of start cycles (e.g. two for Paired End reads)
for the specified Bustard directory.

B<Parameters:>

    $xml_file_path          - bustard configuration file path
    $bustard_dir            - Name of bustard directory
    $read_start_cycles_ref  - ARRAY MAP REF array of start cycles

B<Returns:> 

    status 0 or 1 

=cut
sub getReadStartCycles($$\@)
{
    my ($xml_file_path, $bustard_dir, $read_start_cycles_ref) = @_;
    
    my $rel_xpath_str = 'RunParameters/Reads';
    my @force_array_keys = ('Reads');
    my %key_attr_map = ('Reads' => 'Index');

    my $curr_hash_ref = getBustardRunInfo($xml_file_path, $bustard_dir,
                                          $rel_xpath_str,
                                          @force_array_keys, %key_attr_map);

    if (!defined($curr_hash_ref)) {
        return 0;
    }

    foreach my $read_key (sort keys %$curr_hash_ref) {
        my $read_hash_ref = $curr_hash_ref->{$read_key};
        my $first_cycle = $read_hash_ref->{'FirstCycle'};
        my $last_cycle = $read_hash_ref->{'LastCycle'};

        if (!(defined($first_cycle))) {
            return 0;
        }

        push(@$read_start_cycles_ref, $first_cycle);
    }

    return 1;
} # sub getReadStartCycles

=pod

=item getReadLengths($xml_file_path, $bustard_dir, $read_lengths_ref)

The procedure a Bustard configuration file
to extract an array of read lengths (e.g. two for Paired End reads)
 for the specified Bustard directory.

B<Parameters:>

    $xml_file_path          - Bustard configuration file path
    $bustard_dir            - Name of bustard directory
    $read_lengths_ref       - ARRAY MAP REF with array of read lengths

B<Returns:> 

    status 0 or 1 

=cut
sub getReadLengths($$\@)
{
    my ($xml_file_path, $bustard_dir, $read_lengths_ref) = @_;

    my $rel_xpath_str = 'RunParameters/Reads';
    my @force_array_keys = ('Reads');
    my %key_attr_map = ('Reads' => 'Index');

    my $curr_hash_ref = getBustardRunInfo($xml_file_path, $bustard_dir,
                                          $rel_xpath_str,
                                          @force_array_keys, %key_attr_map);

    if (!defined($curr_hash_ref)) {
        return 0;
    }

    foreach my $read_key (sort keys %$curr_hash_ref) {
        my $read_hash_ref = $curr_hash_ref->{$read_key};
        my $first_cycle = $read_hash_ref->{'FirstCycle'};
        my $last_cycle = $read_hash_ref->{'LastCycle'};

        if (!(defined($first_cycle) && defined($last_cycle))) {
            return 0;
        }

        my $read_length = ($last_cycle - $first_cycle) + 1;

        push(@$read_lengths_ref, $read_length);
    }

    return 1;
} # sub getReadLengths

=pod

=item getAutoFilename($option_ref, $lane_num, $prefix, $suffix)

Translation of def getAutoFilename (makefile.py). It is not exported

B<Parameters:>

    $option_ref             - 
    $lane_num               - 
    $prefix                 - 
    $suffix                 - 

B<Returns:> 

    filename 

=cut
sub getAutoFilename(\%$$$)
{
    my ($option_ref, $lane_num, $prefix, $suffix) = @_;
    
    my $filename;

    if ($option_ref->{AutoFlag} == 1) {
        if ($option_ref->{AutoLane} == 0) {
            $filename = sprintf("%ss_%02i%s",
                               $prefix, $option_ref->{Cycle}, $suffix);
        } else {
            $filename = sprintf("%ss_%i_%02i%s",
                               $prefix, $option_ref->{AutoLane},
                               $option_ref->{Cycle}, $suffix);
        }
    } elsif ($option_ref->{AutoFlag} == 2) {
        $filename = sprintf("%ss_%i_%02i%s",
                           $prefix, $lane_num, $option_ref->{Cycle}, $suffix);
    } else {
        $filename = $option_ref->{Path};
    }

    return $filename;
} # sub getAutoFilename

=pod

=item getPhasingFilename($phasing_dir, $option_ref, $lane_num)

Translation of def getPhasingFilename (makefile.py).

B<Parameters:>

    $phasing_dir            - name of Phasing directory e.g. Phasing
    $option_ref             - HASH MAP REF to structure from getPhasingOptions->{lane_num}
    $lane_num               - lane id 

B<Returns:> 

    filename 

=cut
sub getPhasingFilename($$;$)
{
    my ($phasing_dir, $option_ref, $lane_num) = @_;

    my $filename = getAutoFilename(%{$option_ref}, $lane_num,
                                   "${phasing_dir}/", "_phasing.xml");

    if (($option_ref->{AutoFlag} == 0) && (!$option_ref->{Path})) {
        $filename = "params${lane_num}.xml";
    }

    return $filename;
} # sub getPhasingFilename

=pod

=head2 Quality values and Probabilities Functions

Utilities for converting between quality values and probabilities

=cut

my $to10Base10=10.0/log(10.0); 

=pod

=item convertProbToQual($p)

The procedure converts prob(base call=X), X=A,C,G,T into a Solexa quality value
0.01 ~= Q20, 0.001 =~ Q30 etc.
Same formula works for all 4 base possibilities
Moreover our QV and the Phred QV round to the same integer for
base calls of Q10 and above

B<Parameters:>

    $p                      - probability

B<Returns:> 

    quality 

=cut
sub convertProbToQual($)
{
    my ($p)=@_;
    my $toBase10=10.0/log(10.0); # NB converts log(x) -> 10*log10(x)
    return $to10Base10*log($p/(1-$p));
} # convertProbToQual

# Convert probability that a called bases is correct into a Phred quality value
# 0.01 = Q20, 0.001 = Q30 etc.
# This formula only makes sense for a called base

=pod

=item convertProbCorrectToPhred($p)

The procedure converts probability that a called bases is correct into a Phred quality value
0.01 = Q20, 0.001 = Q30 etc.
This formula only makes sense for a called base

B<Parameters:>

    $p                      - probability

B<Returns:> 

    quality 

=cut
sub convertProbCorrectToPhred($)
{
    my ($p)=@_;
    my $toBase10=-10.0/log(10.0); # NB converts log(x) -> 10*log10(x)
    return -$to10Base10*log(1-$p);
} # convertProbCorrectToPhred

=pod

=item convertQualToProbCorrect($qv)

The procedure converts a Solexa quality value for X into the prob that base=X

B<Parameters:>

    $qv                     - quality

B<Returns:> 

    probability 

=cut
sub convertQualToProbCorrect($)
{
    my ($qv)=@_;
    return 1.0/(1+10**(-$qv/10.0));
} # convertQualToProbCorrect

=pod

=item convertQualToProbWrong($qv)

The procedure converts a Solexa quality value for X into the prob that base != X

B<Parameters:>

    $qv                     - quality

B<Returns:> 

    probability 

=cut
sub convertQualToProbWrong($)
{
    my ($qv)=@_;
    return convertQualToProbCorrect(-$qv);
} # convertQualToProbWrong

=pod

=item convertPhredToProbCorrect($qv)

The procedure converts a Phred quality value for X into the prob that base=X

B<Parameters:>

    $qv                     - quality

B<Returns:> 

    probability 

=cut
sub convertPhredToProbCorrect($)
{
    my ($qv)=@_;
    return (1.0-(10**(-$qv/10.0)));
} # convertPhredToProbCorrect

=pod

=item convertPhredToProbWrong($qv)

The procedure converts a Phred quality value for X into the prob that base != X

B<Parameters:>

    $qv                     - quality

B<Returns:> 

    probability 

=cut
sub convertPhredToProbWrong($)
{
    my ($qv)=@_;
    return 10**(-$qv/10.0);
} # convertPhredToProbWrong

=pod

=item getQualFromChar($qv)

The procedure converts quality to charater (- 64)

B<Parameters:>

    $qv                     - quality

B<Returns:> 

    character 

=cut
sub getQualFromChar($)
{
    my ($c)=@_;
    my $qual = ord($c)-64;
    return (($qual < 2) ? 2 : $qual);
} # getQualFromChar

=pod

=item expandSymbolicQualityString($qv)

The procedure expands a symbolic format quality string into a string of ASCII numbers

B<Parameters:>

    $qv                     - quality (list of characters)

B<Returns:> 

    string with list of numbers

=cut
sub expandSymbolicQualityString($)
{
    my $qv;
    for (my $i=0;$i<length($_[0]);$i++)
    {
        $qv.=sprintf("%d ",getQualFromChar(substr($_[0],$i,1)));
    } # for
    return $qv;
} # convertSymbolicQualityString($)

=pod

=item getMatchDetails($descriptor, $chrom, $nameSource)

The procedure parses match descriptor into its constituent fields, eg
BAC_plus_vector.fa:45265R18T13C3 goes to
BAC_plus_vector.fa, 45265, R, 18T13C3

B<Parameters:>

    $descriptor            - descriptor e.g.: BAC_plus_vector.fa:45265R18T13C3
    $chrom                 - chromosome name in case the descriptor does not contain it
    $nameSource            - fileName, contigName or undef.

B<Returns:> 

    string with list of numbers

=cut
my %explodeChromIndexMap = (fileName=>0, contigName=>1);

sub getMatchDetails($$$)
{
    my ($descriptor, $chrom, $nameSource)=@_;
    if ($descriptor=~s/^(.+)://)
    {
        $1 =~ /^(.*?)(\/(.*))?$/;
        my @chromContig = ($1,$3);
        if (defined $nameSource)
        {
            # This duplicates explodeChrom because otherwise the circularization
            # does not work for multi-entry fasta alignments.
            $chrom = $chromContig[$explodeChromIndexMap{$nameSource}];
            $chrom = $chromContig[0] unless defined($chrom);
            # the tab must not be appended at this point because PickBests don't expect it
        }
        else
        {
            # This case is exclusively for PickBestAlignmentRNA splice junction alignments.
            # It must put tab in as the string gets directly printed out
            $chrom = "$chromContig[0]\t";
            $chrom .= $chromContig[1] if defined $chromContig[1];
        }
    } # if
    
    if ($descriptor=~/^(-?\d+)([RF])([\^\$0123456789ACGTN]+)$/)
    {
#        $chrom1=$1; $pos1=$2; $strand1=$3; $type1=$4;
        return ($chrom, $1, $2, $3);
    } # if
    else
    {
        die "$0: Unable to parse $descriptor\n";
    } # else
} # sub getMatchDetails

=pod

=item apply_mismatch($new_seq, $num_cycles, $mismatch_str, $line_num)

The procedure applies a mismatch to a sequence.

B<Parameters:>

    $new_seq               - sequence where mismatch will be applied
    $num_cycles            - seqence length
    $mismatch_str          - mismatch string (6 field from SAF) 
    $line_num              - line id

B<Returns:> 

    corrected sequence 

=cut
sub apply_mismatch($$$$)
{
    my ($new_seq, $num_cycles, $mismatch_str, $line_num) = @_;

    my $ref_bases_str = 'ACGTN';

    my @mismatch_fields = split(/([$ref_bases_str])/, $mismatch_str);
    my $pos = 0;

    foreach my $mismatch_field (@mismatch_fields) {
        next if ($mismatch_field eq '');

        if ($pos >= $num_cycles) {
            die("Mismatch string `$mismatch_str' specifies position beyond "
                . "end of sequence (line $line_num).\n");
        }

        if ($mismatch_field =~ /[$ref_bases_str]/) {
            substr($new_seq, $pos++, 1) = $mismatch_field;
        } elsif ($mismatch_field =~ /\d+/) {
            $pos += $mismatch_field;
        } else {
            die("Unexpected field `$mismatch_field' in mismatch "
                . "string `$mismatch_str' (line $line_num)\n");
        }
    }

    return $new_seq;
}

=pod

=item convertToAlignmentDescriptor

The procedure computes a compressed match descriptor from two sequences.
This is a port from aligner.cpp convertToAlignmentDescriptor. IT DOES
NOT PRODUCE GAPS

B<Parameters:>

    $alignedA            - read sequence
    $alignedB            - reference sequence
    ;$refMismatchesCount  - reference which will receive number of mismatches

B<Returns:> 

    compressed match descriptor without gaps 

=cut
sub convertToAlignmentDescriptor($$;\$)
{
    my ($alignedA, $alignedB, $refMismatchesCount) = @_;
    my ($alignedALength, $alignedBLength) = (length($alignedA), length($alignedB));
    
    errorExit "Lengths are different when computing match descriptor read: $alignedA ref: $alignedB" 
        unless ( $alignedALength == $alignedBLength ); 

    $$refMismatchesCount = 0 unless !defined $refMismatchesCount;
    
    my $matches = 0;
    my $ret = '';
    my @arrAlignedA = split(//, $alignedA);
    my @arrAlignedB = split(//, $alignedB);

    for( my $i=0; $i < $alignedALength; ++$i ) 
    {
        if( $arrAlignedA[$i] eq $arrAlignedB[$i] ) 
        {
               ++$matches;
        }
        elsif( '-' eq $arrAlignedA[$i] || '-' eq $arrAlignedB[$i] ) 
        {
               errorExit "Gapped notation found in sequences when calculating non-gapped match descriptor";
        } 
        else
        {
           $ret .= $matches if ($matches);      
               $matches = 0;
               $ret .= $arrAlignedB[$i];
           ++$$refMismatchesCount unless !defined $refMismatchesCount;
        }
    }

    $ret .= $matches if ($matches);      
    return $ret;
}

=pod

=head2 expandMd

Build a pair of read and corresponding reference sequences from compressed 
match descriptor and indel-containing read sequence

I<Parameters:>

=over 4

=item *

$type

Match descriptor

=item *

$read

read sequence

=item *

$qual

quality values

I<Returns:>

Array containing fully expanded reference followed by read followed by the corresponding quality string

I<Exceptions:>

=over 4

=item *

Icorrect format format of the match descriptor

=item *

The resulting sequences happened to be of different lengths. This must not happen if 
everything is correct

=back

=cut


sub expandMd($$$)
{
    my ($type,$read,$qual) = @_;

    my $posInRead = 0;
    my $within_indel = 0;

    $type .= "X"; # otherwise we have special cases for the last element of the MD
    my $parsedPos = 0;

    my ($expandedRead,$expandedRef );

    while ($type =~/(\d*)([XNACGT\^\$])/g)
    {
        my( $save_1,$save_2 ) = ($1,$2);

        if( defined($save_1) )
        {
            $parsedPos += length( $save_1 );
        }
        if( defined($save_2) )
        {
            $parsedPos++;
        }
        

        if( $within_indel == 1 )
        {
            if( $save_2 eq "\$" )
            {
                if( defined($save_1) ) 
                { 
                    $expandedRead .= substr($read,$posInRead,$save_1) unless($save_1 eq ""); 
                    $expandedRef  .= "-" x $save_1 unless($save_1 eq ""); 

                    $posInRead+=$save_1 unless( $save_1 eq ""); 
                    
                }
                $within_indel = 0;
            }
            else
            {
                if( defined($save_1) )
                {
                    if( $save_1 ne "" )
                    {
                        errorExit "ERROR: mix of digits and characters within an indel!";
                    }
                } 

                $expandedRead .= "-";
                $expandedRef  .= $save_2;


            }
        }
        else
        {
            # we're not within an indel
            if( $save_2 eq "\^" )
            {
                if( defined($save_1) ) 
                {
                    $expandedRead .= substr($read,$posInRead,$save_1) unless($save_1 eq ""); 
                    $expandedRef  .= substr($read,$posInRead,$save_1) unless($save_1 eq ""); 
 
                    $posInRead+=$save_1 unless( $save_1 eq ""); 
                }

                $within_indel = 1;
            }
            else
            {
                if( defined($save_1) ) 
                { 
                    $expandedRead .= substr($read,$posInRead,$save_1) unless($save_1 eq ""); 
                    $expandedRef  .= substr($read,$posInRead,$save_1) unless($save_1 eq ""); 

                    $posInRead+=$save_1 unless($save_1 eq ""); 
                }

                if( $save_2 ne 'X' )
                {
                    $expandedRead .= substr($read,$posInRead,1);
                    $expandedRef  .= $save_2;
                }
                
                $posInRead++ if ( $save_2 =~ /[ACGTN]/ );



            }
        }



    } # while

#    print "expandedRead = $expandedRead\n" if $DEBUG == 1 || ($type =~ /N/); 
#    print "expandedRef  = $expandedRef\n" if $DEBUG == 1 || ($type =~ /N/);
    if( length($expandedRead) != length($expandedRef) )
    {
        errorExit "ERROR Expanded read and reference are of different length";
    }

    # now expand the quality
    my $expandedQuality = "";
    my $posInQuality = 0;
    for( my $i=0;$i<length($expandedRead);$i++ )
    {
        if( substr($expandedRead,$i,1) eq "-" ) 
        {
            $expandedQuality .= "-";
        }
        else
        {
            $expandedQuality .= substr($qual,$posInQuality,1);
            $posInQuality++;
        }
    }
        
    if( length($expandedRead) != length($expandedQuality) )
    {
        errorExit "ERROR Expanded read and quality are of different length";
    }

    return ($expandedRef, $expandedRead, $expandedQuality);
}


=pod

=head2 compressMD

Build the compressed match descriptor out of the read and reference pair

I<Parameters:>

=over 4

=item *

$read_aln

read part of the expanded alignment

=item *

$ref_aln

reference part of the expanded alignment

I<Returns:>

String containing the compressed match descriptor of read against reference

I<Exceptions:>

=over 4

=item *

None

=back

=cut

sub compressMd($$)
{
    my( $read_aln,$ref_aln ) = @_;
    
    return "" if( length($read_aln) != length($ref_aln) );

    my $match = 0;
    my $deleted_bases = 0;
    my $escape_mode = 0;

#    my @read = split( //,$read_aln );
#    my @ref = split( //,$ref_aln );
    my $extMD = "";

#    foreach my $cur_read_el(@read )
    for( my $i=0;$i<length($read_aln);$i++ )
    {
#        my $cur_ref_el = shift(@ref);
        my $cur_read_el = substr($read_aln,$i,1);
        my $cur_ref_el = substr($ref_aln,$i,1);
        
        if( $cur_read_el eq $cur_ref_el  ) # match case
        {
            if( $escape_mode == 1 )
            {
                if( $deleted_bases > 0 )
                {
                    $extMD .= $deleted_bases;
                    $deleted_bases = 0;
                }
                $extMD .= "\$";
                $escape_mode = 0;
            }
            $match++;
        }
        elsif( $cur_ref_el eq "-" )
        {
            # deletion in the reference
            if( $match !=0 )
            {
                $extMD .= $match;
            }
            $match = 0;

            if( $escape_mode == 0 )
            {
                $extMD .= "\^";
                $escape_mode = 1;
            }
            $deleted_bases++;
        }
        elsif( $cur_read_el eq '-' )
        {
            # insertion in the reference
            $extMD .= $match if( $match != 0 );
            $match = 0;
            
            if( $escape_mode == 0 )
            {
                $extMD .= "\^";
                $escape_mode = 1;
            }

            if( $deleted_bases > 0 )
            {
                $extMD .= $deleted_bases;
                $deleted_bases = 0;
            }
            
            $extMD .= $cur_ref_el;
        }
        else
        {
            # mismatch case
            if( $escape_mode == 1 )
            {
                if( $deleted_bases > 0 )
                {
                    $extMD .= $deleted_bases;
                    $deleted_bases = 0;
                }
                $extMD .= "\$";
                $escape_mode = 0;
            }
            
            if( $match != 0 )
            {
                $extMD .= $match;
            }
            $match = 0;

            $extMD .= $cur_ref_el;
        }
    }
    
    # cleaning the reset up
    $extMD .= $match if( $match > 0 );
    $extMD .= $deleted_bases if( $deleted_bases > 0 );
    $extMD .= "\$" if( $escape_mode == 1 );
    

    return $extMD;
}

{ # begining scope of cache
    my %cache;
=pod

=item getExactScoreForRead($p)

The procedure gets score to give to read if all bases correct - a log probability.
It uses internall cache to improve performance.
TC 8.8.7 - cache answers to avoid repeated computation

B<Parameters:>

    $p                      - probability

B<Returns:> 

    score 

=cut
    sub getExactScoreForRead($)
    {
        my ($qv) = @_;
        my $score=0;
        my $b;
        for (my $i=0;$i<length($qv);$i++)
        {
            $b=substr($qv,$i,1);
            unless (defined($cache{$b}))
            {
                ##$cache{$b}=log(convertQualToProbCorrect(getQualFromChar($b)))
                $cache{$b}=log(convertPhredToProbCorrect(getQualFromChar($b)))
            } # unless
            $score+=$cache{$b};
        } # for        
        return $score;
    } # getExactScoreForRead
} # end scope of cache

# sub adjustScoreForIncorrectBase
# Get adjustment to apply to exact match score if a base of this quality 
# matches incorrectly TBD hash this to avoid repeated computation
# TC 8.8.7 - cache answers to avoid repeated computation

{ # begining scope of cache
    my %cache;

=pod

=item adjustScoreForIncorrectBase($p)

The procedure gets adjustment to apply to exact match score if a base of this quality 
 matches incorrectly TBD hash this to avoid repeated computation.
TC 8.8.7 - cache answers to avoid repeated computation

B<Parameters:>

    $c                      - quality

B<Returns:> 

    probability 

=cut
    sub adjustScoreForIncorrectBase($)
    {
        my ($c)=@_;
        my $p;
        unless (defined($cache{$c}))
        {
            ##$p=convertQualToProbCorrect(getQualFromChar($c))
            $p=convertPhredToProbCorrect(getQualFromChar($c));
            $cache{$c}=log((1-$p)/3)-log($p);
        } # unless
        return $cache{$c};
    } # adjustScoreForIncorrectBase
    
} # end scope of cache

=pod

=item getWorstBases($qv, $readLengthELAND)

The procedure gets the score adjustments for the three lowest quality bases that were 
 aligned.
 
B<Parameters:>

    $qv                     - quality
    $readLengthELAND        - read length from ELAND

B<Returns:> 

    probability 

=cut
sub getWorstBases($$)
{
    my ($qv, $readLengthELAND)=@_;
    # print "$qv $readLengthELAND\n";
    $qv=substr($qv,0,$readLengthELAND);
    my @qvals=sort(split(//,$qv));
    # print "@qvals\n";
    for (my $i=0; $i<3; $i++)
    {
        $qvals[$i]=adjustScoreForIncorrectBase($qvals[$i]);
    }
    return ($qvals[0],$qvals[1],$qvals[2]);
} # sub getWorstBases

=pod

=item getMeanBaseQuality($exactScore, $qv)

The procedure gets the mean base quality for a read - uses output of getExactScoreForRead
 aligned.
 
B<Parameters:>

    $exactScore             - score from getExactScoreForRead
    $qv                     - quality

B<Returns:> 

    probability 

=cut
sub getMeanBaseQuality($$)
{
    my ($exactScore, $qv)=@_;
    return convertProbCorrectToPhred(exp($exactScore/length($qv)));
} # sub getMeanBaseQuality

=pod

=item getMeanBaseQuality($exactScore, $qv)

# Work out alignment score (=log(probability)) of alignment given qvals
# Inputs:
# qv=quality value ASCII string eg ^^^^^^^^^^^^^^^^^ZW^V
# type=alignment descriptor, eg 20A4c1A
# thisScore=log-probability score for exact match
# Output:
# Adjusted log-probability score for inexact match
 
B<Parameters:>

    $qv                     - quality
    $type                   - alignment descriptor, eg 20A4c1A
    $thisScore              - score (log-probability score for exact match)

B<Returns:> 

    Adjusted log-probability score for inexact match 

=cut
sub getAlignmentScore($$$)
{ 
    my ($qv,$type,$thisScore) =@_;
    my $posInRead=0;
    
    if( $type !~ /\$/ ) 
    {
        # no indel in match descriptor
        # proceed as usual
        while ($type=~/(\d*)[ACGT]/g)
        {
            $posInRead++;
            $posInRead+=$1 unless ($1 eq "");
            $thisScore+=
                adjustScoreForIncorrectBase(substr($qv,$posInRead-1,1));
        } # while
    } else
    {
        # we have at least one indel: replace the insertions/deletions
        # by X and Y characters

        # for deletions with respect to the reference, we don't have
        # to count the characters to get the correct positions within
        # the quality string
        # 
        # example :
        # REFE AAGGAAA
        # READ AA--AAA
        #
        # For insertions with respect to reference we have to count
        # the number of inserted bases to get to the correct position
        # within the quality string
        print STDERR "INDEL: posInRead = $posInRead\tqv/type/thisScore = $qv/$type/$thisScore\n" if DEBUG;

        my @non_escaped = split(/\^.*\$/,$type);

        print STDERR "# elements in non_escaped = " . $#non_escaped . "\t @non_escaped\n" if DEBUG;

        my $parsed_elements = 0;
        my $expanded_md = "";
        while ($type=~/\^([\w\d]*)\$/g)
        {
            my $escaped = $1;
            $expanded_md .= $non_escaped[$parsed_elements] unless !(defined $non_escaped[$parsed_elements]);
            print STDERR "escaped sequence = $1\tlength = " . length($1) . "\n" if DEBUG;

            # check if $1 is a digit or not, insert X's or Y's accordingly
            if( $escaped =~ /\d+/ )
            {
                $expanded_md .= ( "X" x $escaped );
            } 
                elsif( $escaped =~ /\w+/ )
            {
                $expanded_md .= ( "Y" x length($escaped) );
            } 
            else {
#               print STDERR "found a mixture between digits and characters, this should not happen.\n";
            }

            $parsed_elements++;
        }
        $expanded_md .= $non_escaped[-1] unless !(defined $non_escaped[$parsed_elements]); # we have to append the last non-escaped piece of the MD

#       print STDERR "expanded MD = $expanded_md\n"; # if DEBUG;

        # X are insertions in the read -> count them and increment posInRead
        # Y are insertions in the reference -> ignore them, do not increment posInRead
        while ($expanded_md =~/(\d*)([ACGTXY])/g)
        {
#           print STDERR "character = $2\n" unless $2 eq "";
#           print STDERR "posInRead = $posInRead\n";

            $posInRead++ unless $2 eq "Y";
            $posInRead+=$1 unless ($1 eq "");

#           print STDERR "posInRead = $posInRead\n";


            if( $2 =~ /([ACGT])/ )
            {
                # adjust the score
#               print STDERR "character: $1 - adjust score\tposInRead = $posInRead\n";
                $thisScore+=
                    adjustScoreForIncorrectBase(substr($qv,$posInRead-1,1));
            }
        } # while

    }


    return $thisScore;
} # sub getAlignmentScore

=pod

=item getChromosomeSizes($genomeSizeFile, $nameAttribute)

The procedure reads in XML file of chromosomes and their size, and return a ref
# to a hash. Uses XML::Simple
# Input: file name
# Output: hash ref
 
B<Parameters:>

    $genomeSizeFile         - path to genomesize xml file
    $nameAttribute          - name | fileName. Indiates the xml attribute to use
                              as name for chromosome entries

B<Returns:> 

    HASH MAP Ref with xml as tree 

=cut
sub getChromosomeSizes($$)
{
    my ($genomeSizeFile, $nameAttribute)=@_;
    errorExit "ERROR: invalid argument $nameAttribute. 'contigName' or 'fileName' expected" 
        if (!grep /^$nameAttribute$/, ('contigName', 'fileName'));
    
    # NOTE!!! the following is used instead of XMLin($genomeSizeFile) because
    # XML::Simple fails to read in files that have been given as shell process
    # substitutions such as: 
    # someScriptThatWantsPathToXmlFile.pl --pathToXml=<( someScriptThatPrintsXml.pl )
    open my $genomeSize, "<$genomeSizeFile" or errorExit "Failed to open genomesize file: $genomeSizeFile";
    my $ref=XMLin($genomeSize, ForceArray => 1, KeyAttr => [$nameAttribute]);
    #print Data::Dumper->Dumper($ref);
    close ($genomeSize);
    croak "ERROR: Missing 'chromosome' element in '$genomeSizeFile'. Perhaps an old genome size file?\n"
          unless exists $ref->{chromosome};

    my %ret;
    foreach my $chrom (keys %{$ref->{chromosome}})
    {
        croak "ERROR: Missing 'totalBases' element in '$genomeSizeFile'. Perhaps an old genome size file?\n"
          unless exists $ref->{chromosome}->{$chrom}->{totalBases};
        $ret{$chrom} = $ref->{chromosome}->{$chrom}->{totalBases};
    }
    #print Data::Dumper->Dumper(\%ret);
    return \%ret;
}

=pod

=item getTotalGenomeSize($genomeSizeRef)

Work out total size of genome (ie sequence aligned to) in bp
 
B<Parameters:>

    $genomeSizeRef         - HASH MAP Ref with xml as tree from getChromosomeSizes

B<Returns:> 

    tatal genome size 

=cut
sub getTotalGenomeSize(\%)
{
    my ($genomeSizeRef)=@_;
    my $genomeSize = 0;
    my $correction;
    for my $chromosome (keys %{$genomeSizeRef})
    {
        $genomeSize+=$genomeSizeRef->{$chromosome};
    } # for
    print STDERR "Obtained total genome size of $genomeSize bp\n";
    return $genomeSize;
} # sub getRestOfGenomeCorrection($$)

=pod

=item getRestOfGenomeCorrection($genomeSize, $readLengt)

The procedure creates a constant to add to Bayes Theorem calculations to correct
for misalignments from rest of genome
 
B<Parameters:>

    $genomeSize            - total genome size
    $readLength            - read length

B<Returns:> 

    correction (very small) probability value 

=cut
sub getRestOfGenomeCorrection($$)
{
    my ($genomeSize, $readLength)=@_;
    my $correction = 0;
    $correction = exp(log(2)+log($genomeSize)-(log(4)*$readLength))
                  if $genomeSize;
    print STDERR "Obtained rest-of-genome correction of $correction\n";
    return $correction;
} # sub getRestOfGenomeCorrection($$)

=pod

=item scoreAlignFromNhood($qv, $nhood, $bestScore,$totalScoreSoFar,$readsUsedSoFar,
$readLengthELAND,$thisScore,$restOfGenome)

The procedure estimates confidence for an alignment based on quality values and
neighbourhood information
Need to adjust the neighbourhood to account for the fact we already
have an alignment for the best match and possibly some others too

B<Parameters:>

    $qv                   - quality value ASCII string eg ^^^^^^^^^^^^^^^^^ZW^V
    $nhood                - neighbourhood string eg 1:0:0
    $bestScore            - highest probability out of reads found so far
    $totalScoreSoFa       - summed probabilities for reads considered so far
    $readsUsedSoFa        - number of reads used so far
    $readLength           - read length used in originating ELAND alignment (usually 32)
    $thisScore            - score (log-probability score for exact match)
    $restOfGenome         - correction (very small) probability value

B<Returns:> 

    Estimated prob of alignment being wrong in 'Phred' format, ie Q20=0.99% 

=cut
sub scoreAlignFromNhood($$$$$$$$)
{
    my ($qv, $nhood, $bestScore, 
        $totalScoreSoFar, 
        $readsUsedSoFar,
        $readLengthELAND,
        $thisScore,
        $restOfGenome)=@_;

    print "scoreAFN $nhood $bestScore $totalScoreSoFar $readsUsedSoFar $thisScore\n"
        if DEBUG;

#    if ($nhood eq '1:0:0')
#    { # handle this case separately for speed
#        my $toBase10=10.0/log(10.0);#

#        print "scoreAFN 1:0:0 ", 
#        floor($toBase10*(log($restOfGenome/$bestScore))), "\n" if DEBUG;
#        return(floor($toBase10*(log($restOfGenome/$bestScore))));
#        return(floor($toBase10*(log($restOfGenome)-$thisScore)));
#    } # if
##    $totalScoreSoFar+=$restOfGenome;

    my @nbors;

    if ($nhood=~/(\d+):(\d+):(\d+)/)
    {
        @nbors=($1,$2,$3);
    } # if


    for my $nbor (@nbors)
    {
        while (($readsUsedSoFar>0)&&($nbor>0))
        {
            $readsUsedSoFar--;
            $nbor--;
        } # while
    } # for @nbors

    if ($nbors[0]+$nbors[1]+$nbors[2]>0)
    {
        my @worstBases=getWorstBases($qv,$readLengthELAND);


        # print "Adjusted nhood $nbors[0]:$nbors[1]:$nbors[2]\n";

        # add in some 3-neighbours to get approx score in (1,0,0),(0,1,0) and 
        # (0,0,1) cases - arbitrary value for time being
#        push @nbors,0; 

        unshift (@worstBases,0);

        print "ScoreAFN: $thisScore ", exp($thisScore), "\n" if DEBUG;

        for (my $i=0; $i<3; $i++)
        {
#        print STDERR "$i $nbors[$i] $worstBases[$i] ", exp($thisScore), "\n";
            $thisScore+=$worstBases[$i];
            $totalScoreSoFar+=$nbors[$i]*exp($thisScore);
            print "ScoreAFN: $i $nbors[$i] $thisScore ", 
            exp($thisScore), "\n" if DEBUG;

        } # for
        
    } # if

    print "scoreAFN: $totalScoreSoFar $bestScore\n" if DEBUG;

    my $toBase10=-10.0/log(10.0);
    $bestScore=floor($toBase10*
                     log($totalScoreSoFar/($totalScoreSoFar+$bestScore)));

    print "scoreAFN: $bestScore\n" if DEBUG;

#    $bestScore=$toBase10*log($totalScoreSoFar/($totalScoreSoFar+$bestScore));


#        # print "SCORE: $bestScore $totalScoreSoFar ";
#        if ($bestScore==$totalScoreSoFar)
#        {
#            die; # should never be here
#        $bestScore=$bestPossibleScore;
#        } # if
#        else
#        {
#        my $toBase10=10.0/log(10.0);
#        print STDERR "$qv $restOfGenome $thisScore\n";
#        print STDERR "Would say: ",(floor($toBase10*(log($restOfGenome)-$thisScore))), "\n";
#        $bestScore=floor($toBase10*(log($totalScoreSoFar)-$bestScore));
        

#        $bestScore
#            =floor(convertProbCorrectToPhred($bestScore/$totalScoreSoFar));
#    } # else


#    print STDERR "In fact: $bestScore\n";
    return $bestScore;

} # scoreAlignFromNhood

=pod

=item bestAlignFromList($qv, $nhood, $details, $readLengthELAND, $exactScore,
$restOfGenome)

The procedure picks best alignment from a list of possibilities and work out its 
# relative probability with respect to the rest

B<Parameters:>

    $qv                   - quality value ASCII string eg ^^^^^^^^^^^^^^^^^ZW^V
    $nhood                - neighbourhood string eg 1:0:0
    $details              - highest probability out of reads found so far
    $readLength           - read length used in originating ELAND alignment (usually 32)
    $exactScore           - score from getExactScoreForRead

B<Returns:> 

    Estimated prob of alignment being wrong in 'Phred' format, ie Q20=0.99%

=cut
sub bestAlignFromList($$$$$$$)
{
    my($qv, $nhood, $details, $readLengthELAND, $exactScore,
       $restOfGenome, $nameSource)=@_;
    my $totalScore=$restOfGenome;
    my $bestScore=-100000000;
#    my $exactScore=getExactScoreForRead($qv);
    my $thisScore;
    my ($chrom, $pos, $strand, $type);
    my ($bestChrom, $bestPos, $bestStrand, $bestType);
    my $bestPosInList;

    print "bestAFL $nhood $details $readLengthELAND $exactScore $restOfGenome\n"
        if DEBUG;

    my @s=split(",",$details);
    my @scores;
    $scores[scalar(@s)-1]=0; # size to right size

    for (my $i=0;$i<@s;$i++)
    {
        ($chrom, $pos, $strand, $type)=getMatchDetails($s[$i], $chrom, $nameSource);
        
#        $thisScore=$exactScore;

        $thisScore=getAlignmentScore($qv,$type,$exactScore);
         # print "$s[$i] $chrom $pos $strand $type $thisScore ";

        print "bestAFL $i $chrom $pos $strand $type $thisScore = ",
        exp($thisScore), "\n" if DEBUG;

        $thisScore=exp($thisScore);

         # print "$thisScore\n";
        $scores[$i]=$thisScore;

#        $totalScore+=$thisScore;
        
#        $bestScore=$thisScore-1 unless (defined($bestScore));

        if ($thisScore>$bestScore)
        {
            $bestScore=$thisScore;
            ($bestChrom, $bestPos, $bestStrand, $bestType)
                =($chrom, $pos, $strand, $type);
            $bestPosInList=$i;
        } # if
    } # for

    print "bestAFL $bestPosInList\n" if DEBUG;

    for (my $i=0;$i<@scores;$i++)
    {
        $totalScore+=$scores[$i] unless ($i==$bestPosInList);
    }

    # print "SCORE: $bestScore $totalScore\n";

    $bestScore=scoreAlignFromNhood
            ($qv, $nhood, $bestScore, 
             $totalScore, scalar(@s), $readLengthELAND, $exactScore,
             $restOfGenome );

    # print "SCORE NHOOD-ADJUSTED: $bestScore\n";
    
    return ($bestChrom, $bestPos, $bestStrand, $bestType, $bestScore); 
} # bestAlignFromList

=pod

=item pickBestSingleReadAlignment($qv, $nhood, $details, $readLengthELAND, $exactScore,
$restOfGenom)

Takes in read quality value string and a set of alignment details
and output the details of the best alignment
TC 19.10.7 - pos==undef is the signal that no match is found
pos<=0 can occasionally occur for wraparounds

B<Parameters:>

    $qv                   - quality value ASCII string eg ^^^^^^^^^^^^^^^^^ZW^V
    $nhood                - neighbourhood string eg 1:0:0
    $details              - match details string
    $readLength           - read length used in originating ELAND alignment (usually 32)
    $exactScore           - score from getExactScoreForRead
    $restOfGenome         - rest of genome correction

B<Returns:> 

    ARRAY with 
    chrom=chromosome
    pos=position in chromosome (must be >0, else no match)
    strand=F or R
    type=match description
    score=estimated probability alignment is correct (Phred format, ie Q20=99%)

=cut
sub pickBestSingleReadAlignment($$$$$$$)
{
    my ($qv, $nhood, $details, $readLengthELAND, $exactScore,
        $restOfGenome, $nameSource)=@_;
    my ($chrom, $pos, $strand, $type, $score);

    # score for exact match
#    $exactScore=getExactScoreForRead($qv); 
    print 
    "pickBSRA $nhood $details $readLengthELAND $exactScore $restOfGenome\n" 
     if DEBUG;
    print 
    "pickBSRA $qv\n" 
     if DEBUG;

    if ($details eq '-')
    { # no match found for read
#        print "$readID aaa$nhood 0 F - 0\n";
#        next;
        $pos=undef;
        return ($chrom, $pos, $strand, $type, $score);
    } # if
    elsif ($details!~/,/)
    { # just one match found for read, score it
        ($chrom, $pos, $strand, $type)=getMatchDetails($details, $chrom, $nameSource);
        
        # score for exact match
#        $score=getExactScoreForRead($qv); 
        $score=$exactScore; 
        # adjust to score for actual match
        $score=getAlignmentScore($qv,$type,$score); 
        # convert from score to prob

        print "pBSRA: score converted to $score = ", exp($score),
        " by mismatches\n" if DEBUG;


        $score=exp($score);
        # convert prob to score relative to neighbours


        $score=scoreAlignFromNhood
            ($qv, $nhood, $score, $restOfGenome, 1, $readLengthELAND,
             $exactScore, $restOfGenome );

#        print "$readID $chrom $pos $strand $type $score\n";
    } # if
    else
    { # multiple matches found for read, pick best and score it
        ($chrom, $pos, $strand, $type, $score)
            =bestAlignFromList( $qv, $nhood, $details, $readLengthELAND,
                                $exactScore, $restOfGenome, $nameSource);
#        print "$readID $chrom $pos $strand $type $score\n";
    } # else
    return ($chrom, $pos, $strand, $type, $score);
} # sub pickBestSingleReadAlignment

#------------------------------------------------------------------------------
# returns chrom or contig followed by tab.
sub explodeChrom($$)
{
    $_[0]=~/^(.*?)(\/(.*))?$/;
    my @chromContig=($1,$3);
    my $chrom = $chromContig[$explodeChromIndexMap{$_[1]}];
    $chrom = $chromContig[0] unless defined($chrom);
    return "$chrom\t";
}

sub explodePartnerOffset($$)
{
    my ($chromContigOffset, $nameSource) = @_;
    $chromContigOffset=~/((.*?)(\/(.*))?:)?(\-?\d+)$/;
    my @chromContig = ($2, $4);
    my $offset = $5;
    $offset="" unless defined($offset);
    my $chrom = $chromContig[$explodeChromIndexMap{$nameSource}];
    $chrom = $chromContig[0] unless defined($chrom);
    $chrom = '' unless defined($chrom);
    return "$chrom\t\t$offset";
}
#------------------------------------------------------------------------------


# Uncomment this and run to demo quality value conversion
#for (my $qv=-5;$qv<=60;$qv++)
#{
#    my $p1=convertQualToProbCorrect($qv);
#    my $p2=convertQualToProbWrong($qv);
#    my $qw=convertProbToQual($p2/3); # assume 3 wrong bases equiprobable
#    printf ("%d %10.7f %10.7f %10.7f %10.7f %d %d %10.7f\n", 
#            $qv, $p1, $p2, convertProbToQual($p2), $qw, floor($qw), 
#            floor($qw)-$qv, convertProbCorrectToPhred($p1) );
#} # for

=pod

=head2 SIGPIPE debugging Functions

 set of functions to do SIGPIPE debugging.

=cut
my $sigpipe_count = 0;

=pod

=item sigpipe_handler($)

The procedure sigpipe_handler

B<Parameters:>

    $                  - ???

B<Returns:> 

    Nothing

=cut
sub sigpipe_handler($)
{
    ++$sigpipe_count;
}

=pod

=item conditionally_init_sigpipe_handler()

The procedure conditionally_init_sigpipe_handler

B<Parameters:>

B<Returns:> 

    Nothing

=cut
sub conditionally_init_sigpipe_handler()
{
    my $catch_pipe_sigs = (defined($ENV{'ILMN_CATCH_PIPE_SIGS'})
                           ? $ENV{'ILMN_CATCH_PIPE_SIGS'}
                           : 0);

    if ($catch_pipe_sigs != 0) {
        $SIG{PIPE} = \&sigpipe_handler;
    }
}

=pod

=item conditionally_report_sigpipes($)

The procedure conditionally_report_sigpipes

B<Parameters:>

B<Returns:> 

    Nothing

=cut
sub conditionally_report_sigpipes()
{
    if ($sigpipe_count != 0) {
        print STDERR "Warning: SIGPIPE count : $sigpipe_count";
    }
}

#------------------------------------------------------------------------------


## Not used
#sub trim_external_Ns($$$)
#{
#    my $read_seq = shift;
#    my $trim_start_pos_ref = shift;
#    my $trim_len_ref = shift;

    # To persuade a simple split to do what is wanted, we append an 'N'
    # to the string, do the split and then subtract 1 one from the trailing
    # N count.
#    my @chunks = split(/[^N]/, ($read_seq . 'N'));
#    $$trim_start_pos_ref = length($chunks[0]);
#    my $trail_n_count = length($chunks[$#chunks]) - 1;

#    if (scalar(@chunks) == 1) {
        # This should not happen as an all-Ns read should have failed quality
        # control and not reached here at all but handle the case anyway.
#        $$trim_start_pos_ref = 0;
#        $$trim_len_ref = 0;
#        return '';
#    }

#    $$trim_len_ref = length($read_seq) - $$trim_start_pos_ref - $trail_n_count;

#    return substr($read_seq, $$trim_start_pos_ref, $$trim_len_ref);
#}

=pod

=head2 Run folder path/name toolkit

 set of functions parse and generate paths within run folder

=cut


=pod

=item lane_prefix($lane_prefix)

The procedure generates lane prefix

B<Parameters:>

    $lane_num              - lane number

B<Returns:> 

    lane prefix e.g.: "s_"

=cut
sub lane_prefix($)
{
    my $lane_num = shift;
    return 's_' . $lane_num;
}

=pod

=item lane_num($lane_prefix)

The procedure parses lane prefix

B<Parameters:>

    $lane_prefix            - lane prefix e.g.: "s_"

B<Returns:> 

    lane number

=cut
sub lane_num($)
{
    
    my $lane_prefix = shift;

    if ($lane_prefix =~ /^s_(\d+)$/) {
        return $1;
    }

    return 0;
}

=pod

=item get_tiles_by_lane($tile_file_path, $lane_prefix, $tile_list_ref)

The procedure gets list of all tiles from a give lane.

B<Parameters:>

    $tile_file_path         - tile file path
    $lane_prefix            - lane prefix e.g.: "s_"
    $tile_list_ref          - ARRAY Ref to tile list

B<Returns:> 

    status 0 or 1

=cut
sub get_tiles_by_lane($$\@)
{
    my ($tile_file_path, $lane_prefix, $tile_list_ref) = @_;
    open(TILE_FILE, "<$tile_file_path") or return 0;

    my $orig_size = @$tile_list_ref;
    while (my $line = <TILE_FILE>) {
        next if ($line =~ m/^\#/);

        my @curr_tiles = split(" ", $line);

        foreach my $tile (@curr_tiles) {
            if ($tile =~ /^${lane_prefix}/) {
                push @$tile_list_ref, $tile;
            }
        }
    }

    close(TILE_FILE);

    return ($orig_size != @$tile_list_ref);
}

=pod

=item get_expected_tile_nums($tile_file_path, $lane_num, $expected_tile_nums_ref)

The procedure gets list of all tiles from a give lane.

B<Parameters:>

    $tile_file_path         - tile file path
    $expected_tile_nums_ref - HASH MAP Ref to expected list of tiles
    $lane_num               - line id

B<Returns:> 

    status 0 or 1

=cut
sub get_expected_tile_nums($$\@)
{
    my ( $tile_file_path, $lane_num, $expected_tile_nums_ref) = @_;

    my @tile_list;

    if (!get_tiles_by_lane($tile_file_path, lane_prefix($lane_num),
                           @tile_list)) {
        return 0;
    }

    @$expected_tile_nums_ref
        = map { my @arr = split('_', $_); $arr[-1] } @tile_list;

    return 1;
}

=pod

=item determine_missing_tile_nums($tile_file_path, $tile_freqs_ref, 
$lane_num, $missing_tile_nums_ref)

The procedure compares expected tile list with existing tile list and returns
difference.

B<Parameters:>

    $tile_file_path         - tile file path
    $tile_freqs_ref         - HASH MAP Ref to existing list of tiles
    $lane_num               - line id
    $missing_tile_nums_ref  - HASH MAP Ref to difference list of tiles

B<Returns:> 

    status 0 or 1

=cut
sub determine_missing_tile_nums($\%$\@)
{
    my ( $tile_file_path, $tile_freqs_ref, $lane_num, $missing_tile_nums_ref) = @_;
    
    my @expected_tile_nums = ();

    if (!get_expected_tile_nums($tile_file_path, $lane_num,
                                @expected_tile_nums)) {
        die("Failed to get expected tile nums from file `$tile_file_path'.");
    }

    my @missing_tile_nums = ();

    foreach my $expected_tile_num_str (@expected_tile_nums) {
        # Strip any leading zeros to match keys in %tile_freqs.
        my $expected_tile_num = int($expected_tile_num_str);

        if (!defined($tile_freqs_ref->{$expected_tile_num})) {
            push(@missing_tile_nums, $expected_tile_num);
        }
    }

    @$missing_tile_nums_ref = @missing_tile_nums;
    return 1;
}

=pod

=item get_paths($tile_file_path, $lane_prefix, $read_num, $suffix,
    $path_list_ref)

This procedure constructs a list of file paths of the form s_L_R_TTTT<suffix>
from those tiles found in the specified tile file that have the specified 
lane prefix using the specified read number and suffix. 

B<Parameters:>

    $tile_file_path         - tile file path
    $lane_prefix            - lane prefix
    $read_num               - read number
    $suffix                 - filename suffix
    $path_list_ref          - ARRAY Ref to path list

B<Returns:> 

    status 0 or 1

=cut
sub get_paths($$$$\@)
{
    my ($tile_file_path, $lane_prefix, $read_num, $suffix, $path_list_ref) 
        = @_;

    my @tile_list;

    if (!get_tiles_by_lane($tile_file_path, $lane_prefix, @tile_list)) {
        die("Failed to get tile list from file `$tile_file_path'.");
    }

    foreach my $tile_name (@tile_list) {
        my ($prefix, $lane_num, $tile_str) = split('_', $tile_name);
        push(@{$path_list_ref},
             (join('_', $prefix, $lane_num, $read_num, $tile_str)
              . $suffix)); 
    }

    return 1;
}

sub getEOL($)
{
    my $handle = shift;
    my $firstTwoLines = undef;
    my $eol = undef;

    $firstTwoLines = <$handle>;
    $firstTwoLines .= <$handle>;

    $eol = $1  if ($firstTwoLines =~ /^[^\n\r]+([\n\r]+)[^\n\r]+/);

    seek $handle, 0, 0;
    return $eol;
}



=pod

=item reallyRealPath($path)

This procedure attempts to convert a path provided by the user on the
command line into an absolute path. It should be able to handle "~"
paths and conventional relative paths using ".." or ".". Resolution of
links should follow the convention of "Cwd::realpath".

B<Parameters:>

    $dirRef         - path (converted to absolute path in place)

B<Returns:>

    returns zero if successful, non-zero otherwise.

=cut
sub reallyRealPath(\$) {
    my ($dirRef) = @_;
    my @tmp=glob($$dirRef);
    return 1 if(scalar(@tmp) != 1);
    my $ret = realpath($tmp[0]);
    return 1 if !$ret && !($ret = File::Spec->rel2abs($tmp[0]));
    $$dirRef = $ret;
    return 0;
}



1;    # says use was ok
__END__
=pod

=back

=cut


