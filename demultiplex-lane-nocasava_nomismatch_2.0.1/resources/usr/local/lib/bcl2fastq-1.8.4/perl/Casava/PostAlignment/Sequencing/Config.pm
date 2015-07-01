package Casava::PostAlignment::Sequencing::Config;

# PROJECT: CASAVA
# MODULE:  $RCSfile: Config.pm,v $
# AUTHOR:  Lukasz Szajkowski, Richard Carter
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

=pod

=head1 NAME

Casava::PostAlignment::Sequencing::Config.pm - Perl utility library for configuring BullFrog.

=head1 SYNOPSIS

The library contains procedures and variables usefully in configuring
BullFrog application. It allows to store constants and variables and external
files with format:
[SECTION_NAME]
key spaces value

# include what functions you need... 
use Casava::PostAlignment::Sequencing::Config.pm qw();

=head1 DESCRIPTION

# Global variable
    $globalAppConfFileName default CASAVA.conf
    %CONF_PROJ
    %CONF_APP
    %CONF_RUNS
    %CONF_RUN
    %chrEnds
    %export
    %doubleExport
    %CASAVA
    %runsConfig

=head1 AUTHORSHIP

Copyright (c) 2008-2009 Illumina
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
    @ISA    = qw(Exporter);
    @EXPORT =
      qw( %CONF_PROJ %CONF_APP %CONF_RUNS %CONF_RUN %chrEnds
      &configure &addExportFiles &writeConfig &readConfig &getRunConf &updateRunConf
      &updateRunConfValue &loadConfiguration &configureSampleDirectories %runsConfig
      &configureTargets &updateProjectConfValue &readSampleDirsConfigXML
      &configureDirectories inputBltQC
      &writeSampleDirsConfigXML &isSpliceJunctionChrom &getSpliceJunctionFileName readProjectParameters);
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime ceil);
use XML::Simple;

use File::Spec qw(catdir splitpath);
use Cwd qw(abs_path);
use Sys::Hostname;
use Carp;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use List::Util qw(max);
use IO::File;

use Casava::Common::Log;
use Casava::Common::IOLib
  qw(writeParameters readParameters createDirs executeCmd readXML);
use Casava::TaskManager qw(configureWorkflow);
use Casava::PostAlignment::Sequencing::SortLib qw(checkSortBin);
#use Casava::PostAlignment::Sequencing::DnaSeqLib qw(configureBuild);

sub inputBltQC(\%);
sub inputQC ($$\%\%\%\%$);
sub addExportFiles(\%\@);
sub configure($$\%\%\@\@$);
sub configureDirectories($$\%);
sub configureChromDirectories($$\%\%\@);
#sub configureRuns(\@\@\@\@\%$);
sub configureSampleDirectories(\@\%$);
sub configureTargets($\@\@$);
sub getRunConf($;\%;\%);
sub loadConfiguration($;\%);
#sub readRunsConfig(\%;$);
#sub readRunsConfigXML(\%;$;$);
sub readSampleDirsConfigXML($\@;\@);
sub readConfig(\%;\%; $);
#sub listRuns(\%);
#sub removeRun(\%;$;$);
sub updateRunConf(\%;$;\%;\%);
sub updateRunConfValue($;$);
sub updateProjectConfValue($;$);
#sub repairNonUniqueRunIds(\%);
sub writeRunsConfig(\%;$);
#sub writeRunsConfigXML(\%;$);
sub writeSampleDirsConfigXML(\%;$);
sub writeConfig(\%\%$);
#sub createExperiments ($;$;\%);
sub isSpliceJunctionChrom ($\%);
sub getSpliceJunctionFileName (\%);
our $globalAppConfFileName = "global.conf";

our %CONF_PROJ             = ();
our %CONF_APP              = ();
our %CONF_RUNS             = ();
our %CONF_RUN              = ();
our %chrEnds               = ();
our %runsConfig            = ();

=pod

=head1 Adds export files found in the aligned sample directories.

=over 4

=item addExportFiles($runsConfig, $exportFilesRef)

Adds export files found in the aligned sample directories.

Parameters:
    $runsConfig     - REF to hash map with all sample configuration
    $exportFilesRef - REF to array containing all of the export files

Returns:
    nothing
=back

=cut

sub addExportFiles(\%\@) {
    croak "ERROR: addExportFiles\n" unless ( @_ == 2 );
    my ($runsConfig, $exportFilesRef) = @_;
    
    # evaluate all of the aligned sample directories
    my @filteredExportFiles = ();
    foreach my $alignedSampleDir (@{$runsConfig->{inputDirectories}->{sampleDirectory}}) {

        # grab all of the export filenames
        opendir(DIR, $alignedSampleDir);
        my @dirExportFiles = map { "$alignedSampleDir/$_" } grep(/(L00\d_R\d_\d\d\d_export.txt$|L00\d_R\d_\d\d\d_export.txt.gz$)/, readdir(DIR));
        closedir(DIR);

        if ($CONF_PROJ{lanes})
        {
            my $allowedLanes = join('|', split(/[\s,]/, $CONF_PROJ{lanes}));
            my $grepRegex = "_L00($allowedLanes)_R\\d_\\d\\d\\d_export.txt";
            print STDERR "$grepRegex\n";
            @dirExportFiles = grep /$grepRegex/, @dirExportFiles;
        }
        if ($CONF_PROJ{barcodes})
        {
            my $allowedBarcodes = join('|', split(/[\s,]/, $CONF_PROJ{barcodes}));
            my $grepRegex = "_($allowedBarcodes)_L00\\d_R\\d_\\d\\d\\d_export.txt";
            print STDERR "$grepRegex\n";
            @dirExportFiles = grep /$grepRegex/, @dirExportFiles;
        }
        # add these elements to another array since grep doesn't append
        push(@filteredExportFiles, @dirExportFiles);
    }
    

    # these are POSSIBLE pair xml file names. The actual ones might very well not exist.
    my %pairXmlsHash= ();
    map {my $ori = $_; (s/_R(\d)_(\d\d\d)_export.txt(\.gz)?/_$2_pair.xml/); 
         $pairXmlsHash{$_}->{PairXml}=$_;$pairXmlsHash{$_}->{"Read$1"}=$ori; } @filteredExportFiles;

    foreach my $pairXmlFile (keys %pairXmlsHash)
    {
        my $readMode = "paired";
        $readMode = 'single' if(!-s $pairXmlFile); #assume non-paired data if pair.xml is not avaialble

        if ($CONF_PROJ{read} && 1 == $CONF_PROJ{read}) {
            $pairXmlsHash{$pairXmlFile}->{Read2} = "";
            $readMode        = "single";
        }
        elsif ($CONF_PROJ{read} && 2 == $CONF_PROJ{read}) {
            $pairXmlsHash{$pairXmlFile}->{Read1} = "";
            $readMode        = "single";
        }
        
        $pairXmlsHash{$pairXmlFile}->{PairXml} = "" if 'single' eq $readMode;
        if (!$pairXmlsHash{$pairXmlFile}->{Read1} && !$pairXmlsHash{$pairXmlFile}->{Read2})
        {
            delete $pairXmlsHash{$pairXmlFile}; #remove entries for which data has been masked out
        }
        else
        {
            errorExit("ERROR: Missing read1 export file for $pairXmlFile") if 'paired' eq $readMode && !$pairXmlsHash{$pairXmlFile}->{Read1};
            errorExit("ERROR: Missing read2 export file for $pairXmlFile") if 'paired' eq $readMode && !$pairXmlsHash{$pairXmlFile}->{Read2};
            errorExit("ERROR: Mixed read mode is not allowed. Current project read mode is '$CONF_PROJ{readMode}'. " . 
                      "Will not add $pairXmlsHash{$pairXmlFile}->{Read1} $pairXmlsHash{$pairXmlFile}->{Read2} " . 
                      "with readMode $readMode") if $CONF_PROJ{readMode} && $CONF_PROJ{readMode} ne $readMode;
            
            updateProjectConfValue('readMode', $readMode) if !$CONF_PROJ{readMode};
            push @{$exportFilesRef}, $pairXmlsHash{$pairXmlFile};
        }
        
    }

#    die Data::Dumper::Dumper($exportFilesRef);

}


=pod

=item inputBltQC($CONF_PROJ_Ref)

Checks validity of configuration options for the blt-based alleleCaller

B<Parameters:>

    $CONF_PROJ_Ref - 
    
B<Returns:> 
    nothing
    
=cut

sub inputBltQC (\%) {
    croak "ERROR: init " . join( ',', @_ ) unless ( @_ == 1 );
    my ( $CONF_PROJ_Ref ) = @_;

    # check for required configuration options:
    my @conf_used = qw(
        binSizeBuild
        dirBuildParsed
        QVCutoff
        QVCutoffSingle
        qualityType
        singleScoreForPE
        dirRefSeq
        snpThreshold
        snpThreshold2
        snpMaxRatio
        variantsMDFilterCount
        variantsMDFilterFlank);

    for my $k (@conf_used) {
        if(not defined $CONF_PROJ_Ref->{$k}){
            errorExit "ERROR: can't find required blt alleleCaller configuration option: $k\n";
        }
    }

    if(not (($CONF_PROJ_Ref->{qualityType} eq "Phred64") or
            ($CONF_PROJ_Ref->{qualityType} eq "Solexa64"))){
        errorExit "ERROR: invalid qualityType for blt alleleCaller: $CONF_PROJ_Ref->{qualityType}\n";
    }
}



=pod

=item inputQC($projectDir,    $buildResultsDir, $chrEndsRef,
        $CONF_PROJ_Ref, $runsConfigRef, $CONF_APP_Ref)

The inputQC initialises AlleleCaller

B<Parameters:>

    $projectDir  - 
    $buildResultsDir  - 
    $chrEndsRef  - 
    $CONF_PROJ_Ref  - 
    $runsConfigRef  - 
    $CONF_APP_Ref - 
    
B<Returns:> 
    nothing
    
=cut

sub inputQC ($$\%\%\%\%$) {

    my (
        $projectDir,    $buildResultsDir, $chrEndsRef,
        $CONF_PROJ_Ref, $runsConfigRef,   $CONF_APP_Ref, $checkExportInputs
      )
      = @_;
    my $sortExportCmd = File::Spec->catfile('/usr/local/libexec/bcl2fastq-1.8.4', $CONF_APP_Ref->{cmdExport});
    my $maxNumberOfBins;
    my $error = 0;
    if ( !-d $projectDir ) {
        errorExit
"ERROR: unable to resolve --projectDir=${projectDir} into a folder path. Path to a directory is expected.";
    }

    printLog( sprintf( "CHECK %-40s [%s]\n", "Directory structure", "OK" ), 0 , 'info');

    if ( !-d $CONF_PROJ_Ref->{dirRefSeq} ) {
        errorExit
"ERROR: unable to resolve --refSequences=$CONF_PROJ_Ref->{dirRefSeq} into a folder path. Path to a directory with .fa files is expected.";
    }

    # foreach my $chrom ( sort keys %$chrEndsRef ) {
    #     my $refSeqFile = $chrEndsRef->{$chrom}->{file};

    #     if ( isSpliceJunctionChrom($chrom, %{$CONF_PROJ_Ref}) )
    #     {
    #         next;
    #     }
    #     unless ( -e $refSeqFile ) {
    #         errorExit "ERROR: reference sequence $refSeqFile does not exist.\n"
    #           . "Copy reference sequences to $projectDir/genomes or provide "
    #           . "--refSequences option. \n";
    #     }    # unless
    # }    # foreach

    printLog(
        sprintf( "CHECK %-40s [%s]\n", "Reference sequence exists", "OK" ), 0 , 'info');

    if ($checkExportInputs)
    {
        # Check if export files exist and if lanes not already in build
        my $cmd = "$sortExportCmd --projectDir=$projectDir --target=check";
        executeCmd( $cmd, 5 );
        printLog( sprintf( "CHECK %-40s [%s]\n", "All required input files exist", "OK" ), 0 , 'info');
    }

    # allow CASAVA to be rerun on older builds:
    #
    if(not defined $CONF_PROJ_Ref->{snpThreshold2} ){
        logWarning("Parameter 'snpThreshhold2' is undefined. Setting to value of 'snpThreshold'");
        $CONF_PROJ_Ref->{snpThreshold2} = $CONF_PROJ_Ref->{snpThreshold};
    }

    # check parameters for blt allele-caller:
    # 
    inputBltQC(%$CONF_PROJ_Ref);

    # check that working sort program is found:
    #
    checkSortBin(%$CONF_APP_Ref);

    unless ( $CONF_PROJ_Ref->{qualityType} =~ "Phred64|Solexa64" ) {
        errorExit
          "ERROR: unsupported quality type $CONF_PROJ_Ref->{qualityType}\n";
    }
    printLog( sprintf( "CHECK %-40s [%s]\n", "Parameters validation", "OK" ),
        0 , 'info');
}


sub configureBuild(\@$$)
{
    my ($refAlignSampleDirs, $buildDir, $checkExportInputs) = @_;

    my $reportsCmd  = File::Spec->catfile( '/usr/local/libexec/bcl2fastq-1.8.4', $CONF_APP{cmdReports} );
    my $runConfPath = File::Spec->catfile( $buildDir, $CONF_APP{dirConf}, $CONF_APP{runConfFileName} );
    
    # configure each of the provided project-sample directories
    configureSampleDirectories(@$refAlignSampleDirs, %runsConfig, $runConfPath);

    inputQC( $buildDir, $CONF_PROJ{dirBuildParsed}, %chrEnds, %CONF_PROJ, %runsConfig, %CONF_APP, $checkExportInputs);
    
    printLog(sprintf("CHECK %-40s [%s]\n", "Creating/updating data directories", "IN-PROGRESS"), 0, 'info');
    configureDirectories( $buildDir, $CONF_PROJ{dirBuildParsed}, %chrEnds);
    printLog(sprintf("CHECK %-40s [%s]\n", "Creating/updating data directories", "DONE"), 0, 'info');

    executeCmd( "$reportsCmd --projectDir=$buildDir --target=configure", 4 );
    executeCmd( "$reportsCmd --projectDir=$buildDir --target=main", 4 );
}

=pod

=head1 The procedure configures targets.

=over 4

=item configureTargets($confStatus, $targetRef, $refAlignSampleDirs, $buildDir)

The procedure configures targets. In case $confStatus == 1
(previous run failed) the procedure is trying provide guidance for run restart.
It also confiugres the build folder layout if the target list includes 'all' or 'sort'
targets

Parameters:
    $confStatus         - status returned by configure procedure
    $targetRef          - REF to array of target names
    $refAlignSampleDirs - REF to array of input sample directories
    $buildDir           - output folder

Returns:
    nothing
=back

=cut

sub configureTargets($\@\@$) {
    my ( $confStatus, $targetRef, $refAlignSampleDirs, $buildDir) = @_;
    my @targets = @{$targetRef};

    # if the previous CASAVA run didn't finish
    if ( $confStatus == -1 ) {
        
        my %PREV_RUN_CON =
          %{ getRunConf( $CONF_PROJ{currentRunId}, %CONF_RUNS, %CONF_PROJ ) };
       
        unless ( defined $targets[0] && $targets[0] eq 'allClean' ) {
            my $runContent =
              "----------RUN_$PREV_RUN_CON{_runId}-----------------\n";
            map { $runContent .= "$_\t\t= " . $PREV_RUN_CON{$_} . "\n" }
              keys %PREV_RUN_CON;
            $runContent .=
              "------------------------------------------------\n";
            my $ignoreCmd = "$0 $PREV_RUN_CON{argv} --force\n";
            my $msg =
                "$runContent\nERROR: previous run "
              . "[$PREV_RUN_CON{_runId}] didn't finish\n"
              . "To ignore this warning use:\n$ignoreCmd";
            errorExit($msg);
        }    # unless
    }    # if
    
    if ( $confStatus == 0 ) {
        push (@targets, 'all') if ( !scalar(@targets) );

        configureBuild(@$refAlignSampleDirs, $buildDir, 1) if grep (/^sort$|^all$/, @targets);
        @{$targetRef} = ();
        foreach my $target (@targets)
        {
            if ('all' eq $target)
            {
                my $defaultStagesName = "defaultStages$CONF_PROJ{applicationType}$CONF_PROJ{readMode}";
                errorExit "ERROR: application configuration is missing for $defaultStagesName" unless defined $CONF_APP{$defaultStagesName};
                push @{$targetRef}, split (',', $CONF_APP{$defaultStagesName});
#                if (defined $CONF_PROJ{applicationType})
#                {
#                    if($CONF_PROJ{applicationType} eq 'RNA')
#                    {
#                        push( @{$targetRef}, 'rnaCounts' );
#                    }
#                    elsif ($CONF_PROJ{applicationType} eq 'DNA'
#                        && $CONF_PROJ{readMode} eq 'paired')
#                    {
#                        my $sortPos = grep { $targetRef->[$_] eq 'sort' } 0..$#$targetRef;
#                        if(defined $sortPos) {
#                            splice(@{$targetRef},$sortPos,0,'assembleIndels');
#                        } else {
#                            push( @{$targetRef}, 'assembleIndels' );
#                        }
#                    }
#                }   
            }
            elsif ($target =~ /^no(.*)/)
            {
                @{$targetRef} = grep (!/^$1$/, @{$targetRef});                                 
            }
            else
            {
                push (@{$targetRef}, $target);
            }
        }

        if (scalar( @{$targetRef} ))
        {
            if(${$targetRef}[0] ne 'continue')
            {
                my $sortPos = grep { $targetRef->[$_] eq 'sort' } 0..$#$targetRef;
                if(defined $sortPos) {
                    splice(@{$targetRef},$sortPos,0,'refSeq');
                } else {
                    unshift( @{$targetRef}, 'refSeq' );
                }

                push( @{$targetRef}, 'refSeqClean' );
    
                if (scalar(@$targetRef) != scalar(keys %{{ map { $_ => 1 } @$targetRef}}))
                {
                    errorExit "The target list must contain each target only once. "
                             ."Please check the CASAVA command line. Current list: "
                             . join( ',', @{$targetRef} );
                }
                else
                {
                    logInfo "Will configure following targets:" . join( ',', @{$targetRef}), 0;
                }
            }
        }
        else
        {
            logInfo "No targets to configure", 0;
        }
        updateRunConfValue( 'targets', join( ',', @{$targetRef} ) );
    }    # elsif
}    # sub configureTargets



my $maxSamtoolsRefSize = 2**29;



=pod

=head2 scanDirRefSeq

The procedure creates HASH with chromosome file information by scanning 
dirRefSeq and scanning all unmasked fasta files to find total and known
reference length.

I<Parameters:>

=over 4

=item *

$projectDir

project directory

=item *

$CONF_PROJ_REF

project configuration hash reference

=item *

$chrEndsRef

reference to hash to where to put the results

I<Returns:>

Nothing

I<Exceptions:>

=over 4

=item *

Given reference sequence location has no reference files

=item *

New reference file found. 
Adding more reference files to a conifugred build is not allowed as the 
most common scenario for this to happen is when the user changes the 
--refSequences to a location with the sequence that uses different file 
chromsome naming convention such as chr instead c

=item *

Content of a reference file is different from the last time checked.
It is not allowed to edit reference sequence once it has been used
for the build.

=back

=cut

sub scanDirRefSeq($\%\%) {
    my ( $projectDir, $CONF_PROJ_REF, $chrEndsRef) = @_;

    my $globMask = File::Spec->catfile( $CONF_PROJ_REF->{dirRefSeq}, $CONF_PROJ_REF->{maskRefSeqFiles});
    my @files = glob($globMask);
    errorExit "ERROR: Could not find a single reference in $globMask" if (!scalar(@files));

    my $filesString = join (' ', map {"'$_'"} @files);
    my $countFastaBasesCmd = File::Spec->catfile('/usr/local/libexec/bcl2fastq-1.8.4', 'countFastaBases') . " $filesString";

    my @contigInfos = `$countFastaBasesCmd`;
    errorExit "ERROR: The following command failed with exit code $?: $countFastaBasesCmd" unless !$?;

    my @scaffolds;
    my $isContigName = ('contigName' eq $CONF_PROJ_REF->{chromNameSource});

    my %scannedContigs;
    foreach my $contigInfo ( @contigInfos ) {
        my ($fileName, $contigName, $knownBases, $totalBases) = split /\s/, $contigInfo;
        errorExit "ERROR: duplicate reference sequence name $contigName in $fileName is prohibited when --chromNameSource=$CONF_PROJ_REF->{chromNameSource}"
            if (exists $scannedContigs{$contigName});
        $scannedContigs{$contigName} = $fileName;

        my $scaffold = ($isContigName ? $contigName : $fileName );
        if ($isContigName && 
            'on' eq $CONF_PROJ_REF->{chromNameValidation} && 
            $scaffold =~ /[$CONF_APP{contigNameForbiddenCharacters}]/) {
            my $msgTxt;
            eval "\$msgTxt = \"ERROR: Chromosome name cannot contain the following characters: $CONF_APP{contigNameForbiddenCharacters}\"";
            errorExit $msgTxt;
        }
        push @scaffolds, $scaffold;
        my $scaffoldFile = File::Spec->catfile( $CONF_PROJ_REF->{dirRefSeq}, $fileName);
        my $singleFile = ($isContigName ? File::Spec->catfile($projectDir, $CONF_PROJ_REF->{refSeqCacheDirName}, $contigName): $scaffoldFile);
        $chrEndsRef->{$scaffold}->{file} = $singleFile;

        my $projConfKnownBasesKey = "knownBases_$scaffold";
        my $projConfTotalBasesKey = "totalBases_$scaffold";
        my $projConfTimeStampKey = "timeStamp_$scaffold";
        my $scaffoldTs = (stat($scaffoldFile))[9];

        # this is not the first time we're scanning the folder and we don't have timestamp for this file
        errorExit "ERROR: unexpected new reference sequence $scaffold found in $CONF_PROJ_REF->{dirRefSeq}" 
            if ($CONF_PROJ_REF->{refSeqLocked} && !defined $CONF_PROJ_REF->{$projConfKnownBasesKey});

        if (!defined $CONF_PROJ_REF->{$projConfTimeStampKey} ||
            $scaffoldTs != $CONF_PROJ_REF->{$projConfTimeStampKey})
        {
            if (defined $CONF_PROJ_REF->{$projConfTimeStampKey} &&
                ($totalBases != $CONF_PROJ_REF->{$projConfTotalBasesKey} ||
                 $knownBases != $CONF_PROJ_REF->{$projConfKnownBasesKey})) {
                errorExit "ERROR: The contents of the reference sequence file changed. "
                         ."Expected: knownBases=$CONF_PROJ_REF->{$projConfKnownBasesKey}, totalBases=$CONF_PROJ_REF->{$projConfTotalBasesKey}"            
                         ."Actual: knownBases=$knownBases, totalBases=$totalBases. ";
            }
            $CONF_PROJ{$projConfKnownBasesKey} = $knownBases;
            $CONF_PROJ{$projConfTotalBasesKey} = $totalBases;
            $CONF_PROJ{$projConfTimeStampKey} = $scaffoldTs;
            logInfo "Reference $scaffold has $knownBases known out of $totalBases bases total", 0;
        }

        if($CONF_PROJ_REF->{$projConfTotalBasesKey} > $maxSamtoolsRefSize) {
            logWarning("Reference $scaffold has length " . $CONF_PROJ{$projConfKnownBasesKey} . ", which exceeds the maximum size of the BAM file index: $maxSamtoolsRefSize");
        }

        $chrEndsRef->{$scaffold}->{knownBases} = $CONF_PROJ_REF->{$projConfKnownBasesKey};
        $chrEndsRef->{$scaffold}->{length} = $CONF_PROJ_REF->{$projConfTotalBasesKey};
    }

    $CONF_PROJ{chromosomes} = join("\t", @scaffolds);

    $CONF_PROJ_REF->{refSeqLocked} = 'yes';
}



=pod

=head2 getChrEnds

The procedure fills in the chrEnds HASH based on configuration data.

I<Parameters:>

=over 4

=item *

$projectDir

project directory

=item *

$CONF_PROJ_REF

project configuration hash reference

=item *

$chrEndsRef

reference to hash to where to put the results

I<Returns:>

Nothing

=cut

sub getChrEnds($\%\%) {
    my ( $projectDir, $CONF_PROJ_REF, $chrEndsRef) = @_;

    my $filePath;
    if('contigName' eq $CONF_PROJ_REF->{chromNameSource}) {
        $filePath = File::Spec->catdir($projectDir, $CONF_PROJ_REF->{refSeqCacheDirName});
    } else {
        $filePath = $CONF_PROJ_REF->{dirRefSeq};
    }

    unless(defined($CONF_PROJ_REF->{chromosomes})){
        errorExit("ERROR: Attempting to run with incompatible configuration file\n");
    }

    my $chromString = $CONF_PROJ_REF->{chromosomes};
    my @chroms = split("\t",$chromString);

    for my $chrom (@chroms) {
        my $projConfKnownBasesKey = "knownBases_$chrom";
        my $projConfTotalBasesKey = "totalBases_$chrom";
        $chrEndsRef->{$chrom}->{knownBases} = $CONF_PROJ_REF->{$projConfKnownBasesKey};
        $chrEndsRef->{$chrom}->{length} = $CONF_PROJ_REF->{$projConfTotalBasesKey};
        $chrEndsRef->{$chrom}->{file} = File::Spec->catfile($filePath,$chrom);
    }
}



=pod

=head1 The procedure loads the application's config.

=over 4

=item loadConfiguration($buildSampleDir; \%PARAMS)

The procedure loads the application's config from $buildSampleDir/conf/. 
Should be run at the beginning of every script which want to use the
configuration generated by configure procedure.

Parameters:
    $buildSampleDir - path to project directory
    $PARAM          - optional hash of properties that override the ones loaded

Returns:
    nothing
=back

=cut

sub loadConfiguration($;\%) {
    my ($buildSampleDir, $PARAMS_Ref) = @_;

    my $logPath = File::Spec->catfile( $buildSampleDir, "CASAVA.log" );
    initLog($logPath, 0,4);

    my $globalConf  = File::Spec->catfile( '/usr/local/etc/bcl2fastq-1.8.4', $globalAppConfFileName );

    readParameters( %CONF_APP, "APP_CONFIGURATION", $globalConf );
    my $projectConf = File::Spec->catfile( $buildSampleDir, $CONF_APP{dirConf}, $CONF_APP{projectConf} );
    errorExit "ERROR: loadConfiguration() CASAVA configuration [$projectConf] cannot be found" 
        unless ( -e $projectConf );
    readConfig( %CONF_PROJ, %CONF_RUNS, $projectConf );

    my $currentAppMajorMinor = 'bcl2fastq-1.8.4';
    $currentAppMajorMinor =~ s/(.*)\..*$/$1/;
    
    my $buildAppMajorMinor = $CONF_PROJ{appVersion};
    $buildAppMajorMinor =~ s/(.*)\..*$/$1/;
    
    if ($buildAppMajorMinor ne $currentAppMajorMinor) {
        errorExit "ERROR: Project was produced by an incompatible version of CASAVA. "
                 ."Expected: $currentAppMajorMinor.x. Actual: $CONF_PROJ{appVersion}."
    }

    errorExit "ERROR: applicationType cannot be changed for a configured project." 
        if defined $PARAMS_Ref && $PARAMS_Ref->{applicationType} 
        && $CONF_PROJ{applicationType} ne $PARAMS_Ref->{applicationType};

    if (defined $PARAMS_Ref) {
        foreach my $k ( keys %$PARAMS_Ref ) {
            $CONF_PROJ{$k} = $PARAMS_Ref->{$k} if ( defined $PARAMS_Ref->{$k} );
        }
    }
    $CONF_PROJ{dirRefSeq} = abs_path($CONF_PROJ{dirRefSeq}) if $CONF_PROJ{dirRefSeq};

    getChrEnds( $buildSampleDir, %CONF_PROJ, %chrEnds);

    my $runConf = File::Spec->catfile( $buildSampleDir, $CONF_APP{dirConf}, $CONF_APP{runConfFileName} );

    my $runsConfigRef = \%runsConfig;
    readSampleDirsConfigXML($runConf,
                            @{$runsConfigRef->{inputDirectories}->{sampleDirectory}}, 
                            @{$runsConfigRef->{exportFiles}}) if ( -s $runConf );

    %CONF_RUN = %{ getRunConf( $CONF_PROJ{currentRunId}, %CONF_RUNS, %CONF_PROJ ) };
    configureWorkflow( $CONF_PROJ{isWorkflow}, $CONF_PROJ{workflowFile} );
    initLog($logPath, $CONF_PROJ{verbose},4);    
}


=pod

=head1 The procedure reads HASH MAP from the conf/project.conf file.

=over 4

=item readProjectParameters( $parametersRef, $sectionName, $prorjectDir)

The procedure reads HASH MAP with parameters from the parameters file.

Parameters:
    $parametersRef  - HASH MAP with parameters
    $sectionName    - parameters section name in the file
    $prorjectDir    - location of conf/project.conf

Returns:
    nothing
=back

=cut

sub readProjectParameters(\%$$) 
{
    my ( $parametersRef, $sectionName, $prorjectDir ) = @_;

    my $projectConf = File::Spec->catfile( $prorjectDir, $CONF_APP{dirConf}, $CONF_APP{projectConf} );
    readParameters(%$parametersRef, $sectionName, $projectConf);
}

=pod

=head1 The procedure configures the application.

=over 4

=item configure

The procedure configures the application. First time CASAVA is run for given
$buildSampleDir the $buildSampleDir/conf/Project.conf file is generated based on
etc/global.conf and current command line parameters. Then each time CASAVA
is run for given projectDir the run configuration is updated in
$buildSampleDir/conf/Project.conf

Parameters:
    $buildSampleDir     - project directory
    $argvStr        - concatenated command line parameters
    $CUR_RUN_Ref    - current run config REF
    $PARAMS_Ref     - other parameters from the command line
    $targetsRef     - 
    $readModeRef    - 
    $usage          - 

Returns:
    0 - success; -1 - previous run problem
=back

=cut
## TODO: Refactor !!!
sub configure($$\%\%\@\@$) {

    my ( $buildDir, $argvStr, $CUR_RUN_Ref, $PARAMS_Ref,
        $targetsRef, $readModeRef, $usage ) = @_;

    my $globalConf = File::Spec->catfile( '/usr/local/etc/bcl2fastq-1.8.4', $globalAppConfFileName );
    readParameters( %CONF_APP, "APP_CONFIGURATION", $globalConf );
    
    unless ( -d $buildDir ) {
        mkdir $buildDir, 0775
          or errorExit "ERROR: createDirs() can't make $buildDir $!\n";
    }    # unless
    
    my %CURRENT_RUN_CONF = %{$CUR_RUN_Ref};
    my %PARAMS           = %{$PARAMS_Ref};    

    $CONF_PROJ{dirBuild} = $buildDir if ( !$PARAMS{dirBuild} );
    $CONF_PROJ{dirBuild} = abs_path( $CONF_PROJ{dirBuild} );
    
    if ( scalar(@$targetsRef) == 0 ) {
        $CURRENT_RUN_CONF{target} = 'all';
    } elsif ( scalar(@$targetsRef) == 1 ) {
        $CURRENT_RUN_CONF{target} = $targetsRef->[0];
    }
    
    # check the forced read options
    if($PARAMS{forceRead1} && $PARAMS{forceRead2}) {
        errorExit "ERROR: Both --forceRead1 and --forceRead2 were specified on the command line. Please use choose only one of the options.\n\n$usage";
    }
    
    $CURRENT_RUN_CONF{argv} = $argvStr;
    my $isCreate    = 0;
    my $fullConfDir = File::Spec->catdir( $buildDir, $CONF_APP{dirConf} );
    my $projectConf = File::Spec->catfile( $fullConfDir, $CONF_APP{projectConf} );

    if ( $PARAMS{samtoolsRefFile} ) {
        errorExit "ERROR: Ambiguous, both --refSequences and --samtoolsRefFile are provided" if $PARAMS{dirRefSeq};
        $PARAMS{samtoolsRefFile} = abs_path($PARAMS{samtoolsRefFile});
        my @pathParts = File::Spec->splitpath($PARAMS{samtoolsRefFile});
        $PARAMS{dirRefSeq} = File::Spec->catpath($pathParts[0], $pathParts[1]); 
        $PARAMS{maskRefSeqFiles} = $pathParts[2];
        $PARAMS{chromNameSource} = 'contigName'; 
    }
    
    unless ( -e $projectConf) {
        # Create configuration
        errorExit "ERROR: $0:configure undefined PARAMS_Ref" unless ( defined($PARAMS_Ref) );
        $isCreate = 1;
        my @dirs = ( $buildDir, $fullConfDir );
        createDirs( "", @dirs );
        readConfig( %CONF_PROJ, %CONF_RUNS, $globalConf );
        $CONF_PROJ{buildId} = Digest::MD5::md5_hex( rand() . localtime() . hostname() );
        $CONF_PROJ{dirBuild} = $buildDir;
        foreach my $keyTmp ( keys %PARAMS ) {
            $CONF_PROJ{$keyTmp} = $PARAMS{$keyTmp} if ( defined $PARAMS{$keyTmp} );
        }

        $CONF_PROJ{dirRefSeq} = abs_path($CONF_PROJ{dirRefSeq}) if $CONF_PROJ{dirRefSeq};
        scanDirRefSeq( $buildDir, %CONF_PROJ, %chrEnds);
    }
    else {
        # load existing and update with %PARAMS
        readParameters( %CONF_PROJ, "PROJECT_CONFIGURATION", $globalConf );
        loadConfiguration($buildDir, %PARAMS);
    }
    
    my $timeStr = strftime $CONF_APP{formatTimeStamp}, localtime;

    if ( !defined $CONF_PROJ{dirRefSeq} && $CURRENT_RUN_CONF{target} !~ /^clean/) {
        errorExit "ERROR: $0: configure -parameter --refSequences not defined\n";
    }

    if ( !defined $PARAMS{workflowFile} ) {
        $CONF_PROJ{workflowFile} = $buildDir . "/tasks.$timeStr.txt";
    }
    $CURRENT_RUN_CONF{'startTime'} = $timeStr;

    ## For CASAVA RNA always calculate exon counts and add spliceJunction as scaffold
    if ( $CONF_PROJ{applicationType} eq 'RNA' )
    {
        if ( $CONF_PROJ{rnaCountMethod} ne 'readBases' )
        {
            errorExit "unsuported value for rnaCountMethod=[$CONF_PROJ{rnaCountMethod}]\n";
        }

        $CONF_PROJ{isDenseAlleleCalls} = 1;
    }

    if ( $CONF_PROJ{isWorkflowAuto} == 1 || $CONF_PROJ{isWorkflowSGE} == 1 || $CONF_PROJ{isWorkflowMake} == 1) {
        $CONF_PROJ{isWorkflow} = 1;
    }
    
    unless ( defined $CONF_PROJ{dirBuildParsed} ) {
        $CONF_PROJ{dirBuildParsed} = File::Spec->catdir( $CONF_PROJ{dirBuild}, 'Parsed');
        $CONF_PROJ{dirBuildParsed} =
          File::Spec->catdir( $CONF_PROJ{dirBuild}, "Parsed_" . ( strftime $CONF_APP{formatDate}, localtime ))
            unless $CONF_PROJ{skipVariableMetadata};
    }
    
    if (   $isCreate == 1
        || !exists $CONF_PROJ{configured}
        || $CONF_PROJ{configured} == 0 )
    {
        $CONF_PROJ{dirBuildExport} = File::Spec->catdir( $buildDir, $CONF_APP{dirDup} );
        $CONF_PROJ{projectId}      = $timeStr;
        $CONF_PROJ{previousRunId}  = "";
        $CONF_PROJ{configured}     = 1;
    }
    else {
        $CONF_PROJ{previousRunId} = $CONF_PROJ{currentRunId};

        if ( defined $CURRENT_RUN_CONF{target}
            && $CURRENT_RUN_CONF{target} eq 'continue' )
        {
            writeConfig( %CONF_PROJ, %CONF_RUNS, $projectConf );
            return 0;
        }
        elsif (   defined $CONF_RUN{currentStatus}
            && $CONF_RUN{currentStatus} ne $CONF_APP{stateFinished}
            && $CURRENT_RUN_CONF{force} != 1 )
        {
            return -1;
        }    # Check if target == continue
    }

    if (defined $CONF_PROJ{dirBuildTemp})
    {
        $CONF_PROJ{dirBuildTemp} =~ s/^\$buildDir/$buildDir/;
        printLog "Using $CONF_PROJ{dirBuildTemp} for scratchpad\n", 0;
        unless ( -d $CONF_PROJ{dirBuildTemp} ) {
            mkdir $CONF_PROJ{dirBuildTemp}, 0775
              or errorExit "ERROR: createDirs() can't make $buildDir $!\n";
        }
    }
    
    if( $CONF_PROJ{readSampleRateStart} > $CONF_PROJ{readSampleRateEnd} ) {
        errorExit "ERROR: readSampleRateStart must be = readSampleRateEnd\n";
    }

    if( defined $CONF_PROJ{readSampleRateInput} and 
        ($CONF_PROJ{readSampleRateInput} <= 0.0 or $CONF_PROJ{readSampleRateInput} > 1.0) ) {
        errorExit "ERROR: readSampleRateInput must be in (0;1] range \n";
    }
    
    configureWorkflow( $CONF_PROJ{isWorkflow}, $CONF_PROJ{workflowFile} );
    updateRunConf( %CURRENT_RUN_CONF, $timeStr, %CONF_RUNS, %CONF_PROJ );
    writeConfig( %CONF_PROJ, %CONF_RUNS, $projectConf );
    %CONF_RUN = %{ getRunConf( $CONF_PROJ{currentRunId}, %CONF_RUNS, %CONF_PROJ ) };

    return 0;
}

=pod

=head1 The procedure reads HASH MAP with config from the config file.

=over 4

=item readConfig($confProjRef, $confRunsRef, $fileName)

The procedure reads HASH MAPs with project configuration,
and all runs configuration name from the config file.

Parameters:
    $confProjRef - project configuration
    $confRunsRef - all runs configuration
    $fileName    - config file name

Returns:
    nothing
=back

=cut

sub readConfig(\%;\%; $) {
    croak "ERROR: readConfig\n" unless ( @_ == 3 );
    my ( $confProjRef, $confRunsRef, $fileName ) = @_;
    my $sectionName    = "";
    my $currentConfRef = "";
    if ( -z $fileName ) {
        errorExit "ERROR: readConfig() empty file $fileName\n";
    }
    open( FILE, "<$fileName" )
      || errorExit "ERROR: readConfig() Couldn't open file handle for $fileName $!\n";
    while (<FILE>) {
        my $line = $_;
        chomp $line;

        #        print $line . " ";
        my $comment = "";
        my $key     = "";
        my $value   = "";
        if ( $line =~ /(#.+)$/ ) {
            $comment = $1;
        }    # if
        elsif ( $line =~ /^\[(\S+)\]$/ ) {
            $sectionName = $1;
            if ( $sectionName eq "PROJECT_CONFIGURATION" ) {
                $currentConfRef = $confProjRef;
            }    # if
            elsif ( $sectionName =~ /^RUN_(.+)$/ ) {
                my $runId        = $1;
                my %NEW_RUN_CONF = ();

                #print "[$runId = " . $runId . "]\n";
                $NEW_RUN_CONF{_runId} = $runId;
                $currentConfRef = \%NEW_RUN_CONF;
                ${$confRunsRef}{$runId} = $currentConfRef;
            }    # elsif
            else {
                $currentConfRef = "";
            }
        }    # elsif
        elsif ( $line =~ /^(\S+)\s*(.*)$/ ) {    # if
            if ( $currentConfRef ne "" ) {
                $key   = $1;
                $value = $2;
                ${$currentConfRef}{$key} = $value;
            }                                    # if
        }    # elsif
    }    # while
    close FILE;
}

=pod

=head1 The procedure writes HASH MAP with config to the config file.

=over 4

=item writeConfig($confRef, $confRunsRef, $fileName)

The procedure writes HASH MAPs with project configuration,
and all runs configuration to the config file.

Parameters:
    $confProjRef - project configuration
    $confRunsRef - all runs configuration
    $fileName    - config file name

Returns:
    nothing
=back

=cut

sub writeConfig(\%\%$) {
    my ( $confRef, $confRunsRef, $fileName ) = @_;
    my $content = "";
    my $sectionNameTmp;
    my $overwrite = 0;
    my $found     = 0;
    my $padWidth = 32;
    my $contentProject .= "[PROJECT_CONFIGURATION]\n";
    foreach my $key ( sort keys %{$confRef} ) {
        my $padding = max($padWidth, length($key) + 1);
        my $value = "";
        if ( defined( ${$confRef}{$key} ) ) {
            $value = ${$confRef}{$key};
        }    # if
        $contentProject .= ( sprintf "%-$padding" . 's', $key ) . "$value\n";
    }    # foreach
    $contentProject .= "\n";
    my $contentRuns = "";
    foreach my $runId ( sort { $b cmp $a } keys %{$confRunsRef} ) {
        my $value = "";
        $contentRuns .= "[RUN_$runId]\n";
        my %RUN_CONFIG = %{ ${$confRunsRef}{$runId} };
        foreach my $key ( sort keys %RUN_CONFIG ) {
            my $padding = max($padWidth, length($key) + 1);
            my $value = "";
            if ( defined( $RUN_CONFIG{$key} ) ) {
                $value = $RUN_CONFIG{$key};
            }    # if
            $contentRuns .= ( sprintf "%-$padding" . 's', $key ) . "$value\n";
        }    # foreach
        $contentRuns .= "\n";
    }    # foreach
    $contentRuns .= "\n";
    if ( -e $fileName ) {    # Merge with file
        open( FILE, "<$fileName" )
          || errorExit
"ERROR: writeConfig() $0::Couldn't open file handle for $fileName $!\n";
        while (<FILE>) {
            my $line = $_;
            chomp $line;

            #        print $line . " ";
            if ( $line =~ /^\[(\S+)\]$/ ) {
                $sectionNameTmp = $1;
                if ( $sectionNameTmp eq 'PROJECT_CONFIGURATION' ) {
                    $overwrite = 1;
                    $found     = 1;
                    $content .= $contentProject . $contentRuns . "\n";
                }    # if
                elsif (/^\[RUN_(\S+)\]$/) {
                }    # elsif
                elsif ( $overwrite == 1 ) {
                    $overwrite = 0;
                }    # elsif
            }    # if
            if ( $overwrite == 0 ) {
                $content .= $line . "\n";
            }    # if
        }    # while
        close FILE;
    }    # if
    if ( $found == 0 ) {    # Add to the file (new section)
        $content .=
            "\n# This section contains CASAVA"
          . " parameters for PROJECT_CONFIGURATION\n\n";
        $content .= $contentProject . $contentRuns . "\n";
    }    # if
    my $file = IO::File->new( ">${fileName}.tmp" )
      || errorExit "ERROR: writeConfig Couldn't create file handle for ${fileName}.tmp $!\n";
    print $file $content;
    close $file;
    rename ("${fileName}.tmp", $fileName) || errorExit "ERROR: Cannot rename ${fileName}.tmp into $fileName. $!";

    #exit(0);
}

=pod

=head1 The procedure configures run id and run directories stored in project.conf.

=over 4

=item configureRuns($runIdsRef, $runDirsRef, $runLanesRef, $confRef, $fileName )

The procedure configures run id and run directories stored in project.conf.
 - produture adds or updates runs 

Parameters:
    $runIdsRef    - reference to ARRAY of run ids 
    $runDirsRef   - reference to ARRAY of run directories
    $runLanesRef  - reference to ARRAY of run lanes in each run
    $confRef      - reference to hash map with all runs configuration 
    $fileName     - config file name
       
Returns:
    nothing
    
=back

=cut

#sub configureRuns(\@\@\@\@\%$) {
#    my (
#        $runIdsRef,   $runDirsRef, $runLanesRef,
#        $readModeRef, $confRef,    $fileName
#      )
#      = @_;
#    my @runIds      = @{$runIdsRef};
#    my @runDirs     = @{$runDirsRef};
#    my @runLanes    = @{$runLanesRef};
#    my @readMode    = @{$readModeRef};
#    my $newConfig   = 1;
#    my $defaultMode = "paired";
#    if ( -e $fileName ) {
#
#        #print("Reading run configuration from $fileName\n");
#        readRunsConfigXML( %{$confRef}, $fileName, 0 );
#        $newConfig = 1;
#    }
#    if ( scalar(@readMode) == 0 ) {
#        $defaultMode = "paired";
#    }
#    elsif ( scalar(@readMode) == 1 ) {
#        if ( $readMode[0] ne "paired" && $readMode[0] ne "single" ) {
#            errorExit "ERROR:configureRuns --readMode [$readMode[0]] can be only paired or single\n";
#        }
#        $defaultMode = $readMode[0];
#    }
#    else {
#        my $str = join ",", @readMode;
#        errorExit "ERROR:configureRuns readMode ($str) can be only either paired or single\n";
#    }
#    if (   scalar(@runIds) != scalar(@runDirs)
#        || scalar(@runDirs) != scalar(@runLanes) )
#    {
#        errorExit
#          "ERROR: configureRurepairNonUniqueRunIdsns() inconsistent number of "
#          . "run ids, run directories and run lanes. "
#          . "Check command line parameters\n";
#    }    # if
#
#    for ( my $i = 0 ; $i < scalar(@runIds) ; $i++ ) {
#        my $runId    = $runIds[$i];
#        my $found    = 0;
#        my $setIndex = 1;
#        if ( exists $confRef->{run}->{$runId} ) {
#            my %sets = %{ $confRef->{run}->{$runId}->{set} };
#            foreach my $setId ( sort keys %sets ) {
#                my $runDir = $runsConfig{run}->{$runId}->{set}->{$setId}->{gerald};
#                if ( abs_path($runDir) eq abs_path( $runDirs[$i] ) )
#                {    # update lanes
#                    $found = 1;
#                    my @lanesTemp = split( ',', $runLanes[$i] );
#                    $runsConfig{run}->{$runId}->{set}->{$setId}->{lanes} = \@lanesTemp;
#                    $runsConfig{run}->{$runId}->{set}->{$setId}->{readMode} = $defaultMode;
#                }    # if
#                $setIndex++;
#            }    # foreach
#        }    # if
#        else {
#            $confRef->{run}->{$runId}->{id} = $runId;
#        }    # else
#        if ( $found == 0 ) {
#            my @lanesTemp = split( ',', $runLanes[$i] );
#            $runsConfig{run}->{$runId}->{set}->{$setIndex}->{lanes} = \@lanesTemp;
#            $runsConfig{run}->{$runId}->{set}->{$setIndex}->{gerald} = abs_path( $runDirs[$i] );
#            $runsConfig{run}->{$runId}->{set}->{$setIndex}->{readMode} = $defaultMode;
#        }    # if
#    }    # for
#
#    repairNonUniqueRunIds( %{$confRef} );
#
#    writeRunsConfigXML( %{$confRef}, $fileName );
#}

=pod

=head1 The procedure configures sample directories stored in run.conf.xml

=over 4

=item configureSampleDirectories($alignSampleDirsRef, $sampleDirRef, $fileName)

The procedure configures sample directories stored in run.conf.xml.

Parameters:
    $alignSampleDirsRef - reference to ARRAY of sample directories
    $sampleDirRef       - reference to hash map with sample directories
    $fileName           - config file name

Returns:
    nothing

=back

=cut

sub configureSampleDirectories(\@\%$) {

    # initialization
    my ($alignSampleDirsRef, $sampleDirRef, $xmlFilename) = @_;
    my @alignSampleDirs = @{$alignSampleDirsRef};
    
    # convert the relative paths to absolute paths
    for my $i (0 .. $#alignSampleDirs) {
        my $absDir = abs_path($alignSampleDirs[$i]);
        errorExit("ERROR: The following sample directory does not exist: " . $alignSampleDirs[$i] . ". Please update your --inSampleDir parameter.") unless ($absDir && -e $absDir);
        $alignSampleDirs[$i] = $absDir;
    }
    
    # read in the configuration file if it exists
    if ( -e $xmlFilename ) {
        my @xmlDirs = ();
        readSampleDirsConfigXML($xmlFilename, @xmlDirs);
        push (@alignSampleDirs, @xmlDirs);
    }
    
    # throw everything into a set so that we can remove duplicate sample
    # directories
    my %dirSet = map { $_ => 1 } @alignSampleDirs;
    
    # convert the aligned sample directories to the absolute path
    my $dirIndex = 0;
    foreach my $dir (sort keys %dirSet) {
        $sampleDirRef->{inputDirectories}->{sampleDirectory}->[$dirIndex] = $dir;
        $dirIndex++;
    }

    addExportFiles(%$sampleDirRef, @{$sampleDirRef->{exportFiles}});
    
    # populate the sample directory map and serialize it
    writeSampleDirsConfigXML(%{$sampleDirRef}, $xmlFilename);
}

=pod

=head1 The procedure lists all runs (and lanes) in the configuration.

=over 4

=item listRuns($confRef)

The procedure lists all runs (and lanes) in the configuration.

Parameters:
    $confRef      - reference to hash map with all runs configuration 
       
Returns:
    nothing
    
=back

=cut

#sub listRuns(\%) {
#    croak "ERROR: listRuns\n" unless ( @_ == 1 );
#    my ($confRef) = @_;
#    my $ll =
#      "---------------------------------------------------------------\n";
#    printLog( $ll, 0 );
#    my $content = sprintf( "%-40s%-5s %-20s", 'Run Id', '| Set Id', '| lanes' );
#    printLog( "$content\n", 0 );
#    printLog( $ll, 0 );
#    for my $runIdTmp ( sort keys %{ $runsConfig{run} } ) {
#        my %sets = %{ $runsConfig{run}->{$runIdTmp}->{set} };
#        foreach my $setIdTmp ( keys %sets ) {
#            my @lanes =
#              @{ $runsConfig{run}->{$runIdTmp}->{set}->{$setIdTmp}->{lanes} };
#            my $geraldFullDir =
#              $runsConfig{run}->{$runIdTmp}->{set}->{$setIdTmp}->{gerald};
#            my $lanesStr = join( ',', @lanes );
#            my $content = sprintf( "%-40s| set%-4s| %-20s",
#                $runIdTmp, $setIdTmp, $lanesStr );
#            printLog( "$content\n", 0 );
#        };    # foreach
#    };    # foreach
#    printLog( $ll, 0 );
#}

=pod

=head1 The procedure removes run from run configuration.

=over 4

=item listRuns($confRef)

The procedure removes run from run configuration.

Parameters:
    $confRef       - reference to hash map with all runs configuration 
    $runId         - run id to be removed  
    $runConfigFile - config file name
Returns:
    nothing
    
=back

=cut

#sub removeRun(\%;$;$) {
#    croak "ERROR: removeRun\n" unless ( @_ == 3 );
#    my ( $confRef, $runId, $runConfigFile ) = @_;
#    if ( !exists $runsConfig{run}->{$runId} ) {
#        errorExit "ERROR: removeRun run [$runId] does not exists !$\n";
#    }
#    delete $runsConfig{run}->{$runId};
#    writeRunsConfigXML( %{$confRef}, $runConfigFile );
#}

#=pod
#
#=head1 The procedure reads HASH MAP with runs config from the run XML config file.
#
#=over 4
#
#=item readRunsConfigXML($confRef, $fileName)
#
#The procedure reads HASH MAP with runs config from the run config file.
#
#Parameters:
#    $confRef     - reference to hash map with all runs configuration 
#    $fileName    - config file name
#    $expand      - if 1 then experiments will be expanded with runs
#   
#Returns:
#    nothing
#    
#=back
#
#=cut

#sub readRunsConfigXML(\%;$;$) {
#    croak "ERROR: readRunsConfigXML\n" unless ( @_ == 3 );
#    my ( $confRef, $fileName, $expand ) = @_;
#    my $xs = new XML::Simple(
#        searchpath => ".",
#        KeyAttr    => { run => "+id", set => "+setid", experiment => "+id" },
#        forcearray => 1,
#    );
#    if ( !-e $fileName ) {
#        errorExit "ERROR: readRunsConfigXML file $fileName doesn't exist.\n";
#    }
#    my $ref = $xs->XMLin($fileName);
#    my ( $key, $value );
#    foreach my $key ( keys %{$ref} ) {
#        ${$confRef}{$key} = ${$ref}{$key};
#    }
#    if ( $expand == 1 ) {
#        for my $runId ( sort keys %{ $runsConfig{run} } ) {
#            my %sets = %{ $runsConfig{run}->{$runId}->{set} };
#            foreach my $setId ( keys %sets ) {
#                my @lanes =
#                  @{ $runsConfig{run}->{$runId}->{set}->{$setId}->{lanes} };
#                my $runDir =
#                  $runsConfig{run}->{$runId}->{set}->{$setId}->{gerald};
#                my $expId =
#                  $runsConfig{run}->{$runId}->{set}->{$setId}->{expid};
#                ${$confRef}{experiment}->{$expId}->{run}->{$runId}->{set}
#                  ->{$setId}->{gerald} =
#                  $runsConfig{run}->{$runId}->{set}->{$setId}->{gerald};
#                ${$confRef}{experiment}->{$expId}->{run}->{$runId}->{set}
#                  ->{$setId}->{expid} =
#                  $runsConfig{run}->{$runId}->{set}->{$setId}->{expid};
#                my @tmp =
#                  @{ $runsConfig{run}->{$runId}->{set}->{$setId}->{lanes} };
#                ${$confRef}{experiment}->{$expId}->{run}->{$runId}->{set}
#                  ->{$setId}->{lanes} = \@tmp;
#            }
#        }    # for
#    }
#
#    #   use Data::Dumper;
#    #   print Dumper($ref);
#}

=pod

=head1 The procedure reads HASH MAP with sample directories from the run.conf.xml file.

=over 4

=item readSampleDirsConfigXML($fileName, $sampleDirRef, $exportSetsRef)

The procedure reads HASH MAP with sample directories from the run.conf.xml file.

Parameters:
    $sampleDirRef - reference to hash map with sample directories
    $fileName     - config file name

Returns:
    nothing

=back

=cut

sub readSampleDirsConfigXML($\@;\@) {
    my ( $xmlFilename, $sampleDirRef, $exportSetsRef) = @_;
    
    my $xs = new XML::Simple(
           KeyAttr    => [],
           ForceArray => ['sampleDirectory', 'exportFiles'],
    );
    
    if ( !-e $xmlFilename ) {
        errorExit "ERROR: Could not find the samples configuration file ($xmlFilename)\n";
    }
    
    # read in the sample directories from the XML file
    my $xmlData = $xs->XMLin($xmlFilename);
    
    # convert the relative paths to absolute paths
    foreach my $dir (@{$xmlData->{inputDirectories}->{sampleDirectory}}) {
        my $absDir = abs_path($dir);
        errorExit("ERROR: The following sample directory does not exist: " . $dir . ". Please check your run.conf.xml configuration file.") unless (-e $absDir);
        push @{$sampleDirRef}, $absDir;
    }
    if (defined $exportSetsRef)
    {
        foreach my $exportSet (@{$xmlData->{exportFiles}}) {
            push @{$exportSetsRef}, $exportSet;
        }
    }
    
}

=pod

=head1 The procedure writes HASH MAP with runs config to the run config file.

=over 4

=item writeRunsConfig($confRef, $fileName)

The procedure writes HASH MAP with runs config to the run config file.

Parameters:
    $confRef     - reference to hash map with all runs configuration 
    $fileName    - config file name
   
Returns:
    nothing
    
=back

=cut

sub writeRunsConfig(\%;$) {
    croak "ERROR: writeRunsConfig\n" unless ( @_ == 2 );
    my ( $confRef, $fileName ) = @_;
    my %runDirs = ();

    #    print %{$confRef};
    foreach my $runId ( keys %{$confRef} ) {
        $runDirs{$runId} = ${$confRef}{$runId}->{set1}->{gerald};
    }    # foreach
    my %runLanes = ();
    foreach my $runId ( keys %{$confRef} ) {
        $runLanes{$runId} =
          join( ',', @{ ${$confRef}{$runId}->{set1}->{lanes} } );
    }    # foreach
    writeParameters( %runDirs,  'CONF_RUN_DIRECTORIES', $fileName, 0, 50 );
    writeParameters( %runLanes, 'CONF_RUN_LANES',       $fileName, 0, 50 );
}

=pod

=head1 The procedure reads runs configuration and fixes run id if necessary.

=over 4

=item repairerNonUniqueRunIds($confRef)

The procedure reads runs configuration and fixes run id if necessary.
Checks if all run ids are unique and non empty and replaces run id 
with int(check sum(s_1_export.txt)) if necessary

Parameters:
    $confRef     - reference to hash map with all runs configuration 
   
Returns:
    nothing
    
=back

=cut

#sub repairNonUniqueRunIds(\%) {
#    croak "ERROR: repairNonUniqueRunIds\n" unless ( @_ == 1 );
#
#    my ($confRef) = @_;
#    my %uniqueRunIds = ();
#    foreach my $runIdTmp ( sort keys %{ $confRef->{run} } ) {
#        my $machineName         = "";
#        my $exportRunId         = "";
#        my $fixId               = 'none';
#        my %sets                = %{ $confRef->{run}{$runIdTmp}{set} };
#        my $firstExportFilePath = "";
#        if (   !defined $confRef->{run}{$runIdTmp}{machine}
#            || !defined $confRef->{run}{$runIdTmp}{exportRunId}
#            || !$confRef->{run}{$runIdTmp}{fixId} )
#        {
#            foreach my $setIdTmp ( keys %sets ) {
#
#                #print "$runIdTmp $setIdTmp\n";
#                my $geraldFullDir =
#                  $confRef->{run}{$runIdTmp}{set}{$setIdTmp}{gerald};
#                my @lanes =
#                  sort { $a <=> $b }
#                  @{ $confRef->{run}{$runIdTmp}{set}{$setIdTmp}{lanes} };
#                my $readMode = 'paired';
#                if (
#                    exists $confRef->{run}{$runIdTmp}{set}{$setIdTmp}
#                    {readMode} )
#                {
#                    $readMode =
#                      $confRef->{run}{$runIdTmp}{set}{$setIdTmp}{readMode};
#                }
#                ## check machine name and run id from export files based on 1th export
#                my $laneTmp        = shift @lanes;
#                my $exportFilePath = "";
#                if ( $readMode eq 'paired' ) {
#                    $exportFilePath =
#                      File::Spec->catfile( $geraldFullDir,
#                        "s_$laneTmp\_1_export.txt" );
#                }
#                else {
#                    $exportFilePath =
#                      File::Spec->catfile( $geraldFullDir,
#                        "s_$laneTmp\_export.txt" );
#                }
#                if ( ! -s $exportFilePath and -s "$exportFilePath.gz" ) {
#                    $exportFilePath = "$exportFilePath.gz";
#                }
#                if ( $firstExportFilePath eq "" ) {
#                    $firstExportFilePath = $exportFilePath;
#                }
#
#                my $exportFile = IO::File->new("gunzip -cf $exportFilePath |")
#                  || errorExit "ERROR repairNonUniqueRunIds() Couldn't "
#                  . "open file handle for $exportFilePath $!\n";
#                my $line = <$exportFile>;
#                if ( defined $line && $line ne '' ) {
#                    my @row = split "\t", $line;
#                    $machineName = $row[0];
#                    $exportRunId = $row[1];
#                }    # if
#                close($exportFile);
#                if ( $exportRunId eq "" && $machineName =~ /_/ ) {
#                    my @tmp = split "_", $machineName;
#                    if ( $tmp[1] =~ /^(\d+)$/ ) {
#                        $fixId       = "parsed";
#                        $machineName = $tmp[0];
#                        $exportRunId = $tmp[1];
#                    }
#                    else {
#                    }
#                }
#                if ( defined $confRef->{run}{$runIdTmp}{machine} ) {
#                    if (   $confRef->{run}{$runIdTmp}{machine} ne $machineName
#                        || $confRef->{run}{$runIdTmp}{exportRunId} ne
#                        $exportRunId )
#                    {
#                        errorExit
#"ERROR: run $machineName:$exportRunId not unique for $runIdTmp\n";
#                    }
#                }
#            }    # foreach
#        }    # if
#        else {
#            $machineName = $confRef->{run}{$runIdTmp}{machine};
#            $exportRunId = $confRef->{run}{$runIdTmp}{exportRunId};
#            $fixId       = $confRef->{run}{$runIdTmp}{fixId};
#        }
#        ## check for unique run ids
#        my $uniqueRunId = $machineName . "_" . $exportRunId;
#        if ( $exportRunId eq "" || defined $uniqueRunIds{$uniqueRunId} ) {
#            $fixId = 'generated';
#            my $exportFile = IO::File->new( "gunzip -cf $firstExportFilePath |" )
#              || errorExit "ERROR repairNonUniqueRunIds() Couldn't "
#              . "open file handle for $firstExportFilePath $!\n";
#            my $content  = "";
#            my $headSize = 100;
#            for ( my $index = 0 ; $index < $headSize ; $index++ ) {
#                my $line = <$exportFile>;
#                $content .= $line;
#            }    # if
#            close($exportFile);
#            my $digest = md5_hex($content);
#            $exportRunId = hex( substr( $digest, 0, 8 ) );
#        }
#        else {
#            $uniqueRunIds{$uniqueRunId} = $confRef->{run}{$runIdTmp};
#        }
#        $confRef->{run}{$runIdTmp}{machine}     = $machineName;
#        $confRef->{run}{$runIdTmp}{exportRunId} = $exportRunId;
#        $confRef->{run}{$runIdTmp}{fixId}       = $fixId;
#    }    # foreach
#}    # repairerNonUniqueRunIds

=pod

=head1 The procedure writes HASH MAP with runs config to the xml run config file.

=over 4

=item writeRunsConfigXML($confRef, $fileName)

The procedure writes HASH MAP with runs config to the run xml config file.

Parameters:
    $confRef     - reference to hash map with all runs configuration 
    $fileName    - config file name
   
Returns:
    nothing
    
=back

=cut

#sub writeRunsConfigXML(\%;$) {
#    croak "ERROR: writeRunsConfig\n" unless ( @_ == 2 );
#    my ( $confRef, $fileName ) = @_;
#    unless (!keys %$confRef or !keys %{$confRef->{run}})
#    {
#        my $xs = new XML::Simple(
#            searchpath => ".",
#            KeyAttr    => { run => "+id", set => "+setid", experiment => "+id", file => "+id" },
#            forcearray => 1,
#        );
#        my $xml  = $xs->XMLout($confRef);
#        my $file = IO::File->new( ">" . $fileName )
#          || errorExit
#    "ERROR: writeRunsConfig() Couldn't create/open file handle for $fileName $!\n";
#        print ($file $xml);
#        close $file;
#    }
#}

=pod

=head1 The procedure writes HASH MAP with sample directories to the run.conf.xml file.

=over 4

=item writeSampleDirsConfigXML($sampleDirRef, $fileName)

The procedure writes HASH MAP with sample directories to the run.conf.xml file.

Parameters:
    $sampleDirRef - reference to hash map with input sample directories
    $fileName     - config file name

Returns:
    nothing

=back

=cut

sub writeSampleDirsConfigXML(\%;$) {
    croak "ERROR: writeSampleDirsConfigXML\n" unless ( @_ == 2 );
    my ( $sampleDirRef, $fileName ) = @_;
    
    unless (!keys %$sampleDirRef) {

        my $xs = new XML::Simple(
           XMLDecl    => 1,
           KeyAttr    => [],
           ForceArray => 1,
#           KeepRoot   => 1
        );
        
        my $xml  = $xs->XMLout($sampleDirRef);
        my $file = IO::File->new( ">" . $fileName ) || errorExit "ERROR: writeSampleDirsConfigXML() couldn't open $fileName for writing. $!\n";
        print ($file $xml);
        close $file;
    }
}

=pod

=head1 The procedure adds run config to hash map with run configs.

=over 4

=item updateRunConf($confNewRunRef, $newRunId, $confRunsRef, $confProjRef)

The procedure adds run config to hash map with run configs.

Parameters:
    $confNewRunRef - new run config
    $newRunId      - new run config id
    $confRunsRef   - all runs configuration
    $confProjRef   - project configuration
   
Returns:
    nothing
    
=back

=cut

sub updateRunConf(\%;$;\%;\%) {
    croak "ERROR: updateRunConf [@_]\n" unless ( @_ == 4 );
    my ( $confNewRunRef, $newRunId, $confRunsRef, $confProjRef ) = @_;
    ${$confNewRunRef}{_runId}     = $newRunId;
    ${$confNewRunRef}{appVersion} = 'bcl2fastq-1.8.4';
    ${$confRunsRef}{$newRunId}    = $confNewRunRef;
    ${$confProjRef}{currentRunId} = $newRunId;
}

=pod

=head1 The procedure add/update current run configuration value.

=over 4

=item updateRunConfValue($key, $value)

The procedure add/update current run configuration  value.

Parameters:
    $key   - key name in run configuration
    $value - value name in run configuration
   
Returns:
    nothing
    
=back

=cut

sub updateRunConfValue($;$) {
    croak " ERROR: updateRunConfValue [@_] \n " unless ( @_ == 2 );
    my ( $key, $value ) = @_;
    $CONF_RUN{$key} = $value;
    my $projectConf = File::Spec->catfile( $CONF_PROJ{dirBuild}, $CONF_APP{dirConf}, $CONF_APP{projectConf} );
    writeParameters( %CONF_RUN, 'RUN_' . $CONF_PROJ{currentRunId}, $projectConf, 0 );
}

=pod

=head1 The procedure add/update current project configuration value.

=over 4

=item updateProjectConfValue($key, $value)

The procedure add/update current project configuration  value.

Parameters:
    $key   - key name in run configuration
    $value - value name in run configuration
   
Returns:
    nothing
    
=back

=cut

sub updateProjectConfValue($;$) {
    croak " ERROR: updateProjectConfValue [@_] \n " unless ( @_ == 2 );
    my ( $key, $value ) = @_;
    $CONF_PROJ{$key} = $value;
    my $projectConf = File::Spec->catfile( $CONF_PROJ{dirBuild}, $CONF_APP{dirConf}, $CONF_APP{projectConf} );
    writeParameters( %CONF_PROJ, 'PROJECT_CONFIGURATION', $projectConf, 0 );
}

=pod

=head1 The procedure gets current run configuration.

=over 4

=item updateProjectConfValue($runId, $confRunsRef, $confProjRef)

The procedure gets current run configuration.

Parameters:
    $runId       -
    $confRunsRef -
    $confProjRef -
Returns:
    RUN configu value
    
=back

=cut

sub getRunConf($;\%;\%) {
    croak "ERROR: getRunConf \n" unless ( @_ == 3 );
    my ( $runId, $confRunsRef, $confProjRef ) = @_;
    return ${$confRunsRef}{$runId};
}

=pod

=head1 The procedure configures/creates directories.

=over 4

=item configureDirectories( $dirBuild, $dirBuildParsed, $chromsSizesRef)

The procedure configures/creates directories based on Chromosomes names and
bin sizes

Parameters:
    $dirBuild       - project directory
    $dirBuildParsed - build directory
    $chromsSizesRef - mapRef with chromosomes sizes
Returns:
    nothing
    
=back

=cut

# TODO refactor: this looks ...
sub configureDirectories($$\%) {
    my ( $dirBuild, $dirBuildParsed, $chromsSizesRef) = @_;
    printLog( sprintf("DEBUG: configureDirectories\n"), 10 );
    my @dirs = ();
    my $dirStructRef = \@dirs;
    
    # initialize
    my %numberOfBinsInChromsProject = ();
    my %numberOfBinsInChromsBuild   = ();
    my @buildDirStructure           = ();
    my $projectBinSize              = $CONF_PROJ{binSizeProject};
    my $buildBinSize                = $CONF_PROJ{binSizeBuild};
    my $dirBuildExport              = $CONF_PROJ{dirBuildExport};
    
    # define our directories
    my $dirBuildExportSets = File::Spec->catdir( $dirBuildExport,     'sets' );
    my $templateDir        = File::Spec->catdir( $dirBuildExportSets, 'template' );
    my $templateBuildDir   = File::Spec->catdir( $dirBuildExportSets, 'template.build' );
    my $html               = File::Spec->catdir( $dirBuild,           'html' );
    my $stats              = File::Spec->catdir( $dirBuild,           'stats' );
    my $htmlSets           = File::Spec->catdir( $html,               'sets' );
    my $statsSets          = File::Spec->catdir( $stats,              'sets' );
    
    # add these directories to our directory creation array
    push @{$dirStructRef}, $CONF_PROJ{dirBuild};
    push @{$dirStructRef}, $dirBuildExport;
    push @{$dirStructRef}, $dirBuildParsed;
    push @{$dirStructRef}, $html;
    push @{$dirStructRef}, $stats;
    push @{$dirStructRef}, $dirBuildExportSets;
    push @{$dirStructRef}, $templateDir;
    push @{$dirStructRef}, $templateBuildDir;
    push @{$dirStructRef}, $htmlSets;
    push @{$dirStructRef}, $statsSets;
    
    # add reference sequence-specific directories to our directory creation array
    configureChromDirectories( $templateDir, $projectBinSize,
        %{$chromsSizesRef}, %numberOfBinsInChromsProject, @{$dirStructRef} );
    configureChromDirectories( $dirBuildParsed, $buildBinSize,
        %{$chromsSizesRef}, %numberOfBinsInChromsBuild, @{$dirStructRef} );
    configureChromDirectories( $templateBuildDir, $buildBinSize,
        %{$chromsSizesRef}, %numberOfBinsInChromsBuild, @{$dirStructRef} );

    # add the NMNM directories to our directory creation array if we're keeping all reads
    my $isKeepAllReads      = ((defined $CONF_PROJ{sortKeepAllReads}) and $CONF_PROJ{sortKeepAllReads});
    my $isKeepUnmappedReads = $isKeepAllReads;
    if($isKeepUnmappedReads) {
        push @{$dirStructRef},
          File::Spec->catdir( $templateDir, $CONF_APP{TAG_NMNM} );
        push @{$dirStructRef},
          File::Spec->catdir( $templateBuildDir, $CONF_APP{TAG_NMNM} );
    }

    # add the RNA features directory to our directory creation array
    if ( defined $CONF_PROJ{applicationType}
        && $CONF_PROJ{applicationType} eq 'RNA' )
    {
        push @{$dirStructRef}, File::Spec->catdir( $dirBuild, 'features' );
    }
    
    
    if ( $CONF_PROJ{readMode} eq 'single' )
    {    # at thsi point only one mode can be supported
=cut	

		my $singleReadDir =
		  File::Spec->catdir( $dirBuildExport, 'single-read' );
		push @{$dirStructRef}, $singleReadDir;
		my $singleReadNMNMDir =
		  File::Spec->catdir( $singleReadDir, $CONF_APP{TAG_NMNM} );
		push @{$dirStructRef}, $singleReadNMNMDir;
		configureChromDirectories( $singleReadDir, $buildBinSize,
			%{$chromsSizesRef}, %numberOfBinsInChromsBuild, @{$dirStructRef} );
		## Remove when in hybrid mode
		%numberOfBinsInChromsProject = ();
		my @tmp = ();
		configureChromDirectories( $dirBuildExport, $projectBinSize,
			%{$chromsSizesRef}, %numberOfBinsInChromsProject, @tmp );
		if ( defined $CONF_APP{applicationType}
			&& $CONF_APP{applicationType} eq 'RNA' )
		{
			my $singleReadRMDir =
			  File::Spec->catdir( $singleReadDir, $CONF_APP{TAG_RM} );
			push @{$dirStructRef}, $singleReadNMNMDir;
		}
		else {
			errorExit "ERROR RNA\n";
		}
=cut

    } else {
        # $CONF_PROJ{readMode} != 'single'
        %numberOfBinsInChromsProject = ();
        configureChromDirectories( $dirBuildExport, $projectBinSize,
            %{$chromsSizesRef}, %numberOfBinsInChromsProject,
            @{$dirStructRef} );
    }
    
    # add another NMNM directory to our directory creation array
    if($isKeepUnmappedReads) {
        my $NMNM_DIR = File::Spec->catdir( $dirBuildParsed, $CONF_APP{TAG_NMNM} );
        push @{$dirStructRef}, $NMNM_DIR;
    }

    my $projectConf = File::Spec->catfile( $CONF_PROJ{dirBuild}, $CONF_APP{dirConf}, $CONF_APP{projectConf} );
    writeParameters( %numberOfBinsInChromsProject, "PROJECT_BIN_SIZES", $projectConf, 0 );
    writeParameters( %numberOfBinsInChromsBuild, "BUILD_BIN_SIZES", $projectConf, 0 );
    
    # add the temporary split directories
    # N.B. this creates one temporary directory for each export set.
    for my $i (1 .. scalar(@{$runsConfig{exportFiles}})) {
        push @{$dirStructRef}, File::Spec->catdir($dirBuildExportSets, $i);
        push @{$dirStructRef}, File::Spec->catdir($htmlSets,           $i);
        push @{$dirStructRef}, File::Spec->catdir($statsSets,          $i);
    }

    createDirs( "", @{$dirStructRef} );
    my $confDir = File::Spec->catdir( $dirBuild , $CONF_APP{dirConf} );
    my $dirsPath = File::Spec->catfile( $confDir , "dirs.tar" );
    my $cmd = " tar -C $templateBuildDir -cf $dirsPath ./";
    executeCmd( $cmd, 5 );
    my $pdirsPath = File::Spec->catfile( $confDir , "project.dirs.tar" );
    $cmd = "tar -C $templateDir -cf $pdirsPath ./";
    executeCmd( $cmd, 5 );
    executeCmd( "rm -fr $templateDir", 5 );
    executeCmd( "rm -fr $templateBuildDir", 5 );
}

=pod

=head1 The procedure configures/creates chromosome directories .

=over 4

=item configureChromDirectories( $location, $binSize, $chromName, $chromsSizesRef,
        $numberOfBinsInChromsRef, $dirStructureRef)

The procedure configures/creates directories based on Chromosomes names and
bin sizes

Parameters:
    $location - path where the directory shuld be created
    $binSize - size of bins           
    $chromsSizesRef - mapRef with chromosomes sizes
    $numberOfBinsInChromsRef - mapRef where the number of bins in each chromosome should be stored
    $dirStructureRef - listRef where the structure of directories should be stored

Returns:
    nothing
    
=back

=cut

sub configureChromDirectories($$\%\%\@) {
    my ( $location, $binSize, $chromsSizesRef, $numberOfBinsInChromsRef,
        $dirStructureRef )
      = @_;
    my $i;
    
    foreach my $chromKey ( keys %{$chromsSizesRef} ) {
        my $chromDir     = "$location/$chromKey/";
        my $numberOfBins =
          ceil( ${$chromsSizesRef}{$chromKey}->{length} / $binSize );
        push @{$dirStructureRef}, $chromDir;
        ${$numberOfBinsInChromsRef}{$chromKey} = $numberOfBins;
        for ( $i = 0 ; $i < $numberOfBins ; $i++ ) {
            my $binName = sprintf "%04d", $i;
            my $binDir = "$chromDir$binName/";
            push @{$dirStructureRef}, $binDir;
        }    # for
    }    # foreach
    return;
}

=pod

=item createExperiments($runsConfigRef)

The createExperiments converts export files to fastq

B<Parameters:>

    $buildSampleDir     - path to project directory
    $runsConfigRef  - HASH MAP Ref to runs.conf.xml
    
B<Returns:> 
    nothing
    
=cut

#sub createExperiments ($;$;\%) {
#    croak "ERROR: createExperiments " . join( ',', @_ ) unless ( @_ == 3 );
#    my ( $buildSampleDir, $runConfigFile, $runsConfigRef ) = @_;
#
#    for my $runIdTmp ( sort keys %{ $runsConfigRef->{run} } ) {
#        my %sets = %{ $runsConfigRef->{run}{$runIdTmp}{set} };
#        foreach my $setIdTmp ( keys %sets ) {
#            my @lanes =
#              @{ $runsConfigRef->{run}{$runIdTmp}{set}{$setIdTmp}{lanes} };
#            foreach my $laneTmp (@lanes) {
#                my @reads      = ();
#                my $reportPath =
#                  File::Spec->catdir( $buildSampleDir, "stats", 'runs', $runIdTmp,
#                    $setIdTmp );
#                my $geraldFullDir =
#                  $runsConfigRef->{run}->{$runIdTmp}->{set}->{$setIdTmp}
#                  ->{gerald};
#                my $localLanePath =
#                  File::Spec->catdir( $buildSampleDir, "export", 'lanes', $runIdTmp,
#                    $laneTmp );
#
#                my $readMode =
#                  $runsConfigRef->{run}->{$runIdTmp}->{set}->{$setIdTmp}
#                  ->{readMode};
#                my %sraConf = ();
#                my $expId   = "unknown";
#
#                if ( $readMode eq 'single' ) {
#                    @reads = (1);
#                }
#                else {
#                    @reads = ( 1, 2 );
#                    my $pairsPath =
#                      File::Spec->catfile( $reportPath,
#                        "s_$laneTmp\_pair.xml" );
#                    my %pairs;
#                    if (readXML( %pairs, $pairsPath, 0) )
#                    {
#                        if (!defined $pairs{InsertSize}[0]{Median}[0])
#                        {
#                            errorExit "ERROR: No median insert size defined in $pairsPath";
#                        }
#                        else
#                        {
#                            my $insertSize = int( $pairs{InsertSize}[0]{Median}[0] / 100 ) * 100;
#                            my $readLength1 = $pairs{Length}[0]{Read1}[0]{Total}[0];
#                            my $readLength2 = $pairs{Length}[0]{Read2}[0]{Total}[0];
#                            $expId = sprintf( "%dx%dx%d", $insertSize, $readLength1, $readLength2 );
#                            #$runsConfigRef->{experiment}->{$expId}->{id} = $expId;    
#                            $runsConfigRef->{experiment}->{$expId}->{id} = $expId;
#                            $runsConfigRef->{experiment}->{$expId}->{readLength1} = $readLength1;
#                            $runsConfigRef->{experiment}->{$expId}->{readLength2} = $readLength2;
#                            $runsConfigRef->{experiment}->{$expId}->{insertSize} = $insertSize;
#                            $runsConfigRef->{experiment}->{$expId}{readType} = 'double';
#                            $runsConfigRef->{run}->{$runIdTmp}->{set}->{$setIdTmp}->{expid} = $expId;
#                        }
#                    }
#                }
#            };    # foreach
#        };    # foreach
#    };    # foreach
#    writeRunsConfigXML( %{$runsConfigRef}, $runConfigFile );
#
#}     # createExperiments

=pod

=item isSpliceJunctionChrom($chrom, $CONF_PROJ_Ref)

Checks if chromosome matches the splice junction file name

B<Parameters:>

    $chrom           - Chromosome name
    $CONF_PROJ_Ref   - HASH MAP ref to project configuration
    
B<Returns:> 
    1 if chromosome matches the splice junction file name
    
=cut

sub isSpliceJunctionChrom ($\%) {
    my ( $chrom, $CONF_PROJ_Ref ) = @_;
    
    return ($chrom eq getSpliceJunctionFileName(%$CONF_PROJ_Ref));
}

=pod

=item getSpliceJunctionFileName($CONF_PROJ_Ref)

Returns the name of configured splice junction file

B<Parameters:>

    $CONF_PROJ_Ref   - HASH MAP ref to project configuration
    
B<Returns:> 
    file name
    
=cut

sub getSpliceJunctionFileName (\%) {
    my ( $CONF_PROJ_Ref ) = @_;

    return (File::Spec->splitpath ($CONF_PROJ_Ref->{spliceJunction}))[2];
}

1;    # says use was ok
__END__

