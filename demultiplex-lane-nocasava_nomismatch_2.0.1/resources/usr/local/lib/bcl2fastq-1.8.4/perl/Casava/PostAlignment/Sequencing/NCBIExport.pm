package Casava::PostAlignment::Sequencing::NCBIExport;

# PROJECT: NCBISubmission
# MODULE:  $RCSfile: NCBIExport.pm,v $
# AUTHOR:  Come Raczy, Lukasz Szajkowski
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#
# The library contains procedures and variables usefull in operating on Illumina
# genomic data e.g.: export to SRF and generating SRA xml files

=pod

=head1 NAME

Casava::PostAlignment::Sequencing::NCBIExport.pm - The library - NCBI exported library
	

=head1 SYNOPSIS

The library contains procedures and variables usefull in operating on Illumina
# genomic data e.g.: export to SRF and generating SRA xml files.
 
use Casava::PostAlignment::Sequencing::NCBIExport.pm qw();  

=head1 DESCRIPTION

Exports: 
	generateSrfXML($;\%)
		    
# Global variable


=head1 AUTHORSHIP

Copyright (c) 2008, 2009 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).


=cut

BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT =
      qw(&generateSrfXML &export2fastq &fastq2zip &createSraExperiments);
    @EXPORT_OK = qw();
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);
use Carp;

use File::Spec;

use Casava::Common::Log;
use Casava::Common::IOLib qw(executeCmd
  readXML writeXML string2file);
use Casava::Common::Template;
use Casava::PostAlignment::Sequencing::Config qw(%CONF_APP);
sub generateSrfXML($;\%;\%);
sub export2fastq ($;\%;$;$;$);
sub fastq2zip( $;\%;$;$;$ );
sub createSraExperiments ($;\%;\%);
use Cwd;
use Cwd 'abs_path';

=pod

=item export2fastq($msg)

The export2fastq converts export files to fastq

B<Parameters:>

    $projectDir     - path to project directory
    $runsConfigRef  - HASH MAP Ref to runs.conf.xml
    $runId          - run id
    $lane           - lane id
    $setId          - set id
    
B<Returns:> 
    nothing
    
=cut

sub export2fastq ($;\%;$;$;$) {
    croak "ERROR: export2fastq " . join( ',', @_ ) unless ( @_ == 5 );
    my ( $projectDir, $runsConfigRef, $runId, $lane, $setId ) = @_;

    my $export2FastqCmd = "$CONF_APP{dirLibExec}/CASAVA -a Export2Fastq";

    my @reads;
    my $cmd        = "";
    my $reportPath =
      File::Spec->catdir( $projectDir, "stats", 'runs', $runId, $setId );
    my $geraldFullDir =
      $runsConfigRef->{run}->{$runId}->{set}->{$setId}->{gerald};
    my $localLanePath =
      File::Spec->catdir( $projectDir, "export", 'lanes', $runId, $lane );

    my $readMode = $runsConfigRef->{run}->{$runId}->{set}->{$setId}->{readMode};

    if ( $readMode eq 'single' ) {
        @reads = (1);
    }
    else {
        @reads = ( 1, 2 );
    }

    ## confert to fastq
    foreach my $readId (@reads) {
        my $readStr = "";
        if ( scalar(@reads) > 1 ) {
            $readStr = "_$readId";
        }
        my $exportPath =
          File::Spec->catfile( $geraldFullDir, "s_$lane$readStr\_export.txt" );
        if ( ! -s $exportPath and -s "$exportPath.gz" ) {
            $exportPath = "$exportPath.gz";
        }
        my $e2fReportPath =
          File::Spec->catfile( $reportPath, "s_$lane$readStr\_e2f.xml" );

        my $fastq1Path =
          File::Spec->catfile( $localLanePath, "s_$lane$readStr.fastq" );
        if ( -e $fastq1Path ) {
            $cmd = "rm -f $e2fReportPath";
            executeCmd( $cmd, 0 );
            $cmd = "rm -f $fastq1Path";
            executeCmd( $cmd, 0 );
        }
        $cmd =
            "$export2FastqCmd -e $exportPath "
          . "-o $fastq1Path -s $e2fReportPath --purityFilter=YES";
        executeCmd( $cmd, 0 );
    }    # foreach

}    # export2fastq

=pod

=item export2fastq($msg)

The export2fastq converts export files to fastq

B<Parameters:>

    $projectDir     - path to project directory
    $runsConfigRef  - HASH MAP Ref to runs.conf.xml
    $runId          - run id
    $lane           - lane id
    $setId          - set id
    
B<Returns:> 
    nothing
    
=cut

sub fastq2zip ($;\%;$;$;$) {
    croak "ERROR: export2fastq " . join( ',', @_ ) unless ( @_ == 5 );
    my ( $projectDir, $runsConfigRef, $runId, $lane, $setId ) = @_;
    my @reads;
    my $cmd = "";

    my $submissionPath = File::Spec->catdir( $projectDir, "submission" );
    unless ( -d $submissionPath ) {
        mkdir $submissionPath, 0775
          or errorExit "ERROR: $0 can't make $submissionPath $!\n";
    }    # unless
    my $reportPath =
      File::Spec->catdir( $projectDir, "stats", 'runs', $runId, $setId );
    my $geraldFullDir = $runsConfigRef->{run}{$runId}{set}{$setId}{gerald};
    my $expId         = $runsConfigRef->{run}{$runId}{set}{$setId}{expid};
    my $localLanePath =
      File::Spec->catdir( $projectDir, "export", 'lanes', $runId, $lane );
    my $readMode = $runsConfigRef->{run}{$runId}{set}{$setId}{readMode};

    my %sraXml = ();
    if ( $readMode eq 'single' ) {
        @reads = (1);
        my $e2fReportPath =
          File::Spec->catfile( $reportPath, "s_$lane\_e2f.xml" );
    }
    else {
        @reads = ( 1, 2 );
        my $e2fReportPath1 =
          File::Spec->catfile( $reportPath, "s_$lane\_1_e2f.xml" );
        my $e2fReportPath2 =
          File::Spec->catfile( $reportPath, "s_$lane\_2_e2f.xml" );
        my %conf1;
        my %conf2;
        readXML( %conf1, $e2fReportPath1 );
        readXML( %conf2, $e2fReportPath2 );
        if ( $conf1{minReadLen} != $conf1{maxReadLen} ) {
            errorExit "ERROR: Uncosistent read length in $e2fReportPath1 \n";
        }
        if ( $conf2{minReadLen} != $conf2{maxReadLen} ) {
            errorExit "ERROR: Uncosistent read length in $e2fReportPath2 \n";
        }
        if ( $conf1{totalReads} != $conf2{totalReads} ) {
            errorExit "ERROR: Uncosistent number of reads "
              . $conf1{totalReads} . "/"
              . $conf2{totalReads}
              . " in $e2fReportPath1 \n";
        }
        $sraXml{file}{totalReads} = $conf1{totalReads} + $conf2{totalReads};
        $sraXml{file}{totalBases} = $conf1{totalBases} + $conf2{totalBases};
    }

    ## confert to fastq
    my $filesStr = '';
    foreach my $readId (@reads) {
        my $readStr = "";
        if ( scalar(@reads) > 1 ) {
            $readStr = "_$readId";
        }
        my $fastq1Path =
          File::Spec->catfile( $localLanePath, "s_$lane$readStr.fastq" );
        if ( -e $fastq1Path ) {
            $filesStr .= "$fastq1Path ";
        }
    }    # foreach
    my $tarFileName = "$expId-$runId-$lane.tar.gz";
    my $tarFilePath = File::Spec->catfile( $submissionPath, $tarFileName );
    $cmd = "tar -C $localLanePath -czf $tarFilePath $filesStr";

    executeCmd( $cmd, 0 );

    my $tarReportFile =
      File::Spec->catfile( $reportPath, "$expId-$runId-$lane\_verify.txt" );
    $cmd = "tar -ztf $tarFilePath 2> $tarReportFile";
    executeCmd( $cmd, 0 );

    my $command = "/usr/bin/sha1sum -b $tarFilePath | cut -d ' ' -f 1";
    my $sha1sum = "";
    open( FILE, "$command |" )
      || errorExit "ERROR: fastq2zip() Couldn't run $command \n";
    while (<FILE>) {
        chomp($_);
        $sha1sum = $_;
    }    # while
    close(FILE);
    if ( $sha1sum eq "" ) {
        errorExit
          "ERROR: fastq2zip() Couldn't run compute sha1sum for $tarFilePath \n";
    }
    $sraXml{file}{expid}           = $expId;
    $sraXml{file}{runid}           = $runId;
    $sraXml{file}{setid}           = $setId;
    $sraXml{file}{lane}            = $lane;
    $sraXml{file}{id}              = $tarFileName;
    $sraXml{file}{filename}        = $tarFileName;
    $sraXml{file}{filetype}        = "tar.gz";
    $sraXml{file}{method} = "SHA1";
    $sraXml{file}{checksum}        = $sha1sum;
    my $sraReportPath = File::Spec->catfile( $reportPath, "s_$lane\_sra.xml" );
    writeXML( %sraXml, $sraReportPath );

}    # fastq2zip

=pod

=item export2fastq($msg)

The createSraExperiments creates list of experiments with their data. 
The procedure reads runs.conf.xml and sra.xml.

B<Parameters:>

    $projectDir            - path to project directory
    $runsConfigRef         - HASH MAP Ref to runs.conf.xml
    $expConfigRef          - HASH MAP Ref to experiments.conf.xml
    
B<Returns:> 
    nothing
    
=cut

sub createSraExperiments ($;\%;\%) {
    croak "ERROR: createSraExperiments " . join( ',', @_ ) unless ( @_ == 3 );
    my ( $projectDir, $runsConfigRef, $expConfigRef ) = @_;

    for my $runIdTmp ( sort keys %{ $runsConfigRef->{run} } ) {
        my %sets = %{ $runsConfigRef->{run}{$runIdTmp}{set} };
        foreach my $setIdTmp ( keys %sets ) {
            my @lanes =
              @{ $runsConfigRef->{run}{$runIdTmp}{set}{$setIdTmp}{lanes} };
            my $expId = $runsConfigRef->{run}{$runIdTmp}{set}{$setIdTmp}{expid};
            $expConfigRef->{experiment}{$expId}{id} =
              $runsConfigRef->{experiment}{$expId};
            foreach my $laneTmp (@lanes) {
                my $reportPath =
                  File::Spec->catdir( $projectDir, "stats", 'runs', $runIdTmp,
                    $setIdTmp );
                my $sraReportPath =
                  File::Spec->catfile( $reportPath, "s_$laneTmp\_sra.xml" );
                my %sraXml = ();
                readXML( %sraXml, $sraReportPath );
                push @{ $expConfigRef->{experiment}{$expId}{id}{file} },
                  $sraXml{file};
            };    # foreach
        };    # foreach
    };    # foreach

    #    use Data::Dumper;
    #    print Dumper($expConfigRef);

}

=pod

=head1 The procedure generates the XML forms for the submission.

=over 4

=item configureReadStorage($configRef)

The procedure generates the XML forms for the submission.
The metadata is read from the <metadata-directory> and the
resulting XML forms are created in the local directory.

Parameters:
    $metadataDirectory  - path to NCBI directory in the CASAVA Build directory
	$parametersRef      - HASH MAP REF with fields (submissionId, 
    						accessionId, sampleId, contactPerson, studyName)
Returns:
    Nothing
        
=back

=cut

sub generateSrfXML($;\%;\%) {
    croak "ERROR: generateSrfXML " . join( "\n", @_ ) . " " . scalar(@_) . "\n"
      unless ( @_ == 3 );
    my ( $projectDir, $parametersRef, $experimentListRef ) = @_;
    my $submissionId = $parametersRef->{submissionId};

    croak
"The submission Id should be an alphanumeric word, possibly including '-' and '_': $submissionId"
      unless $submissionId =~ /^([-_0-9a-zA-Z]+)$/;

    my $submissionDir = $projectDir . "/submission";
    my $runXml = File::Spec->catfile( $submissionDir, "$submissionId.run.xml" );
    my $submissionXml =
      File::Spec->catfile( $submissionDir, "$submissionId.submission.xml" );
    my $experimentXml =
      File::Spec->catfile( $submissionDir, "$submissionId.experiment.xml" );

    foreach my $file ( $runXml, $submissionXml, $experimentXml ) {
        system("rm -f $file")              if -e $file;
        croak "The file $_ already exists" if -e $file;
    }

    my %runData        = ();
    my %experimentData = ();
    my %submissionData = ();

    foreach my $key ( keys %{$parametersRef} ) {
        $runData{$key}        = $parametersRef->{$key};
        $submissionData{$key} = $parametersRef->{$key};
        $experimentData{$key} = $parametersRef->{$key};
    }    # foreach

    my $templPath = File::Spec->catdir( $CONF_APP{dirLibExec}, "templates" );
    Casava::Common::Template::configure($templPath);
    my $runBlock        = '';
    my $experimentBlock = '';
    my $submissionBlock = '';
    my $checksumtxt = '';
    my $checkSumMethod = '';
    for my $expIdTmp ( sort keys %{ $experimentListRef->{experiment} } ) {

        #print Dumper( $experimentListRef->{experiment}{$expIdTmp}{id} );
        # print Dumper($experimentListRef->{experiment}{$expIdTmp}{id});

        my $expRef = $experimentListRef->{experiment}{$expIdTmp}{id};
        my @files  = @{ $expRef->{file} };
        $experimentData{expid}      = $expIdTmp;
        $experimentData{baseCoord}  = $expRef->{readLength1} + 1;
        $experimentData{readLength} = $expRef->{readLength1};
        $experimentData{insertSize} = $expRef->{insertSize};
        $experimentBlock .=
          fromTemplate( %experimentData, "experimentBlock.xml" );

        foreach my $expRef (@files) {
            foreach my $fileRef ( values %$expRef ) {
                my $runId   = $fileRef->{runid};
                my $runDate = '';
                if ( $runId =~ /^(\d{2)(\d{2)(\d{2)_/ ) {
                    $runDate =
                      sprintf( "20%02d-%02d-%02dT00:00:00+01:00", $1, $2, $3 );
                }
                else {
                    my (
                        $sec,  $min,   $hour,    # second, minute, hour
                        $mday, $month, $year,    # day of month, month
                        $wday, $yday,  $dst
                      )
                      = localtime( time() );
                    $runDate = sprintf(
                        "%02d-%02d-%02dT00:00:00+01:00",
                        1900 + $year,
                        $month, $mday
                    );
                }
                $fileRef->{runDate}    = $runDate;
                $fileRef->{readLength} =
                  $experimentListRef->{experiment}{$expIdTmp}{id}{readLength1};
                $fileRef->{totalSpots} = $fileRef->{totalReads} / 2;
                foreach my $key ( keys %{$parametersRef} ) {
                    $fileRef->{$key} = $parametersRef->{$key};
                }
                $runBlock .= fromTemplate( %{$fileRef}, "runBlock.xml" );
                $submissionBlock .= fromTemplate( %{$fileRef}, "submissionBlock.xml" );
                $checksumtxt .= sprintf("%s  %s\n", $fileRef->{checksum}, $fileRef->{filename});
                $checkSumMethod = $fileRef->{method};
            }
        };    # foreach
    };    # foreach
    $runData{runBlock}               = $runBlock;
    $experimentData{experimentBlock} = $experimentBlock;    
    $submissionData{submissionBlock} = $submissionBlock;
    
    my $runXmlContent = fromTemplate( %runData, "run.xml" );
    string2file( $runXmlContent, $runXml );

    my $experimentXmlContent =
      fromTemplate( %experimentData, "experiment.xml" );
    string2file( $experimentXmlContent, $experimentXml );

    my $submissionXmlContent =
      fromTemplate( %submissionData, "submission.xml" );
    string2file( $submissionXmlContent, $submissionXml );

    my $checksumTxt =
      File::Spec->catfile( $submissionDir, "$submissionId\_checksum_$checkSumMethod.txt" );
    string2file( $checksumtxt, $checksumTxt );
    
}
1;                                                            # says use was ok
__END__

