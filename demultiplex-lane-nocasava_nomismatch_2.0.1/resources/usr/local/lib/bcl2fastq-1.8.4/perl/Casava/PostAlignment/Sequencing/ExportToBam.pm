#!/usr/bin/env perl
# PROJECT: CASAVA
# MODULE:  $RCSfile: ExportToBam.pm,v $
# AUTHOR:  Tony Cox
#
# Copyright (c) 2008 Illumina
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).
#
# this script creates SAM and BAM files out of a sorted or export file
#

package Casava::PostAlignment::Sequencing::ExportToBam;

use strict;
use warnings FATAL => 'all';
use Exporter 'import';

our @EXPORT_OK = qw(sortedToBam);

use File::Spec;
use File::Temp;

use Casava::Common::Log;
use Casava::Common::IOLib qw(executeCmd);
use Casava::PostAlignment::Sequencing::Config
  qw(isSpliceJunctionChrom %chrEnds %CONF_APP %CONF_PROJ);
use Casava::PostAlignment::Sequencing::SamLib qw(getSamtoolsBin getSamHeaderFH);



# convert a commmand which produces single-bin sorted.txt content to
# a single bin bam file
#
sub sortedToBam($\@$$$) {
    my ( $projectPath, $chromsRef, $chrom, $binId, $sortInfo ) = @_;

    # Configuration:
    my $libexecDir      = '/usr/local/libexec/bcl2fastq-1.8.4';
    my $bamDirName      = $CONF_APP{dirBam};
    my $dirBuildParsed  = $CONF_PROJ{dirBuildParsed};
    my $isPaired        = ($CONF_PROJ{readMode} eq 'paired');
    my $isSolexaQvals   = ($CONF_PROJ{qualityType} eq "Solexa64");
    my $isSkipVariableMetadata = (defined $CONF_PROJ{skipVariableMetadata} and $CONF_PROJ{skipVariableMetadata});


    my $sorted2Bam      = File::Spec->catfile($libexecDir, 'sortedToBam');

    # write to incomplete
    my $tmpChrBamFilename = 'sorted.bam.incomplete';

    my $chromPath = File::Spec->catdir($dirBuildParsed,$chrom);
    errorExit "ERROR: Can't find chrom directory: $chromPath\n" unless (-d $chromPath);

    # determine if this bin exists and contains a sorted.txt file
    #
    # warn option (fixed for the time being):
    my $isWarnMissingSorted = 1;

    my $binDir = File::Spec->catdir($chromPath, $binId);
    if (not -d $binDir) {
        logWarning("No bin dir $binDir for chr $chrom - skipping") if ($isWarnMissingSorted);
        return;
    }

    # create reference length file required for sam->bam conversion from CASAVA configuration data:
    my $samHeaderFH = getSamHeaderFH(0,@$chromsRef,$isSkipVariableMetadata);
    my $tmpChrBamPath = File::Spec->catfile($binDir,$tmpChrBamFilename);

    my $sorted2BamCmd = "| $sorted2Bam --no-unmapped --bam-file $tmpChrBamPath --header-file $samHeaderFH";
    $sorted2BamCmd .= " --paired" if($isPaired);
    $sorted2BamCmd .= " --qlogodds" if($isSolexaQvals);

    open(my $bamFH, $sorted2BamCmd) || errorExit("ERROR: failure in export to bam conversion process: $sorted2BamCmd\n");
    my $sortFH = $sortInfo->{sortFH};
    print $bamFH $_ while(<$sortFH>);
    close($sortFH) || errorExit("ERROR: failure in sort/merge process: " . $sortInfo->{cmd} ."\n");
    close($bamFH) || errorExit("ERROR: failure in export to bam conversion process: $sorted2BamCmd\n");

    $samHeaderFH = undef;
}


1;
__END__
