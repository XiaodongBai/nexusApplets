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

Casava::PostAlignment::Sequencing::SamLib - Utility library for SAM/BAM file format operations 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SUBROUTINES

=cut

# POD for for each subroutines is right before the implementation
# More POD after __END__

package Casava::PostAlignment::Sequencing::SamLib;

use strict;
use warnings "all";
use Exporter 'import';

our @EXPORT_OK = qw(getSamtoolsBin checkSamtoolsBin changeChrLabel
                    getChangeChrLabelType getRefSizeFH getSamHeaderFH
                    matchDescToCigar);


use File::Temp;

use Casava::Common::Log;
use Casava::Common::IOLib qw(executeCmd);
use Casava::PostAlignment::Sequencing::Config 
  qw(isSpliceJunctionChrom getRunConf %chrEnds %CONF_PROJ %CONF_RUNS);

# The samtools program should be defined in CONF_PROJ{cmdSamtools}, if not this is the fallback:
my $defaultSamtoolsBin = "samtools";


#
# return configured samtools program:
#
sub getSamtoolsBin(\%) {
    errorExit "ERROR: getSamtoolsBin wrong parameters\n" unless ( @_ == 1 );
    my ( $confRef ) = @_;

    return (defined $confRef->{cmdSamtools} ? $confRef->{cmdSamtools} : $defaultSamtoolsBin);
}



#
# tests whether sort program exists and logs which version is found:
#
sub checkSamtoolsBin(\%) {
    errorExit "ERROR: checkSamtoolsBin wrong parameters\n" unless ( @_ == 1 );
    my ( $confRef ) = @_;

    my $samtoolsBin = getSamtoolsBin(%{$confRef});

    #
    # test for existence of sort
    #
    my $testFH = File::Temp->new();
    my $cmd = "which $samtoolsBin > $testFH 2> /dev/null";
    system($cmd);
    errorExit "ERROR: Can't find samtools binary: '$samtoolsBin'\n" unless ( $? == 0 );
    my $samtoolsUsed = <$testFH>;
    chomp $samtoolsUsed;
    logInfo( "Found samtools binary: '$samtoolsUsed'\n", 1 );
}



#
# get chromosome label change code
#
sub getChangeChrLabelType($)
{
    my($changeStr) = @_;

    if      ((not defined $changeStr) or ($changeStr eq "") or ($changeStr eq "OFF")){
        return 0;
    } elsif ($changeStr eq "NOFA"){
        return 1;
    } elsif ($changeStr eq "UCSC"){
        return 2;
    } else {
        errorExit("ERROR: Unknown bam chromosome label change code: $changeStr\n");
    }
}



#
# change chromosome labels:
#
sub changeChrLabel($$)
{
    my ($chrName,$changeType) = @_;

    if      ($changeType == 0){
        return $chrName;
    } elsif ($changeType == 1){
        # remove standard CASAVA fasta file extension if present:
        $chrName =~ s/\.fa$//;
        return $chrName;
    } elsif ($changeType == 2){
        $chrName =~ s/\.fa$//;
        if      ($chrName =~ /(\d+|[XY])$/) {
            return "chr$1";
        } elsif ($chrName =~ /M[Tt]?$/) {
            return "chrM";
        } else {
            return $chrName;
        }
    } else {
        errorExit("ERROR: Unknown bam chromosome label change type: $changeType\n");
    }
}



# get tab-delimited chromosome name and size temporary file required for SAM->BAM conversion
#
sub getRefSizeFH($\@)
{
    my ( $bamChangeChrType, $buildChromsRef ) = @_;

    # get a temporary file-handle to a tab delimited
    # chromosome/chromosome-size file appropriate for bam->sam
    # conversion
    #
    my $refSizeFH = File::Temp->new();
    {
        # order references in BAM file by chromosome size:
        my @tmpChromList = sort { $chrEnds{$b}->{length} <=> $chrEnds{$a}->{length} } @{$buildChromsRef};
        for my $tmpChrom (@tmpChromList) {
            next if ( isSpliceJunctionChrom($tmpChrom, %CONF_PROJ) );
            my $printChrom = changeChrLabel($tmpChrom,$bamChangeChrType);
            printf $refSizeFH "%s\t%i\n", $printChrom, $chrEnds{$tmpChrom}->{length};
        }
    }
    return $refSizeFH;
}



# Get full SAM header for sorted->BAM conversion
#
sub getSamHeaderFH($\@$) {

    my ( $bamChangeChrType, $buildChromsRef, $isSkipVariableMetadata ) = @_;

    my $FH = File::Temp->new();

    # all of the arguments passed to this script must be coordinate
    # sorted:
    print $FH join("\t",'@HD',"VN:1.0","SO:coordinate"),"\n";

    if(not $isSkipVariableMetadata) {
        my $prog_name = 'CASAVA';
        my $prog_version = 'bcl2fastq-1.8.4';

        my @sam_prog_parts = ('@PG', "ID:${prog_name}", "VN:${prog_version}");

        my %CONF_RUN =
          %{getRunConf($CONF_PROJ{currentRunId}, %CONF_RUNS, %CONF_PROJ)};
        my $arg_str = $CONF_RUN{argv};

        my $cmd_line
          = join(' ', File::Spec->catfile('/usr/local/bin', 'configureBuild.pl'),
                 $arg_str);

        push(@sam_prog_parts, "CL:${cmd_line}");

        print $FH join("\t", @sam_prog_parts), "\n";
    }

    # order references in BAM file by chromosome size:
    my @chromList = sort { $chrEnds{$b}->{length} <=> $chrEnds{$a}->{length} } @{$buildChromsRef};
    for my $chrom (@chromList) {
        next if ( isSpliceJunctionChrom($chrom, %CONF_PROJ) );
        my $printChrom = changeChrLabel($chrom,$bamChangeChrType);
        printf $FH "\@SQ\tSN:%s\tLN:%i\n", $printChrom, $chrEnds{$chrom}->{length};
    }

    return $FH;
}


#
# richard shaw's function originally from sorted2sam
#
sub matchDescFragLength($)
{
    my ($match_desc) = @_;
    my $len = 0;

    my @match_desc_fields = split(/([ACGTN]+)/, $match_desc);

    for my $match_desc_field (@match_desc_fields) {
        next if ($match_desc_field eq '');

        $len += (($match_desc_field =~ /(\d+)/)
                 ? $1 : length($match_desc_field));
    }

    return $len;
}



#
# richard shaw's function originally from sorted2sam
#
sub matchDescToCigar($)
{
    my ($match_desc) = @_;

    my @match_desc_parts = split(/(\^.*?\$)/, $match_desc);
    my $cigar_str = '';
    my $cigar_del_ch = 'D';
    my $cigar_ins_ch = 'I';
    my $cigar_match_ch = 'M';

    foreach my $match_desc_part (@match_desc_parts) {
        next if (!$match_desc_part);

        if ($match_desc_part =~ /^\^([ACGTN]+)\$$/) {
            # Deletion
            $cigar_str .= (length($1) . $cigar_del_ch);
        } elsif ($match_desc_part =~ /^\^(\d+)\$$/) {
            # Insertion
            $cigar_str .= ($1 . $cigar_ins_ch);
        } else {
            $cigar_str .= (matchDescFragLength($match_desc_part)
                           . $cigar_match_ch);
        }
    }

    return $cigar_str;
}



# example function pod:
=pod

=head2 configure

Load the 'config.xml' and generate the 'Makefile.config'.

I<Parameters:>

=over 4

=item *

$configFile

full path to the input XML config file for the analysis folder

=item *

$output

reference to the output FILE HANDLE (defaults to \*STDOUT)

I<Returns:>

Nothing

I<Exceptions:>

=over 4

=item *

Config file does not exist

=item *

Missing elements in the XML configuration file

=item *

Failing to write to the indicated output

=back

=cut








1;
__END__

=pod

=head1 CONFIGURATION AND ENVIRONMENT

Name and location of configuration files, with a complete description of the properties that can be set.

Name and description of the relevant environment variables

=head1 DEPENDENCIES

*nix sort program with an interface conforming to that gnu coreutils

=over 4

=item Standard perl modules

strict, warnings, Exporter

=item External perl modules

=item Casava perl modules



=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Chris Saunders

=cut

