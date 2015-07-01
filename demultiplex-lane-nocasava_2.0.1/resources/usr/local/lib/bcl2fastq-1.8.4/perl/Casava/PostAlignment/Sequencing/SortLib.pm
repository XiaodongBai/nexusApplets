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

Casava::PostAlignment::Sequencing::SortLib - Utility library for read sorting and merging.

=head1 SYNOPSIS

Sorts or merges sets of files into a destination file. The functions
assume that files contain column delimited data which need to be
numerically sorted on a particular column.

=head1 DESCRIPTION

These functions are designed to abstract the particular method used
for large scale sorting and merging, but at the moment the underlying
work is done by the *nix sort tool. Note that this will lead to
variable performance across different systems due to different
versions of this tool being available.

=head1 SUBROUTINES

sortFileList()
mergeFileList()

=cut

# POD for for each subroutines is right before the implementation
# More POD after __END__

package Casava::PostAlignment::Sequencing::SortLib;

use strict;
use warnings "all";
use Exporter 'import';

our @EXPORT_OK = qw(checkSortBin getConfiguredSortCmd sortFileList mergeFileList);


use File::Temp;

use Casava::Common::Log;
use Casava::Common::IOLib qw(executeCmd);



# The sort program should be defined in CONF_APP{cmdSort}, if not this is the fallback:
my $defaultSortBin = "sort";



#
# return the maximum available command line length (after taking out the length of environment variable definitions):
#
sub getMaxCmdLen() {
    my $pipe = "getconf ARG_MAX |";
    open(GETC, $pipe)
      or errorExit "ERROR: Can't open pipe: '$pipe'\n";
    my $argMax = <GETC>;
    close(GETC);
    $pipe = "env | wc -c |";
    open(ENVP, $pipe)
      or errorExit "ERROR: Can't open pipe: '$pipe'\n";
    my $envStringLen = <ENVP>;
    close(ENVP);
    chomp $argMax;
    chomp $envStringLen;
    return ($argMax - $envStringLen);
}


#
# return configured sort program:
#
sub getSortBin($) {
    my ( $appRef ) = @_;

    return (defined $appRef->{cmdSort} ? $appRef->{cmdSort} : $defaultSortBin);
}



#
# tests whether sort program exists and logs which version is found:
#
sub checkSortBin(\%) {
    my ( $appRef ) = @_;

    my $sortBin = getSortBin($appRef);

    #
    # test for existence of sort
    #
    my $testFH = File::Temp->new();
    my $cmd = "which $sortBin > $testFH 2> /dev/null";
    system($cmd);
    errorExit "ERROR: Can't find sort binary: '$sortBin'\n" unless ( $? == 0 );
    my $sortUsed = <$testFH>;
    chomp $sortUsed;
    logInfo( "Found sort binary: '$sortUsed'\n", 1 );

    #
    # test that this works too while we're here:
    #
    my $tmp = getMaxCmdLen();
}



#
# values setup by the sortInitialize function:
#
my $initProjRef;
my $isSortFiles0 = 0;
my $baseSortCmd = "";
my $sortCmd = "";
my $mergeCmd = "";
my $maxCmdLen = 0;



# pod TBD
sub sortInitialize($$) {
    my ( $projRef, $appRef ) = @_;

    my $sortBin = getSortBin($appRef);
    $maxCmdLen = getMaxCmdLen();

    #
    # derived sort values:
    #
    my $sortBatchFlag = defined $projRef->{sortMaxBatchSize} ? "--batch-size=$projRef->{sortMaxBatchSize}" : "";

    my $isSortBatchSize = 0;
    {
        my $testFH = File::Temp->new();
        my $cmd = "$sortBin $sortBatchFlag $testFH -o /dev/null 2> /dev/null";
        system($cmd);
        $isSortBatchSize = ( 0 == $? );
    }

    {
        my $testFH = File::Temp->new();
        print $testFH "$testFH\0";
        my $cmd = "$sortBin --files0-from=$testFH -o /dev/null 2> /dev/null";
        system($cmd);
        $isSortFiles0 = ( 0 == $? );
    }

    my $fieldSeparator  = "\t";
    my $sortTempFlag = defined $projRef->{dirBuildTemp} ? "-T $projRef->{dirBuildTemp}" : "";
    my $sortBufferFlag = "";
    if(defined $projRef->{sortBufferSize}){
        $sortBufferFlag = "--buffer-size=" . $projRef->{sortBufferSize} . "M";
    }
    $baseSortCmd = "$sortBin $sortTempFlag $sortBufferFlag";

    if( $isSortBatchSize ) {
        $baseSortCmd .= " $sortBatchFlag";
    }
    $sortCmd = "$baseSortCmd -t \"$fieldSeparator\"";
    $mergeCmd = "$sortCmd -m";

    $initProjRef = $projRef;
}



# pod TBD
sub getConfiguredSortCmd(\%\%) {
    my ( $projRef, $appRef ) = @_;
    if( (not defined $initProjRef) or ( $initProjRef != $projRef ) ) {
        sortInitialize($projRef,$appRef);
    }

    return $baseSortCmd;
}



# this definition allows function to be called recursively:
sub sortFileListInt($$$$$$;$);



# pod TBD
sub sortFileListInt($$$$$$;$) {
    my ( $inFilesRef, $outFile, $isMerge, $sortByField, $projRef, $appRef, $outFilter ) = @_;

    if( (not defined $initProjRef) or ( $initProjRef != $projRef ) ) {
        sortInitialize($projRef,$appRef);
    }

    my $isPipeOut = ($outFile eq "-");
    my $isFilter = ($isPipeOut and (defined $outFilter));
    if((not $isPipeOut) and (defined $outFilter)) {
        errorExit("ERROR: outFilter cannot be applied to file output in sortFileList or mergeFileList\n");
    }

    my $cmdPrefix = 'LC_ALL=C ';
    my $cmdSuffix = '';
    if($isFilter) {
        $cmdPrefix.="bash -o pipefail -c '";
        $cmdSuffix.=" | $outFilter'";
    }
    if($isPipeOut) {
        $cmdSuffix.=" |";
    }
    my $localMaxCmdLen = $maxCmdLen - length($cmdPrefix) - length($cmdSuffix);

    my $cmd;
    my $sortListFH;

    if( scalar(@{$inFilesRef}) > 0 ) {
        my $baseCmd = $isMerge ? $mergeCmd : $sortCmd;
        $baseCmd .= " -nk$sortByField";
        $cmd = "$baseCmd ";
        $cmd .= "-o '$outFile' " unless($isPipeOut);
        if( length($cmd) > $localMaxCmdLen ) {
            errorExit "ERROR: Unable to perform sort/merge operation within the command-line length limit. Base cmd: $cmd\n";
        }

        # specify the input files:
        #
        my $inFileCmd;
        if ( $isSortFiles0 ) {
            $sortListFH = File::Temp->new();
            for my $file (@{$inFilesRef}) {
                print $sortListFH "$file\0";
            }
            $inFileCmd = "--files0-from=$sortListFH";
        } else {
            $inFileCmd = join(' ', (map {"'$_'"} @{$inFilesRef}));

            # Check to see whether the command line will exceed the
            # maximum length. If so then recursively merge from
            # smaller sort/merge batches:
            #
            if( length($cmd) + length($inFileCmd) > $localMaxCmdLen ) {

                my @mergeFiles = (File::Temp->new());
                $cmd = "$baseCmd -o '$mergeFiles[-1]'";
                my $cmdLen = length($cmd);
                my @inFilesSubset = ();
                my $subsetIfcLen = 0;
                for my $inFile (@{$inFilesRef}) {
                    my $ifcLen = length($inFile)+1;
                    if( ($cmdLen + $subsetIfcLen + $ifcLen) > $localMaxCmdLen ) {
                        if( scalar(@inFilesSubset) < 2) {
                            errorExit "ERROR: Unable to perform sort/merge operation within the command-line length limit. Base cmd: $cmd\n";
                        }
                        sortFileListInt(\@inFilesSubset,  $mergeFiles[-1], $isMerge, $sortByField, $projRef, $appRef);
                        push @mergeFiles, File::Temp->new();
                        $cmd = "$baseCmd -o '$mergeFiles[-1]'";
                        $cmdLen = length($cmd);
                        @inFilesSubset = ();
                        $subsetIfcLen = 0;
                    }
                    $subsetIfcLen += $ifcLen;
                    push @inFilesSubset, $inFile;
                    if( ($cmdLen + $subsetIfcLen) > $localMaxCmdLen ) {
                        errorExit "ERROR: Unable to perform sort/merge operation within the command-line length limit. Base cmd: $cmd\n";
                    }
                }
                sortFileListInt(\@inFilesSubset,  $mergeFiles[-1], $isMerge, $sortByField, $projRef, $appRef);
                return sortFileListInt(\@mergeFiles, $outFile, 1, $sortByField, $projRef, $appRef, $outFilter);
            }
        }
        $cmd .= "$inFileCmd";
    }
    else {
        if($isPipeOut) {
            $cmd = "echo -n";
        } else {
            $cmd = "touch '$outFile'";
        }
    }    # if

    if( length($cmd) > $localMaxCmdLen ) {
        errorExit "ERROR: attempting system call which exceeds the maximum command-line length: $cmd\n";
    }

    $cmd = $cmdPrefix . $cmd . $cmdSuffix;

    if($isPipeOut) {
        # note that client doesn't need to use tmpFiles, we simply
        # need all File::Temp objects to have the same lifetime as
        # sortFH or $cmd if they're to be of any use:
        open(my $sortFH,"$cmd") or errorExit("ERROR: sort/merge process failed: '$cmd'");
        return { sortFH => $sortFH,
                 cmd => $cmd,
                 tmpFiles => [ $sortListFH, $inFilesRef ] };
    }
    executeCmd( $cmd, 5 );
    return undef;
}


# returns a file handle of $outFile is '-'. Additionally if $outFilter
# is defined, this text will be inserted in a pipe on the sort output.
#
# pod TBD
sub sortFileList(\@$$\%\%;$) {
    my ( $inFilesRef, $outFile, $sortByField, $projRef, $appRef, $outFilter ) = @_;

    return sortFileListInt($inFilesRef, $outFile, 0, $sortByField, $projRef, $appRef, $outFilter);
}



# pod TBD
sub mergeFileList(\@$$\%\%;$) {
    my ( $inFilesRef, $outFile, $sortByField, $projRef, $appRef, $outFilter ) = @_;

    return sortFileListInt($inFilesRef, $outFile, 1, $sortByField, $projRef, $appRef, $outFilter);
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

