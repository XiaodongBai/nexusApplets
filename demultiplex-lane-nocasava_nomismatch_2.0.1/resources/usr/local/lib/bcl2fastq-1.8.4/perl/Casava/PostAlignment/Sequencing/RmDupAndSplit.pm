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

Casava::PostAlignment::Sequencing::RmDupAndSplit -

=head1 SYNOPSIS

duplicate removal

=head1 DESCRIPTION

duplicate removal

=head1 SUBROUTINES

rmDupAndSplitPairFiles

Process the unsorted exportpair files for the norm, orph and anom
cases to: (1) sort by pair start position (2) remove PCR duplicates
(3) split the paired read record into individual read records,
preserving the sorted state for a subset of the reads where
possible. No intermediate files are written out during these steps,
except temporaries used by the sort program.

=cut

# POD for for each subroutines is right before the implementation
# More POD after __END__

package Casava::PostAlignment::Sequencing::RmDupAndSplit;

use strict;
use warnings "all";
use Exporter 'import';

our @EXPORT_OK = qw(&rmDupAndSplitPairFiles);


use Data::Dumper;
use IO::File;

use Casava::Common::Log;
use Casava::PostAlignment::Sequencing::Config qw(%CONF_PROJ %CONF_APP %runsConfig);
use Casava::PostAlignment::Sequencing::DnaSeqLib qw(
       &splitSortNormInit &splitSortNormRecord &splitSortNormClose
       &splitSortOrphInit &splitSortOrphRecord &splitSortOrphClose
       &splitSortAnomInit &splitSortAnomRecord &splitSortAnomClose);
use Casava::PostAlignment::Sequencing::GenomicIOLib qw(%doubleExportFields);
use Casava::PostAlignment::Sequencing::SortLib qw(sortFileList);

my $fieldsRef = \%doubleExportFields;



sub makeDupReport($\%) {
    my ($label,$dupStatsRef) = @_;
    return sprintf "After processing %s, out of all %d reads:\n"
      . "\t%d were not marked as duplicates (%.1f%%)\n"
      . "\t%d were marked as duplicates (%.1f%%)\n",
        $label, $dupStatsRef->{total},
        $dupStatsRef->{nonDupCount}, ( $dupStatsRef->{nonDupCount} / $dupStatsRef->{total} ) * 100,
        $dupStatsRef->{dupCount}, ( $dupStatsRef->{dupCount} / $dupStatsRef->{total} ) * 100;
}



sub getDupDiff(\%) {
    my ($sref) = @_;
    return $sref->{total} - ( $sref->{nonDupCount} + $sref->{dupCount} );
}



sub checkDupStats($\%$) {
    my ($label,$dupStatsRef,$isRmDup) = @_;
    return if ((not $isRmDup) or ($dupStatsRef->{total} <= 0));
    my $difference = getDupDiff(%$dupStatsRef);
    if($difference != 0){
        my $dupReport = makeDupReport($label,%$dupStatsRef);
        errorExit("ERROR: wrong number of duplicates by $difference\n.$dupReport");
    }
}



sub stringToFile($$) {
    my ($str,$filename) = @_;
    my $FH = IO::File->new( ">$filename" )
      || errorExit "ERROR: Couldn't create/open file handle for $filename $!\n";
    print $FH $str;
    close $FH;
}



# SUSPENDED:
#
# # For back-compatibility, we prevent unanchored read-pairs from being
# # selected when an anchored pair is available, even if the anchored pair
# # has a lower PE score:
#     (($_[1]->[$fieldsRef->{SingleScore}]+$_[1][$fieldsRef->{read2SingleScore}])==0) <=>
#         ($_[0]->[$fieldsRef->{SingleScore}]+$_[0][$fieldsRef->{read2SingleScore}]) ||





# Pick the non-filtered read with the best alignment score. Order
# candidates by other fields as well to be consistent across different
# CASAVA builds.
#
sub norm_sorter($$) {
    $_[1]->[$fieldsRef->{Filter}] cmp $_[0]->[$fieldsRef->{Filter}] ||
    $_[1]->[$fieldsRef->{PairedScore}] <=> $_[0]->[$fieldsRef->{PairedScore}] ||
    $_[0]->[$fieldsRef->{Seq}] cmp $_[1]->[$fieldsRef->{Seq}] ||
    $_[0]->[$fieldsRef->{read2Seq}] cmp $_[1]->[$fieldsRef->{read2Seq}] ||
    $_[0]->[$fieldsRef->{QualityString}] cmp $_[1]->[$fieldsRef->{QualityString}] ||
    $_[0]->[$fieldsRef->{read2QualityString}] cmp $_[1]->[$fieldsRef->{read2QualityString}];
}



sub pickBestSrtNorm(\%\%\%$) {
    my ( $ref2Store, $dupStatsRef, $splitSortNormStateRef, $isKeepDupReads ) = @_;

    for (sort keys %$ref2Store) {
        my $reads = $ref2Store->{$_};
        my $is_first=1;
        for my $read (sort norm_sorter @$reads) {
            if($is_first) {
                $dupStatsRef->{nonDupCount}++;
                splitSortNormRecord($read,$splitSortNormStateRef);
                $is_first=0;
            } else {
                $dupStatsRef->{dupCount}++;
                splitSortNormRecord($read,$splitSortNormStateRef,"DUP") if($isKeepDupReads);
            }
        }
        $dupStatsRef->{dist}{norm} {scalar(@$reads)}++;
    }
    %$ref2Store = ();
}



# Pick the non-filtered read with the best single alignment
# score. Order candidates by other fields as well to be consistent
# across different CASAVA builds.
sub orph_sorter($$) {
    $_[1]->[$fieldsRef->{Filter}] cmp $_[0]->[$fieldsRef->{Filter}] ||
    $_[1]->[$fieldsRef->{SingleScore}] <=> $_[0]->[$fieldsRef->{SingleScore}] ||
    $_[0]->[$fieldsRef->{Seq}] cmp $_[1]->[$fieldsRef->{Seq}] ||
    $_[0]->[$fieldsRef->{QualityString}] cmp $_[1]->[$fieldsRef->{QualityString}]
}



sub pickBestSrtOrph(\%\%\%$) {
    my ( $ref2Store, $dupStatsRef, $splitSortOrphStateRef, $isKeepDupReads ) = @_;

    for (sort keys %$ref2Store) {
        my $reads = $ref2Store->{$_};
        my $is_first=1;
        for my $read (sort orph_sorter @$reads) {
            if($is_first) {
                $dupStatsRef->{nonDupCount}++;
                splitSortOrphRecord($read,$splitSortOrphStateRef);
                $is_first=0;
            } else {
                $dupStatsRef->{dupCount}++;
                splitSortOrphRecord($read,$splitSortOrphStateRef,"DUP") if($isKeepDupReads);
            }
        }
        $dupStatsRef->{dist}{orph} { scalar(@$reads) }++;
    }
    %$ref2Store = ();
}



# Pick the non-filtered read with the best sum single alignment
# score. Order candidates by other fields as well to be consistent
# across different CASAVA builds.
sub anom_sorter($$) {
    $_[1]->[$fieldsRef->{Filter}] cmp $_[0]->[$fieldsRef->{Filter}] ||
    ($_[1]->[$fieldsRef->{SingleScore}]+$_[1][$fieldsRef->{read2SingleScore}]) <=>
      ($_[0]->[$fieldsRef->{SingleScore}]+$_[0][$fieldsRef->{read2SingleScore}]) ||
    $_[0]->[$fieldsRef->{Seq}] cmp $_[1]->[$fieldsRef->{Seq}] ||
    $_[0]->[$fieldsRef->{read2Seq}] cmp $_[1]->[$fieldsRef->{read2Seq}] ||
    $_[0]->[$fieldsRef->{QualityString}] cmp $_[1]->[$fieldsRef->{QualityString}] ||
    $_[0]->[$fieldsRef->{read2QualityString}] cmp $_[1]->[$fieldsRef->{read2QualityString}];
}



sub pickBestSrtAnom(\%\%\%$) {
    my ( $ref2Store, $dupStatsRef, $splitSortAnomStateRef, $isKeepDupReads ) = @_;

    for (sort keys %$ref2Store) {
        my $reads = $ref2Store->{$_};
        my $is_first=1;
        for my $read (sort anom_sorter @$reads) {
            if($is_first) {
                $dupStatsRef->{nonDupCount}++;
                splitSortAnomRecord($read,$splitSortAnomStateRef);
                $is_first=0;
            } else {
                $dupStatsRef->{dupCount}++;
                splitSortAnomRecord($read,$splitSortAnomStateRef,"DUP") if($isKeepDupReads);
            }
        }
        $dupStatsRef->{dist}{anom} { scalar(@$reads) }++;
    }
    %$ref2Store = ();
}



{
    my $sortByField = 1;

    sub getSortFH($$$$) {
        my ( $chrom, $binId, $fileType, $isCompressPair ) = @_;

        my $u_file;
        if       ($fileType eq 'norm'){
            $u_file = $CONF_APP{f_unsort};
        } elsif($fileType eq 'anom') {
            $u_file = $CONF_APP{f_unsort_anom};
        } elsif($fileType eq 'orph') {
            $u_file = $CONF_APP{f_unsort_orph};
        } else {
            errorExit("ERROR: Invalid fileType argument to getSortFH '$fileType'\n");
        }

        my $dirBuildExportSets = File::Spec->catdir( $CONF_PROJ{dirBuildExport}, 'sets' );
        my @fileList;
        for my $i (1 .. scalar(@{$runsConfig{exportFiles}})) {
            my $fileUnsorted = File::Spec->catfile($dirBuildExportSets, $i, $chrom, $binId, $u_file);
            if ( -e $fileUnsorted ) {
                push @fileList, $fileUnsorted;
            }
        }

        my $outFilter;
        if($isCompressPair) {
            my $libexecDir = '/usr/local/libexec/bcl2fastq-1.8.4';
            my $cxp = File::Spec->catfile($libexecDir,'compressXPair');
            $outFilter = "$cxp -d";
        }
        return sortFileList(@fileList,'-',$sortByField,%CONF_PROJ,%CONF_APP,$outFilter);
    }
}



sub rmDupAndSplitPairFiles($$$$) {
    my ($sourceDir, $resultsDir, $chrom, $binId) = @_;

    # conf values
    #
    my $binSize     = $CONF_PROJ{binSizeBuild};
    my $isRmDup     = ( defined $CONF_PROJ{rmDup} ) ? ( $CONF_PROJ{rmDup} eq "YES" ) : 1;
    my $isKeepAllReads = ((defined $CONF_PROJ{sortKeepAllReads}) and $CONF_PROJ{sortKeepAllReads});
    my $isCompressPair = (not ((defined $CONF_PROJ{sortNoCompressPair}) and ($CONF_PROJ{sortNoCompressPair})));
    my $isKeepDupReads = $isKeepAllReads;

    # duplicate statistics:
    #
    my %dupStats = (
                    dupCount => 0,
                    total => 0,
                    nonDupCount => 0,
                    dist => {}
                   );

    #################################################### Norm Case

    my $sortInfo = getSortFH($chrom,$binId,'norm',$isCompressPair);
    if (defined $sortInfo) {
        my %STORE;
        my $lastRead;
        my %splitSortNormState;
        splitSortNormInit($binSize, $resultsDir, $fieldsRef, %splitSortNormState);

        my $sortFH = $sortInfo->{sortFH};
        while (<$sortFH>) {
            $dupStats{total}++;
            chomp;
            my ($thisRead, $thisEnd, @temp) = split("\t");
            ## TODO: Non rmDup case needs to be handled more efficiently
            if ( not $isRmDup ) {
                splitSortNormRecord(\@temp,\%splitSortNormState);
            }
            else {
                if ( ( defined $lastRead ) && ( $thisRead != $lastRead ) ) {
                    errorExit("ERROR: unexpected order") if($thisRead < $lastRead);
                    pickBestSrtNorm( %STORE, %dupStats, %splitSortNormState, $isKeepDupReads );
                }
                push @{ $STORE{$thisEnd} }, \@temp;
                $lastRead = $thisRead;
            }
        }
        close($sortFH) || errorExit("ERROR: sort process failed: $!\n");
        $sortInfo = undef;

        if ( $isRmDup and (defined $lastRead) ) {
            pickBestSrtNorm( %STORE, %dupStats, %splitSortNormState, $isKeepDupReads );
        }
        splitSortNormClose(%splitSortNormState);
    }
    checkDupStats("normals",%dupStats,$isRmDup);


    #################################################### Orph Case

    $sortInfo = getSortFH($chrom,$binId,'orph',$isCompressPair);
    if (defined $sortInfo) {
        my %STORE;
        my $lastRead;
        my %splitSortOrphState;
        splitSortOrphInit($binSize, $resultsDir, $fieldsRef, %splitSortOrphState);

        my $sortFH = $sortInfo->{sortFH};
        while (<$sortFH>) {
            $dupStats{total}++;
            chomp;
            my ($thisRead, $thisEnd, @temp) = split("\t");
            if ( not $isRmDup ) {
                splitSortOrphRecord(\@temp, \%splitSortOrphState);
            }
            else {
                if ( ( defined $lastRead ) && ( $thisRead != $lastRead ) ) {
                    errorExit("ERROR: unexpected order") if($thisRead < $lastRead);
                    pickBestSrtOrph( %STORE, %dupStats, %splitSortOrphState, $isKeepDupReads );
                }
                my $strand  = $temp[ $fieldsRef->{Strand} ];
                my $mateSeq = substr( $temp[ $fieldsRef->{read2Seq} ], 0, 8 );
                my $key = "$thisEnd|$strand|$mateSeq";
                push @{ $STORE{$key} }, \@temp;
                $lastRead = $thisRead;
            }
        }
        close($sortFH) || errorExit("ERROR: sort process failed: $!\n");
        $sortInfo = undef;

        if ( $isRmDup and (defined $lastRead) ) {
            pickBestSrtOrph( %STORE, %dupStats, %splitSortOrphState, $isKeepDupReads );
        }
        splitSortOrphClose(%splitSortOrphState);
    }
    checkDupStats("orphans",%dupStats,$isRmDup);


    ################################################# Anom Case

    $sortInfo = getSortFH($chrom,$binId,'anom',$isCompressPair);
    if (defined $sortInfo) {
        my %STORE;
        my $lastRead;
        my %splitSortAnomState;
        splitSortAnomInit($chrom, $binSize, $resultsDir, $fieldsRef, %splitSortAnomState);
        my $sortFH = $sortInfo->{sortFH};
        while (<$sortFH>) {
            $dupStats{total}++;
            chomp;
            my ($thisRead, $thisEnd, @temp) = split("\t");
            if ( not $isRmDup ) {
                splitSortAnomRecord(\@temp, \%splitSortAnomState);
            }
            else {
                if ( ( defined $lastRead ) && ( $thisRead != $lastRead ) ) {
                    errorExit("ERROR: unexpected order") if($thisRead < $lastRead);
                    pickBestSrtAnom( %STORE, %dupStats, %splitSortAnomState, $isKeepDupReads );
                }
                my $key = $thisEnd;
                if ($key !~ /\:/){
                    my $strand1 = $temp[ $fieldsRef->{Strand} ];
                    my $strand2 = $temp[ $fieldsRef->{read2Strand} ];

                    if ( $strand1 ne $strand2 )
                    {   # if we have an RF or FR treat them as the same thing
                        # they are either an oversized pair or an outie
                        $key .= "|FR";
                    } elsif (($strand1 eq "F") and ($strand2 eq "F")) {
                        $key .= "|FF";
                    } elsif (($strand1 eq "R") and ($strand2 eq "R")) {
                        $key .= "|RR";
                    } else {
                        errorExit "ERROR: Read pair not currently accounted for:\n". join("\t",@temp) . "\n";
                    }
                }
                push @{ $STORE{$key} }, \@temp;
                $lastRead = $thisRead;
            }
        }
        close($sortFH) || errorExit("ERROR: sort process failed: $!\n");
        $sortInfo = undef;

        if ( $isRmDup and (defined $lastRead) ) {
            pickBestSrtAnom( %STORE, %dupStats, %splitSortAnomState, $isKeepDupReads );
        }
        splitSortAnomClose(%splitSortAnomState);
    }

    if ( $dupStats{total} > 0 ) {
        my $dupReport = makeDupReport("anoms",%dupStats);
        my $fileDupCount = File::Spec->catfile($sourceDir, "dupCount.txt");
        stringToFile($dupReport,$fileDupCount);

        my $distReport = Dumper($dupStats{dist});
        my $fileDupDist  = File::Spec->catfile($sourceDir, "dupDist.txt");
        stringToFile($distReport,$fileDupDist);

        my $difference = getDupDiff(%dupStats);
        if ( $difference != 0 and $isRmDup ) {
            errorExit("ERROR: wrong number of duplicates by $difference\n.$dupReport");
        }
    }
}



1;
__END__

=pod

=head1 CONFIGURATION AND ENVIRONMENT

Name and location of configuration files, with a complete description of the properties that can be set.

Name and description of the relevant environment variables

=head1 DEPENDENCIES



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

Richard Carter

=cut
