#!/usr/bin/env perl
# PROJECT: CASAVA
# MODULE:  $RCSfile: Reports.pm,v $
# AUTHOR:  Chris Saunders
#
# Copyright (c) 2008 Illumina
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).

package Casava::PostAlignment::Sequencing::Reports;

use warnings FATAL => 'all';
use strict;
use Exporter 'import';

our @EXPORT_OK = qw(getColIndex htmlFileMenuRegister updateTotalRow);



use File::Copy qw(move);

use Casava::Common::Log;



#
# given the key-value hash derived from the "#$" header lines in some
# CASAVA files, find the "COLUMNS" key and parse its value into a hash
# of column labels which return 0-indexed column numbers.
#
sub getColIndex(\%\%$) {
    my ($colIndex, $fileMD, $file) = @_;

    if(not defined $fileMD->{COLUMNS} ) {
        errorExit("ERROR: no column labels in file: $file\n");
    }

    my $col = 0;
    for my $colLabel (split(/\s+/,$fileMD->{COLUMNS})) {
        $colIndex->{$colLabel} = $col;
        $col++;
    }
}



#
# Attempts to register an html page in a javascript file
# included in the standard html template. The idea is that
# by updating here, you update the links on the navigation
# menu in all CASAVA html pages.
#
sub htmlFileMenuRegister($$$) {

    my ( $htmlPath , $htmlFileName, $label) = @_;

    my $jsPath = File::Spec->catdir( $htmlPath , "js" );
    my $menuFile = File::Spec->catfile( $jsPath , "menu.js" );
    return if(not -e $menuFile);

    my $is_start=0;
    my $is_end=0;
    my $newMenuFile = $menuFile . ".tmp";
    open(my $MFH, "< $menuFile") or errorExit("ERROR: can't open file: $menuFile\n");
    open(my $newFH, "> $newMenuFile") or errorExit("ERROR: can't open file: $newMenuFile\n");
    while(<$MFH>) {
        my $line=$_;
        if(not $is_end) {
            if($is_start) {
                if(/,.*'([^']*)'/){
                    if($1 ge $label){
                        if($1 ne $label) {
                            print $newFH " '$htmlFileName', '$label',\n";
                        }
                        $is_end=1;
                    }
                }
                if((not $is_end) and /;/){
                    print $newFH " '$htmlFileName', '$label',\n";
                    $is_end=1;
                }
            } else {
                if(/var registeredPages/){
                    $is_end=1 if(/;/);
                    $is_start=1;
                }
            }
        }
        print $newFH $line;
    }
    close($MFH);
    close($newFH);
    move($newMenuFile,$menuFile) or errorExit("ERROR: File move failed: $!\n");
}



#
# adds each row of data to a totalrow, skips addition
# in the name column
#
sub updateTotalRow($$$) {

    my ($trow,$row,$namecol) = @_;
    if(scalar(@$trow) != 0) {
        for(my $i=0;$i<scalar(@$row);++$i) {
            next if ( $i == $namecol );
            $trow->[$i] += $row->[$i];
        }
    } else {
        @$trow = @$row;
        $trow->[ $namecol ] = "Total";
    }
}



1;
__END__
