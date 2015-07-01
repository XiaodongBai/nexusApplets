#!/usr/bin/env perl

# File: plotIntensity.pl
# Copyright (c) 2005, 2006 Solexa; 2008 Illumina
# Author: A. J. Cox, L. J. Davies
# Extended by R. J. Shaw.
#
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).

# Description:
# Produces intensity-by-cycle plots from intensity data files.

#------------------------------------------------------------------------------

use warnings;
use strict;
use File::Spec;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Common::Utils;


my $plot_data_ext_str = 'dat';
my $plot_script_ext_str = 'txt';
my $tmpdir  = "Temp";
my $plotdir = "Plots";


sub write_intensity_data_file($$$$$)
{
    my ($plot_tmp_name, $intensities_ref, $num_rows,
        $max_is_set_ref, $max_ref) = @_;

    my $max_is_set = $$max_is_set_ref;
    my $max = $$max_ref;

    my $data_filename = join('.', $plot_tmp_name, $plot_data_ext_str);

    open (OUT, ">$data_filename")
    || die "Problem opening $data_filename for writing: $!";

    my $num_cols = 5;

    if ($num_rows > 0) {
        for (my $row_ind = 0; $row_ind < $num_rows; ++$row_ind) {
            for (my $col_ind = 0; $col_ind < $num_cols; ++$col_ind) {
                my $curr_intensity
                    = $intensities_ref->[$row_ind][$col_ind + 1];
                print OUT $curr_intensity, "\t";  

                if ((!$max_is_set) && ($curr_intensity > $max)) {
                    $max = $curr_intensity;
                }
            }

            print OUT "\n";
        }
    } else {
        # No data - either no clusters or none that passed quality filtering.
        # Write a single row of zeros to keep downstream utilities happy.

        for (my $col_ind = 0; $col_ind < $num_cols; ++$col_ind) {
            print OUT "0\t";
        }

        print OUT "\n";
    }

    close (OUT);

    if (!$max_is_set) {
        $$max_ref = (int($max * 1.5 / 1000) + 1) * 1000;

        # Now that all tiles in a lane are handled in the same instance of this
        # script, they could all be plotted on the same scale (set by the first
        # tile) as appeared to be the original purpose of the setmax variable.
        # Uncomment the line below to enable this behaviour.
        # ++$$max_is_set;
    }
}

#------------------------------------------------------------------------------

sub draw_intensity_graph($$$$)
{
    my ($plot_tmp_name, $plot_image_filename, $title, $max) = @_;

    my $plot_script_name = join('.', $plot_tmp_name, $plot_script_ext_str);
    my $data_filename = join('.', $plot_tmp_name, $plot_data_ext_str);

    my $width = 3;

    open (OUT, ">$plot_script_name") 
        || die "Problem opening $plot_script_name for writing: $!";

    print OUT << "END";
#set terminal png color;
set terminal postscript portrait color solid;
set yrange [0:$max];
set ytics 500;
    
set size 1.0, 1.0;
set origin 0.0, 0.0;
set title "$title";
set xlabel "Cycle Number";
set ylabel "Intensity";
plot '$data_filename' using 1:2 title 'A' with linespoints lt 9 lw $width pt 7 ps 2, '$data_filename' using 1:3 title 'C' with linespoints lt 1 lw $width pt 7 ps 2,    '$data_filename' using 1:4 title 'G' with linespoints lt 3 lw $width pt 7 ps 2, '$data_filename' using 1:5 title 'T' with linespoints lt 2 lw $width pt 7 ps 2;  

END

    close (OUT);
    gnuplotImage($plot_script_name, $plot_image_filename);
    unlink $data_filename, $plot_script_name;
}

#------------------------------------------------------------------------------

my $app_name_len = length($0);

my $usage
    = "Usage: $0 <base_dir> <tiles.txt> <s_1> <intensity_path> <_all.txt> <_all.png>\n\n"
    . "For each tile in the specified tile list file that matches the lane\n"
    . "prefix, reads the corresponding intensity file under <intensity_path>\n"
    . "(with suffix <_all.txt>) and produces the intensity graph (in an\n"
    . "image file with suffix <_all.png>)\n";

die $usage unless (@ARGV == 6);

my ($base_dir, $tile_file_path, $lane_prefix, $intensity_dir_str,
    $all_text_suffix, $all_png_suffix) = @ARGV;

if ($base_dir ne ".") {
    $tile_file_path = File::Spec->catfile($base_dir,$tile_file_path);
    $intensity_dir_str = File::Spec->catfile($base_dir,$intensity_dir_str);
    $tmpdir  = File::Spec->catdir($base_dir,"Temp");
    $plotdir = File::Spec->catdir($base_dir,"Plots");
}

my @tile_list;

if (!get_tiles_by_lane($tile_file_path, $lane_prefix, @tile_list)) {
    die("Failed to get tile list from file `$tile_file_path'.");
}

my $plot_tmp_name = File::Spec->catfile($tmpdir, "tmp-plotIntensity-$>-$$");
my $max_is_set = 0;

foreach my $tile (@tile_list) {
    my $all_text_filename = $tile . $all_text_suffix;
    my $plot_image_filename = File::Spec->catfile($plotdir,
                                                  $tile . $all_png_suffix);

    my $all_text_path_str 
        = File::Spec->catfile($intensity_dir_str, $all_text_filename);
    my @intensities;

    open(IN, $all_text_path_str)
        || die "Problem opening $all_text_path_str: $!";

    while (<IN>) {
        next if ($_ =~ m/\#/);
        push (@intensities, [split('\t', $_)]);
    }

    close(IN);

    my $num_rows = scalar(@intensities);
    my $max = 0;
    write_intensity_data_file($plot_tmp_name, \@intensities, $num_rows,
                              \$max_is_set, \$max);

    my $title = $tile . "_all";
    draw_intensity_graph($plot_tmp_name, $plot_image_filename, $title, $max);
}

#------------------------------------------------------------------------------










