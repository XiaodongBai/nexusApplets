#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2007-2010 Illumina, Inc.

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

transformFasta.pl - transform reference sequence

=head1 SYNOPSIS

transformFasta.pl <options> file1|- [file2]...

Options:

=over 4

=item *

--help

brief help message

=item *

--man

full documentation

=item *

--targetPath=s

Path where to put the result files. if '-', the result data is printed
into standard output

=item *

--chromNameSource=fileName|contigName

Controls how the result contig and file names are formed.

=item *

--removeBreaks

If set, will not preserve line breaks within the data.

=item *

--removeHeaders

If set, will remove fasta header lines from the output. Indendex for use with
--removeBreaks to produce one separate line of output per each contig.

=back

=head1 DESCRIPTION

Transforms input collection of fasta files.

=head1 DIAGNOSTICS

=head2 Exit status

0: successful completion
1: abnormal completion
2: fatal error

=head2 Errors

All error messages are prefixed with "ERROR: ".

=head2 Warnings

All warning messages generated by CASAVA are prefixed with "WARNING: ".

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Roman Petrovski

=cut


use warnings "all";
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Path qw(mkpath);

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl';
use Casava::Common::Log qw(errorExit initLog);
use Casava::PostAlignment::Sequencing::Config qw(loadConfiguration);


sub checkMakeDir($) {
    my $dir = shift;
    unless (-e $dir) {
        mkpath($dir) || errorExit("ERROR: Can't create directory '$dir'\n");
    } else {
        errorExit "ERROR: Path is not a directory '$dir'\n"  unless -d $dir;
    }
}


initLog( "", 0, 5 );

$ENV{PATH} = $1 if $ENV{PATH}=~/^(.*)$/;

my $man = 0;
my $help = 0;
my $CFG_TARGET_PATH = undef;
my $chromNameSource = undef;
my $removeBreaks;
my $removeHeaders;

GetOptions('targetPath=s'           => \$CFG_TARGET_PATH,
           'chromNameSource=s'      => \$chromNameSource,
           'removeBreaks'           => \$removeBreaks,
           'removeHeaders'           => \$removeHeaders,
           'help|?' => \$help, 
           'man' => \$man) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2,  -input => $1) if ($man and $0 =~ /(.*)/);
pod2usage("ERROR: --targetPath not specified") if ( !defined $CFG_TARGET_PATH );
pod2usage("ERROR: --chromNameSource not specified") if ( !defined $chromNameSource );


my $inSequence = 1;
my $seqName = '';
my $lastOutputPath = '';
my $sequenceBreak = ($removeBreaks) ? '' : "\n";
my $headerBreak = '';

checkMakeDir($CFG_TARGET_PATH) if($CFG_TARGET_PATH ne '-');

while (@ARGV)
{
    open (SEQINPUT, "<$ARGV[0]") or errorExit "ERROR: Unable to open input stream: $!";
    my $inputFileName = (File::Spec->splitpath(shift (@ARGV)))[2];

    while (<SEQINPUT>)
    {
        chomp;
        my $line = $_;
        if ($line =~ m/^>([^\s]*)/)
        {
            next if (!$inSequence);
            $seqName = $1;
            $seqName = $inputFileName if ('fileName' eq $chromNameSource);
    
            errorExit "ERROR: Expected fasta header entry with name. got: $line" unless length($seqName);
    
            $inSequence = 0;
            my $outputPath = ('-' eq $CFG_TARGET_PATH) ? '-' 
                : File::Spec->catfile($CFG_TARGET_PATH, $seqName);
            errorExit "ERROR: Output file $outputPath already exists." 
                if ('contigName' eq $chromNameSource and '-' ne $CFG_TARGET_PATH and -e $outputPath);
            if ($outputPath ne $lastOutputPath)
            {
                open (SEQCACHE, ">$outputPath") or errorExit "ERROR: Unable to create output stream on: $outputPath: $!";
                $lastOutputPath = $outputPath;
                $headerBreak = ''; 
            }
            else
            {
                #if we remove breaks from data, we need to print them before second header
                $headerBreak = "\n" if $removeBreaks; 
            }
            if (!$removeHeaders)
            {
                print SEQCACHE ('fileName' eq $chromNameSource) ? "${headerBreak}>${inputFileName}\n" : "${headerBreak}${line}\n"
                    || errorExit "ERROR: Unable to write the results: $!";
            }
            elsif (${headerBreak})
            {
                print SEQCACHE "${headerBreak}" || errorExit "ERROR: Unable to write the results: $!";
            }
        }
        else
        {
            errorExit "ERROR: Expected fasta header entry with name. got: $line"
                unless length($seqName);
            
            $inSequence = 1;
            print (SEQCACHE "${line}${sequenceBreak}") || errorExit "ERROR: Unable to write the results: $!";
        }
    
    }
}
close (SEQCACHE);
1;
__END__