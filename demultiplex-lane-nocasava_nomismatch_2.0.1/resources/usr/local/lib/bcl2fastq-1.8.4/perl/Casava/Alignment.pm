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

Casava::Alignment - Utility library for the management of the alignment folder.

=head1 SYNOPSIS

To generate the sample configuration:

  use Casava::Alignment;
  Casava::Alignment::configureSamples($sampleSheet);

To generate the analysis configuration for each dataset:

  use Casava::Alignment;
  Casava::Alignment::configureAnalysis($configFile, $project, $sample, $lane, $barcode);

=head1 DESCRIPTION

This script is a rewrite of the original GERALD library originally
implemented for the Pipeline.

Library for the creation, initialization and configuration of the
Alignment folder.

Provides a number of alignment-related services for CASAVA;

=head1 SUBROUTINES

=cut

# POD for for each subroutines is right before the implementation
# More POD after __END__

package Casava::Alignment;

use strict;
use warnings "all";
use Carp;
use Exporter 'import';

use Data::Dumper;
use File::Basename;
use XML::Simple;

use lib '/usr/local/lib/bcl2fastq-1.8.4/perl'; # substituted by CMake during install
use Casava::Alignment::Config;
use Casava::Common::Log;
use Casava::Demultiplex::DemultiplexConfig;
use Casava::Demultiplex::SampleSheet;
use Casava::Demultiplex::SampleSheet::Csv;

our @EXPORT_OK = qw(configureSamples configureDataset initialize selfTest);

sub _makeShortName($$$\%$;$)
{
    my ($shortNameKey, $elandGenome, $elandGenomeMask, $knownShortnames, $defaultName, $suffix) = @_;
    $suffix = '' unless defined $suffix;
    my $shortName = $knownShortnames->{$shortNameKey};
    unless (defined $shortName)
    {
        my $uniqueNumber =keys(%$knownShortnames);
        my $humanReadableName = $elandGenomeMask;
        $humanReadableName = (File::Spec->splitpath($elandGenome))[2] if $humanReadableName =~ /[^\w\-\.]/;

        $humanReadableName = $defaultName if $humanReadableName =~ /[^\w\-\.]/;
        $shortName = sprintf("$humanReadableName$suffix\_%02d", $uniqueNumber);
        $knownShortnames->{$shortNameKey} = $shortName;
    }
    return $shortName;
}

sub _makeSpliceShortName($$$$$\%)
{
    my($lane, $barcode, $read, $sampleSheet, $config, $knownShortnames) = @_;
    my($project, $sample, $speices) = ($sampleSheet->projectId($lane, $barcode), 
                                       $sampleSheet->sampleId($lane, $barcode), 
                                       $sampleSheet->species($lane, $barcode));

    my $genomeAnnotation = $config->selectValue('ELAND_RNA_GENOME_ANNOTATION', $project, $sample, $lane, $barcode, $speices);
    my $elandGenome = $config->selectValue('ELAND_GENOME', $project, $sample, $lane, $barcode, $speices);
    my $elandGenomeMask = $config->selectValue('ELAND_GENOME_MASK', $project, $sample, $lane, $barcode, $speices);
    my $squashGenomeParams = $config->selectValue('SQUASH_GENOME_PARAMS', $project, $sample, $lane, $barcode, $speices);
    my $spliceSitesParams = $config->selectValue('SPLCSTS_PARAMS', $project, $sample, $lane, $barcode, $speices);
    my $chromNameSource = $config->selectValue('CHROM_NAME_SOURCE', $project, $sample, $lane, $barcode, $speices);
    my $readLength = $config->selectValue("READ_LENGTH$read", $project, $sample, $lane, $barcode, $speices);
    my $shortNameKey = "SquashSplice $elandGenome/$elandGenomeMask $genomeAnnotation $readLength $chromNameSource $squashGenomeParams $spliceSitesParams ";
    my $spliceFlank = $readLength - 1;
    return _makeShortName($shortNameKey, $elandGenome, $elandGenomeMask, %$knownShortnames, '', "Splice$spliceFlank");
}

sub _makeContamShortName($$$$\%)
{
    my($lane, $barcode, $sampleSheet, $config, $knownShortnames) = @_;
    my($project, $sample, $speices) = ($sampleSheet->projectId($lane, $barcode), 
                                       $sampleSheet->sampleId($lane, $barcode), 
                                       $sampleSheet->species($lane, $barcode));

    my $elandGenome = $config->selectValue('ELAND_RNA_GENOME_CONTAM', $project, $sample, $lane, $barcode, $speices);
    my $elandGenomeMask = $config->selectValue('ELAND_RNA_GENOME_CONTAM_MASK', $project, $sample, $lane, $barcode, $speices);
    my $squashGenomeParams = $config->selectValue('SQUASH_GENOME_PARAMS', $project, $sample, $lane, $barcode, $speices);
    my $shortNameKey = "SquashContam $elandGenome/$elandGenomeMask $squashGenomeParams";
    return _makeShortName($shortNameKey, $elandGenome, $elandGenomeMask, %$knownShortnames, 'contam');
}

sub _makeGenomeShortName($$$$\%)
{
    my($lane, $barcode, $sampleSheet, $config, $knownShortnames) = @_;
    my($project, $sample, $speices) = ($sampleSheet->projectId($lane, $barcode), 
                                       $sampleSheet->sampleId($lane, $barcode), 
                                       $sampleSheet->species($lane, $barcode));

    my $elandGenome = $config->selectValue('ELAND_GENOME', $project, $sample, $lane, $barcode, $speices);
    my $elandGenomeMask = $config->selectValue('ELAND_GENOME_MASK', $project, $sample, $lane, $barcode, $speices);
    my $squashGenomeParams = $config->selectValue('SQUASH_GENOME_PARAMS', $project, $sample, $lane, $barcode, $speices);
    my $shortNameKey = "SquashGenome $elandGenome/$elandGenomeMask $squashGenomeParams";
    return _makeShortName($shortNameKey, $elandGenome, $elandGenomeMask, %$knownShortnames, 'genome');
}


=pod

=head2 configureAnalysis($sampleSheet, $output)

Load the sample sheet and generate the configuration for analysis of the specified dataset.

  use Casava::Alignment;
  Casava::Alignment::configureAnalysis($sampleSheet, $output);

=over 4

=item *

Parameters:

  $configFile - full path to the input XML config file for the analysis folder
  $output     - reference to the output FILE HANDLE (defaults to \*STDOUT)

=item *

Returns:

  Nothing

=item *

Exceptions

  Incorrect number of arguments
  Config file does not exist
  XML config file can't be read (no read access or invalid XML)
  Missing elements in the XML configuration file
  Failing to write to the indicated output

=back

=cut

sub configureAnalysis($$;$)
{
    my ($sampleSheetFile, $configFile, $output) = @_;
    errorExit("ERROR: $sampleSheetFile: file does not exist") unless -e $sampleSheetFile;
    errorExit("ERROR: $configFile: file does not exist") unless -e $configFile;
    use Casava::Demultiplex::SampleSheet;
    my $sampleSheet = Casava::Demultiplex::SampleSheet::Csv->new();
    open my $sampleSheetHandle, '<', "$sampleSheetFile" or errorExit("ERROR: Failed to open $sampleSheetFile: $!");
    $sampleSheet->load($sampleSheetHandle, 0) or errorExit("Couldn't load $sampleSheetFile");
    close $sampleSheetHandle or errorExit("ERROR: Failed to close $sampleSheetFile: $!");
    # Reverse the order of the elements
    my $reversed = {};
    foreach my $lane ($sampleSheet->laneNumbers)
    {
        foreach my $barcode ($sampleSheet->barcodes($lane))
        {
            my $sample = $sampleSheet->sampleId($lane, $barcode);
            # Empty sample names in SampleSheet.csv are an instruction to not write out the data for thos samples
            if (length($sample))
            {
                my $project = $sampleSheet->projectId($lane, $barcode);
                $reversed->{$project} = {} unless exists $reversed->{$project};
                $reversed->{$project}->{$sample} = {} unless exists $reversed->{$project}->{$sample};
                $reversed->{$project}->{$sample}->{$lane} = [] unless exists $reversed->{$project}->{$sample}->{$lane};
                push @{ $reversed->{$project}->{$sample}->{$lane} }, $barcode;
            }
        }
    }
    $output = \*STDOUT unless $output;
    my $config = Casava::Alignment::Config->new;
    $config->readXml($configFile);
    my %knownGenomeShortnames;
    my @configuration;
    push @configuration, 'FLOWCELL:=' . $config->getVariable(name => 'FLOWCELL') . "\n";
    my @projectList;
    foreach my $project (sort keys %$reversed)
    {
        my $projectDir = ('Undetermined_indices' eq $project) ? $project : "Project_$project";
        my $projectRef = $reversed->{$project};
        my @sampleList;
        my @allProjectTiles=();
        my @allProjectPaddedTiles=();
        foreach my $sample (sort keys %$projectRef)
        {
            my $sampleRef = $projectRef->{$sample};
            my @laneList;
            foreach my $lane (sort keys %$sampleRef)
            {
                my $tilesString = $config->getVariable(lane=>$lane, name=>'Tiles');
                next unless defined $tilesString;
                my @tiles = split /\s+/, $tilesString;
                next unless @tiles;
                my $laneRef = $sampleRef->{$lane};
                my @barcodes = map {$_ ? $_ : 'empty'} @$laneRef;
                push @configuration, "${project}_${sample}_${lane}_BARCODES:=" . join(" ", sort @barcodes);
                foreach my $barcode (sort @barcodes)
                {
                    my $species = $sampleSheet->species($lane, $barcode);
                    my $analysis = $config->selectValue('ANALYSIS', $project, $sample, $lane, $barcode, $species);
                    my @reads = split(' ', $config->selectValue('READS', $project, $sample, $lane, $barcode, $species));
                    push @configuration, "${project}_${sample}_${lane}_${barcode}_REFERENCE:=" . $sampleSheet->species($lane, $barcode);
                    push @configuration, "${project}_${sample}_${lane}_${barcode}_GENOME_SHORTNAME:=" . _makeGenomeShortName($lane, $barcode, $sampleSheet, $config, %knownGenomeShortnames);
                    if ('eland_rna' eq $analysis)
                    {
                        push @configuration, "${project}_${sample}_${lane}_${barcode}_CONTAM_SHORTNAME:=" . _makeContamShortName($lane, $barcode, $sampleSheet, $config, %knownGenomeShortnames);
                        foreach my $read (@reads)
                        {
                            push @configuration, "${project}_${sample}_${lane}_${barcode}_SPLICE_SHORTNAME_R${read}:=" . _makeSpliceShortName($lane, $barcode, $read, $sampleSheet, $config, %knownGenomeShortnames);
                        }
                    }
                }
                push @configuration, "${project}_${sample}_${lane}_TILES:=" . $tilesString;
                push @configuration, "${project}_${sample}_${lane}_PADDED_TILES:=" . 
                    join(' ', map ({sprintf "%04d", $_} @tiles));
                push @allProjectTiles, @tiles;
                push @allProjectPaddedTiles, map ({sprintf "%04d", $_} @tiles);
                push @laneList, $lane;
            }
            next unless @laneList;
            push @configuration, "${project}_${sample}_LANES:=" . join(" ", sort @laneList) . "\n";
            push @configuration, "${project}_${sample}_OUT_DIR:=" . File::Spec->catdir($config->getVariable(name => 'OUT_DIR'), $projectDir, "Sample_$sample");
            push @configuration, "";
            push @sampleList, $sample;
        }
        next unless @sampleList;
        push @configuration, "${project}_SAMPLES:=" . join(" ", sort @sampleList) . "\n";
        push @configuration, "${project}_TILES:=" . join(" ", keys %{{map{$_ => 1} @allProjectTiles}}) . "\n";
        push @configuration, "${project}_PADDED_TILES:=" . join(" ", keys %{{map{$_ => 1} @allProjectPaddedTiles}}) . "\n";
        push @configuration, "${project}_OUT_DIR := " . File::Spec->catdir($config->getVariable(name => 'OUT_DIR'), $projectDir). "\n";
        # absolute paths in makefiles make it difficult to debug targets
        push @configuration, "${project}_SUMMARY_DIR :=" . File::Spec->catdir($projectDir, 'Summary_Stats_$(FLOWCELL)'). "\n";
        push @projectList, $project;
    }
    push @configuration, "PROJECTS:=" . join(" ", sort @projectList) . "\n";
    push @configuration, '', '# Global variables from config file', '';
    push @configuration, 'EXPT_DIR := ' . $config->getVariable(name => 'EXPT_DIR');
    push @configuration, 'EMAIL_LIST := ' . $config->getVariable(name => 'EMAIL_LIST');
    push @configuration, 'EMAIL_SERVER := ' . $config->getVariable(name => 'EMAIL_SERVER');
    push @configuration, 'EMAIL_DOMAIN := ' . $config->getVariable(name => 'EMAIL_DOMAIN');
    push @configuration, 'WEB_DIR_ROOT := ' . $config->getVariable(name => 'WEB_DIR_ROOT');
    push @configuration, 'NUM_LEADING_DIRS_TO_STRIP := ' . $config->getVariable(name => 'NUM_LEADING_DIRS_TO_STRIP');
    print $output join("\n", @configuration, '') or errorExit("ERROR: Failed to write the configuration of the analysis: $!");
}

=pod

=head2 configureDataset($configFile, $project, $sample, $lane, $barcode, $output)

Load the 'config.xml' and generate the configuration for analysis of the specified dataset.

  use Casava::Alignment;
  Casava::Alignment::configureDataset($configFile, $output);

=over 4

=item *

Parameters:

  $configFile - full path to the input XML config file for the analysis folder
  $project    - name of the project
  $sample     - name of the sample
  $lane       - lane number
  $barcode    - barcode sequence
  $output     - reference to the output FILE HANDLE (defaults to \*STDOUT)

=item *

Returns:

  Nothing

=item *

Exceptions

  Incorrect number of arguments
  Config file does not exist
  XML config file can't be read (no read access or invalid XML)
  Missing elements in the XML configuration file
  Failing to write to the indicated output

=back

=cut

sub configureDataset($$$$$$;$)
{
    my ( $configFile, $project, $sample, $lane, $barcode, $reference, $output) = @_;
    errorExit("ERROR: $configFile: file does not exist") unless -e $configFile;
    #my $configRef = XMLin($configFile, SuppressEmpty => 1) or croak "ERROR: couldn't load $configFile: $!";
    #errorExit("ERROR: no 'Defaults' parameters in config file $configFile") unless exists $configRef->{'Defaults'};
    #my $defaults = $configRef->{'Defaults'};
    my $config = Casava::Alignment::Config->new;
    $config->readXml($configFile);
    $config->validateDataset($project, $sample, $lane, $barcode, $reference);
    my @content = ("##\n## Auto-generated file. Do not edit.\n##\n");
    push @content, 'dataset := $(project)_$(sample)_$(lane)_$(barcode)';
    push @content, '$(dataset)_PROJECT := $(project)';
    push @content, '$(dataset)_SAMPLE := $(sample)';
    push @content, '$(dataset)_LANE := $(lane)';
    push @content, '$(dataset)_TILES := $($(project)_$(sample)_$(lane)_TILES)';
    push @content, '$(dataset)_BARCODE := $(barcode)';
    if ('Undetermined_indices' eq $project)
    {
        push @content, '$(dataset)_SOURCE := $(project)/Sample_$(sample)';
        push @content, '$(dataset)_DESTINATION := $(project)/Sample_$(sample)';
    }
    else
    {
        push @content, '$(dataset)_SOURCE := Project_$(project)/Sample_$(sample)';
        push @content, '$(dataset)_DESTINATION := Project_$(project)/Sample_$(sample)';
    }
    push @content, '$(dataset)_TEMP_DIR := $(TEMP_DIR)/$($(dataset)_DESTINATION)';
    push @content, '$(dataset)_PLOTS_DIR := $($(project)_SUMMARY_DIR)/Plots/Sample_$(sample)';
    push @content, '$(dataset)_STATS_DIR := $($(project)_SUMMARY_DIR)/Stats/Sample_$(sample)';
    push @content, '$(dataset)_FILE_PREFIX := $(sample)_$(barcode)_L00$(lane)';
    push @content, '$(dataset)_SOURCE_PREFIX := $($(dataset)_SOURCE)/$($(dataset)_FILE_PREFIX)';
    push @content, '$(dataset)_DESTINATION_PREFIX := $($(dataset)_DESTINATION)/$($(dataset)_FILE_PREFIX)';
    push @content, '$(dataset)_TEMP_PREFIX := $($(dataset)_TEMP_DIR)/$($(dataset)_FILE_PREFIX)';
    push @content, '$($(dataset)_DESTINATION_PREFIX)%: project := $(project)';
    push @content, '$($(dataset)_DESTINATION_PREFIX)%: sample := $(sample)';
    push @content, '$($(dataset)_DESTINATION_PREFIX)%: lane := $(lane)';
    push @content, '$($(dataset)_DESTINATION_PREFIX)%: barcode := $(barcode)';
    push @content, '$($(dataset)_DESTINATION_PREFIX)%: dataset := $(dataset)';
    push @content, '$($(dataset)_DESTINATION_PREFIX)_R1%: read := 1';
    push @content, '$($(dataset)_DESTINATION_PREFIX)_R2%: read := 2';
    push @content, '$($(dataset)_TEMP_PREFIX)%: project := $(project)';
    push @content, '$($(dataset)_TEMP_PREFIX)%: sample := $(sample)';
    push @content, '$($(dataset)_TEMP_PREFIX)%: lane := $(lane)';
    push @content, '$($(dataset)_TEMP_PREFIX)%: barcode := $(barcode)';
    push @content, '$($(dataset)_TEMP_PREFIX)%: dataset := $(dataset)';
    push @content, '$($(dataset)_TEMP_PREFIX)_R1%: read := 1';
    push @content, '$($(dataset)_TEMP_PREFIX)_R2%: read := 2';
    push @content, '$(dataset)_TARGET := $($(dataset)_DESTINATION_PREFIX)_finished.txt';
    push @content, '$(dataset)_POSTRUN := $($(dataset)_DESTINATION_PREFIX)_postrun.done';
    push @content, '$($(dataset)_TARGET): config.xml';
    push @content, 'DATASET_TARGETS += $($(dataset)_POSTRUN)';
    my %dictionary = $config->getDefaultVariables;
    my %isMacro = (DATASET_POST_RUN_COMMAND=>1);
    foreach my $variable (sort keys %dictionary)
    {
        my $assignment = ' := ';
        $assignment = ' = ' if $isMacro{$variable};
        my $value = $config->selectValue($variable, $project, $sample, $lane, $barcode, $reference);
        my $line = '$(dataset)_' . $variable . $assignment . $value;
        push @content, '', 'ifneq (,$($(dataset)_TILES))' if 'ANALYSIS' eq $variable;
        push @content, $line;
        push @content, 'else', '$(dataset)_ANALYSIS := none', 'endif', '' if 'ANALYSIS' eq $variable;
    }
    push @content,
    '',
    '# Include the analysis-specific makefile for the dataset',
    '#ifneq (,$($(dataset)_TILES))',
    'include $(MAKEFILES_DIR)/$($(dataset)_ANALYSIS).mk',
    '#endif',
    '';
    $output = \*STDOUT unless $output;
    print $output join("\n", @content, "\n") or errorExit("ERROR: Failed to print the configuration: $!");
}

sub _createDirectoryStructure
{
    my ($outputDirectory) = @_;
    if (!-d $outputDirectory)
    {
        mkdir $outputDirectory or errorExit("ERROR: Failed to create output directory $outputDirectory: $!");
    }
}

sub _copyConfigFile
{
    my ($configFile, $outputDirectory) = @_;
    errorExit("ERROR: missing config file") unless $configFile;
    errorExit("ERROR: missing output directory") unless $outputDirectory;
    my $outputFile = File::Spec->catfile($outputDirectory, 'config.txt');
    errorExit("ERROR: output config.txt already exists: $outputFile") if -e $outputFile;
    # copy the content to support any type of file (on-disk, in-memory, etc.)
    open my $input, '<', $configFile or errorExit("ERROR: Failed to open input config file $configFile: $!");
    open my $output, '>', $outputFile or errorExit("ERROR: Failed to open output config file $outputFile: $!");
    while (<$input>)
    {
        print $output $_ or errorExit("Failed to write into the output config file $outputFile: $!");
    }
    close $output or errorExit("ERROR: Failed to close output config file $outputFile: $!");
    close $input or errorExit("ERROR: Failed to close input config file $configFile: $!");
}

sub _generateSampleSheet
{
    my ($inputDirectory, $sampleSheetFile, $outputDirectory) = @_;
    my $demultiplexConfigFile = File::Spec->catfile($inputDirectory, 'DemultiplexConfig.xml');
    my $outputFile = File::Spec->catfile($outputDirectory, 'SampleSheet.csv');
    open my $output, '>', $outputFile or errorExit("ERROR: Failed to open output sample sheet $outputFile: $!");
    if ($sampleSheetFile || !-e $demultiplexConfigFile)
    {
        $sampleSheetFile = File::Spec->catfile($inputDirectory, 'SampleSheet.csv') unless $sampleSheetFile;
        open my $input, '<', $sampleSheetFile or errorExit("ERROR: Failed to open input sample sheet $sampleSheetFile: $!");
        while (<$input>)
        {
            print $output $_ or errorExit("ERROR: Failed to write into output sample sheet $outputFile: $!");
        }
        close $input or errorExit("ERROR: Failed to close input sample sheet $sampleSheetFile: $!");
    }
    elsif (-e $demultiplexConfigFile)
    {
        my $dc = Casava::Demultiplex::DemultiplexConfig->new();
        $dc->load($demultiplexConfigFile);
        my $ss = Casava::Demultiplex::SampleSheet::create($outputFile);
        $ss->clone($dc);
        $ss->save($output);
    }
    close $output or errorExit("ERROR: Failed to close output sample sheet $outputFile: $!");
}

sub _printSupportTxt($$)
{
    my ($outputDirectory, $commandLine) = @_;
    open my $support, '>', File::Spec->catfile($outputDirectory, 'support.txt') or errorExit("ERROR: Failed to open 'support.txt': $!");

    my $data = Data::Dumper->new([{
        'OS' => "$^O", 
        'PID' => "$$",
        'PERL_VERSION' => "$^V", 
        'PERL_EXECUTABLE' => "$^X",
    }, 'bcl2fastq-1.8.4', $0, $commandLine, \@INC, [sort(values %INC)], \%ENV],
    [qw(_System _ID-string _Program _Command-line _Locations _Modules _Environment)]);
    $data->Pair(' : '); $data->Indent(1);
    my $supportContent = join "\n", $data->Dump;
    print $support <<"SUPPORT_END" or errorExit("ERROR: Failed to write into 'support.txt': $!");
# This file contains information about your system that is useful for
# diagnosing a problem you may have. For technical assistance, please 
# contact the Illumina Customer Support team <techsupport\@illumina.com>\n"&&
# and send them this file.

$supportContent

SUPPORT_END

    close $support or errorExit("ERROR: Failed to close 'support.txt': $!");
}

=pod

=head2 initialize(%optionsMap)

Load the default GERALD.xml, load the config.txt, load the command line options,
write the config.txt, write the config.xml and copy the default Makefile.

  use Casava::Alignment;
  Casava::Alignment::initialize('config-file' => $configFile,
                                'make' => $make,
                                'experiment-dir' => $experimentDir,
                                'output-dir' => $outputDir);

=over 4

=item *

Parameters:

  %optionsMap - HASH MAP to the command line parameters

=item *

Returns:

    Nothing

=back

=cut

sub initialize
{
    my ($config, $userSampleSheetFile, $configFile, $commandLine, $skipVariableMetadata) = @_;
    my $inputDirectory = $config->getVariable(name => 'EXPT_DIR');
    my $outputDirectory = $config->getVariable(name => 'OUT_DIR');

    _createDirectoryStructure($outputDirectory);
    _printSupportTxt($outputDirectory, $commandLine);

    _copyConfigFile($configFile, $outputDirectory);

    # store the software version history
    $config->software({Name => (File::Spec->splitpath($0))[2], Version => 'bcl2fastq-1.8.4'});
    $config->cmdAndArgs($commandLine) unless $skipVariableMetadata;
    my $demultiplexConfigFile = File::Spec->catfile($inputDirectory, 'DemultiplexConfig.xml');
    if (-e $demultiplexConfigFile)
    {
        my $dc = Casava::Demultiplex::DemultiplexConfig->new();
        $dc->load($demultiplexConfigFile);
        $config->demuxSoftware($dc->software());
    }

    _generateSampleSheet($inputDirectory, $userSampleSheetFile, $outputDirectory);

    $config->writeXml(File::Spec->catfile($outputDirectory, 'config.xml'));

    my $inputMakefile = File::Spec->catfile( '/usr/local/share/bcl2fastq-1.8.4', "makefiles", 'Makefile');
    my $outputMakefile = File::Spec->catfile( $outputDirectory, 'Makefile');

    use File::Copy;
    copy($inputMakefile, $outputMakefile) or die "Failed to copy $inputMakefile to $outputMakefile: $!";
}

sub selfTest
{
    my ($geraldDir) = @_;
    use Cwd;
    my $currDir = getcwd();
    chdir $geraldDir;
    my $selfTestCommand = "make self_test";
    logInfo("Running self tests: '$selfTestCommand'", 0);
    my $ENV_REF = \%ENV;
    $ENV_REF->{'PATH'} = $1 if $ENV_REF->{'PATH'} =~ /^(.*)$/;    # untaint PATH variable
    my $selfTestOutput = `$selfTestCommand 2>&1`;
    my $selfTestResult = $?;
    my $exitValue = ($selfTestResult >> 8);
    if ( $exitValue != 0 )
    {
        my $signal = ($selfTestResult & 127);
        my $colorizeSh = File::Spec->catfile( '/usr/local/libexec/bcl2fastq-1.8.4', 'colorize.sh');
        open (COLOREDOUTPUT, "|($colorizeSh 'ERROR.*\$') 1>&2 ") or
            errorExit "ERROR: Could not execute $colorizeSh: $!";
        print COLOREDOUTPUT "\noutput of '$selfTestCommand':\n\n$selfTestOutput\n\n" .
             "Self test command exited with error $exitValue (signal $signal)\n\n" .
             "Investigate and fix errors, then run the command again\n\n";
        close COLOREDOUTPUT;
        chdir $currDir;
        errorExit("Self test command exited with error $exitValue (signal $signal)");
    }
    chdir $currDir;
    logInfo("Running self tests on $geraldDir completed with no problems\n", 0);
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

strict, warnings, Carp, GetOpt::Long, Pod::Usage, 

=item External perl modules

XML::Simple

=item Casava perl modules

Casava::Alignment

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Come Raczy

=cut
