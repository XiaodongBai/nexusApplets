# PROJECT: CASAVA
# MODULE:  $RCSfile: IOLib.pm,v $
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#

=pod

=head1 NAME

Casava::Common::IOlib.pm - Perl utility library for accessing system resources.

=head1 SYNOPSIS

# include what functions you need... 
use Casava::Common::IOlib qw();  

=head1 DESCRIPTION

Exports: 
    bufferedPrint($;\$;$);
    createDirs ($;\@);
    createFilesInDirs ($;$;\%;\@);
    executeCmd ($;$);
    genId($;$;\$;\$);
    getProcessorsCount ();
    writeParameters(\%;$;$;$,$);
    testIfSupports($);
    readParameters(\%$$);
    
# Global variable

=head1 AUTHORSHIP

Copyright (c) 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

Created by Lukasz Szajkowski <lszajkowski@illumina.com>

=head1 SOURCE CODE

The most current release of the Perl source code 
for this module is available through CVS at genie01 a.k.a. 10.44.0.81 
cvs co BullFrog

=cut

package Casava::Common::IOLib;

#
# Place functions/variables you want to *export*, ie be visible
# from the caller package into @EXPORT_OK
#
BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT =
      qw($VERSION &executeCmd &executeCmdEx &executeCmdOut2Str &genId &readParameters);
    @EXPORT_OK = qw(&createFilesInDirs &createDirs createDir &string2file
      &getProcessorsCount &writeParameters &readXML &writeXML &wcl 
      &bufferedPrint &testIfSupports &table2file &table2fileEx &file2table &file2tableEx &formatValueWithDeviation);
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);
use XML::Simple;
use Carp;
use List::Util qw(max);
use IO::File;

use Casava::Common::Log;

sub bufferedPrint($;\$;$);
sub createDir ($);
sub createDirs ($\@);
sub createFilesInDirs ($;$;\%;\@);
sub executeCmd ($;$);
sub executeCmdEx ($$;$);
sub executeCmdOut2Str ($$);
sub genId($;$;\$;\$);
sub writeParameters(\%;$;$;$;$);
sub testIfSupports($);
sub readXML(\%;$;$);
sub writeXML(\%;$);
sub readParameters(\%$$);
sub getProcessorsCount ();
sub file2table(\@;$);
sub file2tableEx(\@;\%;$);
sub table2file(\@;$);
sub table2fileEx(\@;$;$);
sub string2file($;$);
sub wcl($);
sub formatValueWithDeviation($;);


my %buffer            = ();
my $bufferSize        = 10000;
my $currentBufferSize = 0;

# Global variable

=pod

=head1 The procedure parses the cvs string.

=over 4

=item genId($cvs_id, $cvs_version, $idRef, $verRef)

The procedure parses the cvs string into id and version strings.

Parameters:
    $cvs_id         - cvs id string 
    $cvs_version    - cvs version string
    $idRef          - out id variable
    $verRef         - out version variable

Returns:
    nothing
    
=back

=cut

sub genId($;$;\$;\$) {
    croak "ERROR: genId\n" unless ( @_ == 4 );
    my ( $cvs_id, $cvs_version, $idRef, $verRef ) = @_;
    ${$verRef} = $cvs_version;
    ${$verRef} =~ s/\$//g;
    ${$verRef} =~ s/Revision//;
    ${$idRef} = $cvs_id;
    ${$idRef} =~ s/\$//g;
}

=pod

=head1 The procedure executes the command.

=over 4

=item executeCmd($command, $verbose)

The procedure executes the command but exits with error code when 
command failed.

Parameters:
    command - command to be executed
    verbose - verbose level (optional)
Returns
    Nothing 
=back

=cut
sub executeCmd ($;$) {
    my ( $command, $verbose ) = @_;
    executeCmdEx($command, 1, $verbose);
}

=pod

=head1 The procedure executes the command via system() call.

=over 4

=item executeCmdEx($command, $dieOnFailure, $verbose)

The procedure executes the command. If $dieOnFailure is set, exits the process
with error code if command fails.

Parameters:
    command - command to be executed
    dieOnFailure - whether to terminate process in case the command fails
    verbose - verbose level (optional)
Returns
    process exit code if dieOnFailure is not set. 0 if successful
=back

=cut
sub executeCmdEx ($$;$) {
    my ( $command, $dieOnFailure, $verbose ) = @_;
    logInfo( $command, $verbose );
    system($command);
    my $ret = $?;
    my $exit_value  = $ret >> 8;
    my $signal_num  = $ret & 127;
    my $dumped_core = $ret & 128;
    my $message = ( $signal_num > 0 ) 
        ? "Process killed by SIG($signal_num): $command - See errors above\n"
        : ( $dumped_core > 0 )
            ? "$command dumped core: $command - See errors above\n"
            : ( $exit_value != 0 )
                ? "In $command - See errors above\n"
                : undef;

    if (defined $message){
        errorExit ("ERROR:$message") unless !$dieOnFailure;
        logWarning $message;
    }
    return $ret;
}

=pod

=head1 The procedure executes the command.

=over 4

=item executeCmd($command, $verbose)

The procedure executes the command but exits with error code when 
command failed.

Parameters:
    command - command to be executed
    verbose - verbose level (optional)

=back

=cut

sub executeCmdOut2Str ($$) {
    my ( $command, $verbose ) = @_;
    printLog( $command . "\n", $verbose );
    my $count = 0;
    my $returnCode = open( CMD, "$command 2>/dev/null |" );
    if ( !$returnCode ) {
        die "ERROR: Couldn't run $command $! \n";
    }
    my $resultContent            = '';
    my $resultContentNoEndOfLine = '';
    while ( defined( my $line = <CMD> ) ) {
        $resultContent .= $line;
        chomp($line);
        $resultContentNoEndOfLine .= $line;
    }
    my $exit_value  = $? >> 8;
    my $signal_num  = $? & 127;
    my $dumped_core = $? & 128;
    close(CMD);
    if ( $signal_num > 0 ) {
        die
"ERROR: executeCmdOut2Str() Process killed by SIG($signal_num): $command : $resultContentNoEndOfLine\n";
    }    # if
    if ( $dumped_core > 0 ) {
        die
"ERROR: executeCmdOut2Str() $command dumped core: $command : $resultContentNoEndOfLine\n";
    }    # if
    if ( $exit_value != 0 ) {
        die
"ERROR: executeCmdOut2Str() EXIT_CODE($exit_value): $command : $resultContentNoEndOfLine $!\n";
    }    # if
    printLog( $resultContent, $verbose );
    return $resultContent;
}

=pod

=head1 The procedure creates files 

=over 4

=item createFilesInDirs($fileNamePref, $dirLocation, $filesRef, $dirsRef)

The procedure creates files named $fileNamePref.txt 
    in the @{$dirsRef} directories
Parameters:
    $fileNamePref   - file name pref used to create each file (eg.: wrong)
    $dirLocation    - temp directory location
    $filesRef       - Hash map where handles to files will be stored
    $dirsRef        - list of directories

Returns:
    nothing
=back

=cut

sub createFilesInDirs ($;$;\%;\@) {
    croak "ERROR: createFilesInDirs wrong parameters \n"
      unless ( @_ == 4 );
    my ( $fileNamePref, $dirLocation, $filesRef, $dirsRef ) = @_;
    foreach my $dir ( @{$dirsRef} ) {
        unless ( defined $filesRef->{$dir} ) {
            my $outChrDir = File::Spec->catdir( $dirLocation, $dir );
            my $file =
              File::Spec->catfile( $outChrDir, $fileNamePref . '.txt' );
            $filesRef->{$dir} = IO::File->new(">>$file")
              or die "ERROR: createFilesInDirs() Couldn't open " . "$file $!\n";
        }    # unless
    }    # foreach
}

=pod

=head1 The procedure creates directories 

=over 4

=item createDirs($dir)

The procedure creates a directory

Parameters:
    $dirLocation    - dir to create
    $dirsRef        - list of directories

Returns:
    nothing
=back

=cut

sub createDir ($) {
    my ( $dir ) = @_;
    unless ( -d $dir ) {
        mkdir $dir, 0775
          or errorExit "ERROR: createDir() can't make $dir $!";
    }
}


=pod

=head1 The procedure creates directories 

=over 4

=item createDirs($dirLocation, $dirsRef)

The procedure creates directories from @{$dirsRef} list

Parameters:
    $dirLocation    - temp directory location
    $dirsRef        - list of directories

Returns:
    nothing
=back

=cut

sub createDirs ($\@) {
    my ( $dirLocation, $dirsRef ) = @_;
    foreach my $dir ( @{$dirsRef} ) {
        my $outChrDir = File::Spec->catdir( $dirLocation, $dir );
        createDir ( $outChrDir );
    }
}

=pod

=head1 Procedure to read the number of processors 

=over 4

=item getProcessorsCount()

The procedure to count the number of processors on current machine.

Returns:
    number of processors
=back

=cut

sub getProcessorsCount () {
    croak "ERROR: getProcessorsCount wrong parameters \n"
      unless ( @_ == 0 );
    my $command = "cat /proc/cpuinfo | grep processor | wc -l 2> /dev/null";
    my $count   = 0;
    open( SGE, "$command |" )
      || die
"ERROR: getProcessorsCount() Couldn't run cat /proc/cpuinfo | grep processor | wc $! \n";
    while (<SGE>) {
        chomp($_);
        $count = $_;
    }    # while
    close(SGE);
    return $count;
}

=pod

=head1 The procedure reads HASH MAP from the file.

=over 4

=item readParameters( $parametersRef, $sectionName, $fileName)

The procedure reads HASH MAP with parameters from the parameters file.

Parameters:
    $parametersRef  - HASH MAP with parameters
    $sectionName    - parameters section name in the file
    $fileName       - file name

Returns:
    nothing
=back

=cut

sub readParameters(\%$$) {
    my ( $parametersRef, $sectionName, $fileName ) = @_;
    my $sectionNameTmp = "";
    my $currentParfRef = "";
    open( FILE, "<$fileName" )
      || die
      "ERROR: readParameters() Couldn't open file handle for $fileName $!\n";
    while (<FILE>) {
        my $line = $_;
        chomp $line;
        my $comment = "";
        my $key     = "";
        my $value   = "";
        if ( $line =~ /(#.+)$/ ) {
            $comment = $1;
        }    # if
        elsif ( $line =~ /^\[(\S+)\]$/ ) {
            $sectionNameTmp = $1;
            if ( $sectionNameTmp eq $sectionName ) {
                $currentParfRef = $parametersRef;
            }    # if
            elsif ( $currentParfRef ne "" ) {
                $currentParfRef = "";
            }
        }    # if
        elsif ( $line =~ /^(\S+)\s*(.*)$/ && $currentParfRef ne "" ) {
            $key   = $1;
            $value = $2;
            ${$currentParfRef}{$key} = $value;
        }    # elsif
    }    # while
    close FILE;
}

=pod

=head1 The procedure writes/appends HASH MAP the file.

=over 4

=item writeParameters( $parametersRef, $sectionName, $fileName, $sortType, $widthPar )

The procedure writes/appends HASH MAP with parameters to a parameters file

Parameters:
    $parametersRef  - HASH MAP with parameters
    $sectionName    - parameters section name in the file
    $fileName       - file name
    $sortType       - sort type ( 0 - key sort lexical, 1 - key sort num
                                  2 - value sort lexical, 3 -value sort num)
    $widthPar          - optional width (default 20)
Returns:
    0 if success; -1 if error
=back

=cut

sub writeParameters(\%;$;$;$;$) {
    croak "ERROR: writeParameters\n" unless ( @_ == 4 || @_ == 5 );

    #print "writeParameters(". join (',', @_).")\n";
    my ( $parametersRef, $sectionName, $fileName, $sortType, $widthPar ) = @_;
    my $sectionNameTmp = "";
    my $overwrite      = 0;
    my $found          = 0;
    my $content        = "";
    my $width          = 32;
    if ( defined $widthPar ) {
        $width = $widthPar;
    }
    my @keys = ();
    if ( $sortType == 0 ) {
        @keys = sort keys %{$parametersRef};
    }    # if
    elsif ( $sortType == 1 ) {
        @keys =
          sort { ${$parametersRef}{$a} <=> ${$parametersRef}{$b} }
          keys %{$parametersRef};
    }    # elsif
    elsif ( $sortType == 2 ) {
        @keys =
          sort { ${$parametersRef}{$a} cmp ${$parametersRef}{$b} }
          keys %{$parametersRef};
    }    # elsif
    elsif ( $sortType == 3 ) {
        @keys = sort { ${$parametersRef}{$a} <=> ${$parametersRef}{$b} }
          keys %{$parametersRef};
    }    # else
    my $contentBody .= "[$sectionName]\n";
    foreach my $key (@keys) {
        my $value = "";
        if ( defined( ${$parametersRef}{$key} ) ) {
            $value = ${$parametersRef}{$key};
        }    # if                           # if
        $contentBody .= ( sprintf "%-" . max($width, length($key) +1) . "s", $key ) . $value . "\n";
    }    # foreach
    $contentBody .= "\n";
    if ( -e $fileName ) {    # Merge with file
        open( FILE, "<$fileName" )
          || die
"ERROR: writeParameters() Couldn't open file handle for $fileName $!\n";
        while (<FILE>) {
            my $line = $_;
            chomp $line;

            #        print $line . " ";
            my $comment = "";
            my $key     = "";
            my $value   = "";
            if ( $line =~ /^\[(\S+)\]$/ ) {
                $sectionNameTmp = $1;
                if ( $sectionNameTmp eq $sectionName ) {
                    $overwrite = 1;
                    $found     = 1;
                    $content .= $contentBody . "\n";
                }    # if
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
          . " parameters for $sectionName\n\n";
        $content .= $contentBody . "\n";
    }    # if
    my $file = IO::File->new( ">${fileName}.tmp" ) || errorExit "ERROR: Failed to create ${fileName}.tmp. $!";

    print $file $content;
    close $file;
    rename ("${fileName}.tmp", $fileName) || errorExit "ERROR: Cannot rename ${fileName}.tmp into $fileName. $!";
    return 0;
}

=pod

=head1 Buffered write to the file

=over 4

=item bufferedPrint($fileHandle, $contentRef, $isFlush)

The procedure writes the content to the file as bufferSize chunks of text;


Parameters:
    $fileHandle    - file to store the content 
    $contentRef    - reference to the content
    $content      - flush the buffer

Returns:
    nothing
        
=back

=cut

sub bufferedPrint($;\$;$) {
    croak "ERROR: bufferedPrint " . join( ',', @_ ) unless ( @_ == 3 );
    my ( $fileHandle, $contentRef, $isFlush ) = @_;
    if ( $fileHandle > -1 ) {
        print $fileHandle ${$contentRef};
        if ( $isFlush ) {
            $fileHandle->flush();
        }
    }
}

=pod

=head1 The procedure tests if current machine supports a command.

=over 4

=item testIfSupports($command)

The procedure tests if current machine supports given command(executes it ).

Parameters:
    $command - some command

Returns:
    1 if machine supports the command
=back

=cut

sub testIfSupports($) {
    croak "ERROR: testIfSupports\n" unless ( @_ == 1 );
    my ($command) = @_;
    my $res = 1;
    system( $command . " 1>/dev/null 2> /dev/null" );
    my $exit_value  = $? >> 8;
    my $signal_num  = $? & 127;
    my $dumped_core = $? & 128;
    if ( $signal_num > 0 ) {
        $res = 0;
    }    # if
    if ( $dumped_core > 0 ) {
        $res = 0;
    }    # if
    if ( $exit_value != 0 ) {
        $res = 0;
    }    # if
    return $res;
}

=pod

=item string2file($content, $fileName)

string2file procudure saves a string to a file.

B<Parameters:>

    $content            - content to be saved
    $fileName           - filename where content will be saved
    
B<Returns:> 
    nothing
    
=cut

sub string2file($;$) {
    croak "$0::ERROR: string2file " . join( ',', @_ ) unless ( @_ == 2 );
    my ( $content, $fileName ) = @_;
    my $file = IO::File->new( ">" . $fileName )
      || die
"ERROR: string2file() Couldn't create/open file handle for $fileName $!\n";
    print $file $content;
    close $file;
}

=pod

=head1 The procedure stores table to a file

=over 4

=item table2file($tableRef, $fileName)

The procedure stores table to a file.

Parameters:
    $tableRef     - ARRAYref of ARRAYS reference to table  
    $fileName     - file name where table should be stored 

Returns:
    table as text
    
=back

=cut

sub table2file(\@;$) {
    croak "$0:ERROR: table2file\n" unless ( @_ == 2 );
    my ( $tableRef, $fileName ) = @_;

    return table2fileEx(@$tableRef, $fileName, undef);
}
=pod

=head1 The procedure stores table to a file

=over 4

=item table2fileEx($tableRef, $fileName, $header)

The procedure stores table to a file.

Parameters:
    $tableRef     - ARRAYref of ARRAYS reference to table  
    $fileName     - file name where table should be stored
    $header       - header line(s) to be printed at the top of the file. 

Returns:
    table as text
    
=back

=cut

sub table2fileEx(\@;$;$) {
    croak "$0:ERROR: table2file\n" unless ( @_ == 3 );
    my ( $tableRef, $fileName, $header ) = @_;

    #print "table2file $fileName\n";
    my $rowId   = 0;
    my $content = '';
    foreach my $rowRef ( @{$tableRef} ) {
        my $cellId = 0;
        my $found  = 0;
        foreach my $cell ( @{$rowRef} ) {
            $found = 1;
            if ( $cellId == 0 ) {
                $content .= "$cell";
            }    # if
            else {
                $content .= "\t$cell";
            }    # else
            $cellId++;
        }    # foreach
        if ( $found == 1 ) {
            $content .= "\n";
            $rowId++;
        }
    }    # foreach
    my $file = IO::File->new( ">" . $fileName )
      || die
      "ERROR: table2file() Couldn't create/open file handle for $fileName $!\n";

    if (defined $header)
    {
        print $file $header;
    }
    print $file $content;
    close $file;

    return $content;
}

=pod

=head1 The procedure loads table from a file. It ignores file lines beginning with #

=over 4

=item file2table($tableRef, $fileName)

The procedure loads table from a file.

Parameters:
    $tableRef     - ARRAYref of ARRAYS reference to table  
    $fileName     - file name where table should be stored 

Returns:
    table as text
    
=back

=cut

sub file2table(\@;$) {
    croak "ERROR: file2table\n" unless ( @_ == 2 );
    my ( $tableRef, $fileName ) = @_;
    if ( -d $fileName ) {
        die
"ERROR: file2table() $fileName is directory. Path to a file expected.\n";
    }
    open( FILE, "<$fileName" )
      || die "ERROR: file2table() Couldn't open file handle for $fileName $!\n";
    while (<FILE>) {
        my $line = $_;
        next if $line =~ /^\#/;
        chomp $line;
        my @row = split "\t", $line;
        push @{$tableRef}, \@row;
    }    # while
    close FILE;
}

=pod

=head1 

Loads a table from a file.  It ignores file lines
beginning with "#" which occur in a continuous block at the start of
the file, except for lines beginning with "#$", where the remainder of
the line is treated as a key-value pair which is parsed and stored in
the metaDataRef hash.

=over 4

=item file2tableEx($tableRef, $metaDataRef, $fileName)

The procedure loads table from a file.

Parameters:
    $tableRef     - ARRAYref of ARRAYS reference to table
    $metaDataRef  - HASHref of header meta-data
    $fileName     - file name where table should be stored 

Returns:
    table as text
    
=back

=cut

sub file2tableEx(\@;\%;$) {
    croak "ERROR: file2tableEx\n" unless ( @_ == 3 );
    my ( $tableRef, $metaDataRef, $fileName ) = @_;
    if ( -d $fileName ) {
        die
"ERROR: file2tableEx() $fileName is directory. Path to a file expected.\n";
    }
    open( FILE, "<$fileName" )
      || die "ERROR: file2tableEx() Couldn't open file handle for $fileName $!\n";

    my $isHeader = 1;
    while (<FILE>) {
        my $line = $_;
        chomp $line;
        if( $isHeader ) {
            if($line =~ /^#/) {
                if($line =~ /^#\$/) {
                    $line =~ s/^#\$\s+//;
                    my ($k,$v) = split(/\s+/,$line,2);
                    $$metaDataRef{$k} = $v;
                }
                next;
            }
            $isHeader = 0;
        }
        my @row = split "\t", $line;
        push @{$tableRef}, \@row;
    }    # while
    close FILE;
}

=pod

=head1 The procedure reads HASH MAP XML config file.

=over 4

=item readXML($confRef, $fileName, $dieOnFailure=0)

The procedure reads HASH MAP XML config file..

Parameters:
    $confRef        - reference to hash map with all runs configuration 
    $fileName       - config file name
    $dieOnFailure   - if set, terminates the process in case file does not exist
                      default - 0
   
Returns:
    nothing
    
=back

=cut

sub readXML(\%;$;$) {
    croak "ERROR: readXML\n" unless ( @_ == 2 || @_ == 3 );
    my ( $confRef, $fileName, $dieOnFailure ) = @_;
    my %keyAttr = ();
    my $xs      = '';
    if ( -d $fileName ) {
        errorExit
"ERROR: readXML() $fileName is directory. Path to xml file expected.\n";
    }
    $xs = new XML::Simple(
        searchpath => ".",
        forcearray => 1,
    );

    if ( !-e $fileName ) {
        my $message = "readXML() file $fileName doesn't exist.";
        errorExit "ERROR:$message\n" unless ( defined $dieOnFailure and !$dieOnFailure );
        logWarning $message;
        return 0;
    }
    else{
        my $ref = $xs->XMLin($fileName);
        my ( $key, $value );
        foreach my $key ( keys %{$ref} ) {
            ${$confRef}{$key} = ${$ref}{$key};
        }
    }
    return 1;
}

=pod

=head1 The procedure writes HASH MAPto the xml run file.

=over 4

=item writeXML($confRef, $fileName)

The procedure writes HASH MAP with runs config to the run xml config file.

Parameters:
    $hamMapRef     - reference to hash map with all runs configuration 
    $fileName    - config file name
   
Returns:
    nothing
    
=back

=cut

sub writeXML(\%;$) {
    croak "ERROR: writeXML\n" unless ( @_ == 2 || @_ == 3 );
    my ( $hamMapRef, $fileName, $keyAttrRef ) = @_;
    my %keyAttr = ();
    my $xs      = '';
    if ( -d $fileName ) {
        errorExit
"ERROR: writeXML() $fileName is directory. Path to xml file expected.\n";
    }
    if ( defined $keyAttrRef ) {
        %keyAttr = %{$keyAttrRef};
        $xs      = new XML::Simple(
            searchpath => ".",
            KeyAttr    => %keyAttr,
            forcearray => 1,
        );
    }
    else {
        $xs = new XML::Simple(
            searchpath => ".",
            forcearray => 1,
        );
    }
    my $xml  = $xs->XMLout($hamMapRef);
    my $file = IO::File->new( ">" . $fileName )
      || errorExit
"ERROR: writeRunsConfig() Couldn't create/open file handle for $fileName $!\n";
    print $file $xml;
    close $file;
}

=pod

=head1 The procedure counts number of lines from command ouput

=over 4

=item wcl($cmd)

The procedure counts number of lines from command oupu

Parameters:
    $cmd - command to be executed

Returns:
    nothing
=back

=cut

sub wcl($) {
    my ($cmd)   = @_;
    my $command = "$cmd | wc -l";
    my $count   = 0;
    open( CMD, "$command |" )
      || errorExit "Could not $command $! \n";

    while (<CMD>) {
        chomp($_);
        $count = $_;
    }
    close(CMD);
    return $count;
}

=pod

=head1 Returns a string which is a concatenation of mean +/- stdev

=over 4

=item formatValueWithDeviation($xml, $element)

The procedure counts number of lines from command oupu

Parameters:
    $xml - hashref with mean and stdev keys

Returns:
    nothing
=back

=cut
sub formatValueWithDeviation($;) {
    croak "ERROR: formatValueWithDeviation " unless ( @_ == 1 );
    my ( $element ) = @_;
    
    my $ret = defined ($element->{mean}) ?  
                $element->{mean} . " +/- " . $element->{stdev} : 0;

    return $ret;    
}


1;    # says use was ok
__END__

