# PROJECT: CASAVA
# MODULE:  $RCSfile: Log.pm,v $
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008, 2009 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#

=pod

=head1 NAME

Casava::Common::Log.pm - Perl utility library for reporting messages and errors.

=head1 SYNOPSIS

# include what functions you need... 
use Casava::Common::Log qw();  

=head1 DESCRIPTION

=head1 AUTHORSHIP

Copyright (c) 2008, 2009 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

Created by Lukasz Szajkowski <lszajkowski@illumina.com>

=cut

package Casava::Common::Log;

#
# Place functions/variables you want to *export*, ie be visible
# from the caller package into @EXPORT_OK
#
BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT = qw($VERSION
      &initLog &printLog logDebug
      &logInfo &logWarning &getWarningsCount &errorExit);
    @EXPORT_OK = qw();
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);
use XML::Simple;
use Carp;
use File::Spec;
use Term::ANSIColor;
use IO::File;

sub initLog ($$$;$;$);
sub printLog ($$;$);
sub logWarning ($);
sub logInfo ($;$);
sub logDebug ($;$);
sub errorExit ($);

# Global variable
my $defaultVerbose    = 0;
my $defaultLogVerbose = 0;
my $defaultLogPath    = undef;
my $warningsCount = 0;
my $defaultColorMode = 0;
my $defaultPrintTimeStamps = 1;
if (defined $ENV{CASAVA_LOG_LEVEL} && $ENV{CASAVA_LOG_LEVEL} =~ /^\d+$/)
{
    $defaultPrintTimeStamps = 0;
    $defaultVerbose = $ENV{CASAVA_LOG_LEVEL};
}

=pod

=head1 Initialize log file to be writable for all processes

=over 4

=item initLog($message, $verbose)

B<Parameters:>
    $logPath        - path to a log file
    $verbose        - default verbose level
    $logVerbose     - default verbose level for a log file

B<Returns:>
    nothing
=back

=cut

sub initLog ($$$;$;$) {
    my ( $logPath, $verbose, $logVerbose, $colorMode, $printTimeStamps) = @_;

    # CASAVA_LOG_LEVEL is an environment variable that is set by workflow makefiles.
    # disregard everything in favour of CASAVA_LOG_LEVEL if it is set!
    unless (defined $ENV{CASAVA_LOG_LEVEL} && $ENV{CASAVA_LOG_LEVEL} =~ /^\d+$/)
    {
        $defaultVerbose    = $verbose;
        $defaultLogVerbose = $logVerbose;
        $defaultLogPath    = $logPath;
        $defaultColorMode = $colorMode if defined $colorMode;
        $defaultPrintTimeStamps = $printTimeStamps if defined $printTimeStamps;
    }
}


=pod

=item sub getWarningsCount()

Returns number of warnings logged by L<logWarnings>.

B<Parameters:>

    none
    
B<Returns:> 
    number of warnings
    
=cut

sub getWarningsCount () {
    return $warningsCount;
}


=pod

=head1 The prints message to a stdout or a file.

=over 4

=item printLog($message, $verbose)

The procedure prints the message to standart output or to a file. 
Global variable $defaultVerbose controls printing.

B<Parameters:>
    $message - message to be printed
    verbose - verbose level (optional)

B<Returns:> 
    nothing
=back

=cut

my %AnsiColorFromMessageType = (
    info => 'dark green',
    warning => 'dark red',
    error => 'bold red',
    message => 'blue',
);

sub printLog ($$;$) {

    my ( $message, $verbose, $color ) = @_;
    my $timeStr     = '';
    $timeStr = ("[" . ( strftime "%Y-%m-%d %H:%M:%S", localtime ) . "]\t") if $defaultPrintTimeStamps;
    my $scriptStr     = "[".(File::Spec->splitpath ($0))[2]."]";
    my $timeScriptStr ="$timeStr$scriptStr\t"; 
    my $fullMessage = "$timeScriptStr$message";

    # print "printLog: [$defaultLogPath] $fullMessage\n";
    if ( !defined $verbose || $defaultVerbose >= $verbose ) {
        if ($defaultColorMode && defined $color)
        {
            print STDERR $timeScriptStr;
            
            use Term::ANSIColor qw (:constants);
            print STDERR color $AnsiColorFromMessageType{$color};
            $message =~ s/\n*$//;
            print STDERR $message;
            print STDERR RESET;
            print STDERR "\n";
        }
        else
        {
            print STDERR "$fullMessage";
        }
    }

    if ( !defined $verbose || $defaultLogVerbose >= $verbose ) {
        if ( $defaultLogPath ) {
            my $file = IO::File->new( ">>" . $defaultLogPath )
              || die
"ERROR: printLog() Couldn't create/open file handle for $defaultLogPath $!\n";
            print $file "$fullMessage";
            close $file;
        }

    }
}

=pod

=item logWarning($msg)

The procedure logs a warning message. 
Global variable $defaultVerbose controls printing.

B<Parameters:>

    $warning - warning to be printed
    
B<Returns:> 
    nothing
    
=cut

sub logWarning ($) {
    croak "ERROR: logWarning " . join( ',', @_ ) unless ( 1 == @_ );
    my ( $warning ) = @_;

    printLog( "WARNING: ".$warning."\n", 0, 'warning');
    ++$warningsCount;
}

=pod

=item logDebug

The procedure prints the debug log message. 
Global variable $defaultVerbose controls printing.

B<Parameters:>

    $warning - warning to be printed
    verbose  - verbose level (default is 2)
    
B<Returns:> 
    nothing
    
=cut

sub logDebug ($;$) {
    my ( $message, $verbose ) = @_;
    $verbose = 2 if !defined $verbose;
    printLog( $message."\n", $verbose, 'info');
}

=pod

=item logInfo

The procedure prints the information log message. 
Global variable $defaultVerbose controls printing.

B<Parameters:>

    $warning - warning to be printed
    verbose  - verbose level (default is 1)
    
B<Returns:> 
    nothing
    
=cut

sub logInfo ($;$) {
    my ( $message, $verbose ) = @_;
    $verbose = 1 if !defined $verbose;
    printLog( "INFO: ".$message."\n", $verbose, 'info');
}

=pod

=item errorExit($msg)

The procedure prints the error to standart error output or to a file. 
Next the procedure exists
Global variable $defaultVerbose controls printing.

B<Parameters:>

    $errorMessage - message to be printed
    
B<Returns:> 
    nothing
    
=cut

sub errorExit ($) {
    croak "ERROR: errorExit " . join( ',', @_ ) unless ( @_ == 1 );
    my ($errorMessage) = @_;
#    my $timeStr     = "[" . ( strftime "%d/%m/%y %H:%M:%S", localtime ) . "] ";
#    my $scriptStr     = "[".(File::Spec->splitpath ($0))[2]."]";
#    my $fullMessage = "$timeStr\t$scriptStr\t$errorMessage\n";

#    if ( $defaultLogPath ) {
#        my $file = IO::File->new( ">>" . $defaultLogPath )
#          || die "FATAL ERROR: printLog() Couldn't create/open "
#          . "file handle for $defaultLogPath $!\n";
##        my $backtraceMessage .= "$timeStr "
##          . Carp::longmess("Stack:")
##          . $fullMessage;
#        print $file $backtraceMessage;
#        close $file;
#    }
    printLog "$errorMessage\n", 0, 'error';
    printLog "BACKTRACE: " . Carp::longmess(), 0;
    die;
}
1;    # says use was ok
__END__
