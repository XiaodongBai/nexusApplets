# PROJECT: CASAVA
# MODULE:  $RCSfile: Template.pm,v $
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#

=pod

=head1 NAME

Casava::Common::Template.pm - Perl utility library for generating content 
    (html) from templates; 

=head1 SYNOPSIS

# include what functions you need... 
use Casava::Common::Templates qw();  

=head1 DESCRIPTION

Exports: 

    
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

package Casava::Common::Template;

#
# Place functions/variables you want to *export*, ie be visible
# from the caller package into @EXPORT_OK
#
BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT = qw($VERSION &configure &fromTemplate &fromHtmlTemplate
      &fromHtmlTemplate2HtmlFile &htmlTable &htmlTableWithTotals);
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);
use Carp;
use IO::File;

use Casava::Common::Log;
sub configure($);
sub fromTemplate(\%;$);
sub fromHtmlTemplate(\%;$);
sub fromHtmlTemplate2HtmlFile(\%;$;$);
sub htmlTable(\@;$;);
sub htmlTableWithTotals(\@;\@;$);

# Global variable
our %templConfig   = ();
our %defaultValues = ();
$defaultValues{userName} = 'CASAVA';

=pod

=head1 The procedure configures template library.

=over 4

=item configure($templPath)

The procedure configures template library. e.g.: sets the path to template files.

Parameters:
    $templPath      - path to templates files 

Returns:
    nothing
    
=back

=cut

sub configure($) {
    croak "$0:ERROR: configure\n" unless ( @_ == 1 );
    my ($templPath) = @_;
    $templConfig{templPath} = $templPath;
}

=pod

=head1 The procedure configures template library.

=over 4

=item configure($templPath)

The procedure configures template library. e.g.: sets the path to template files.

Parameters:
    $valuesRef    - HASH MAP ref with key value pairs to replace keys in the template 
    $templateName - template name

Returns:
    some content
    
=back

=cut

sub fromTemplate(\%;$) {
    croak "$0:ERROR: fromTemplate\n" unless ( @_ == 2 );
    my ( $valuesRef, $templateName ) = @_;
    my $templFile =
      File::Spec->catfile( $templConfig{templPath}, $templateName );
    unless ( defined $templConfig{templPath} ) {
        errorExit "ERROR: fromTemplate() - Template library not configured.\n";
    }
    my $templContent = '';
    open( FILE, "<$templFile" )
      || errorExit "ERROR: $0::Couldn't open file handle for $templFile $!\n";
    while (<FILE>) {
        my $line = $_;
        $templContent .= $line;
    }    # while
    close FILE;
    foreach my $templKey ( keys %{$valuesRef} ) {
        $templContent =~ s/\%$templKey/$valuesRef->{$templKey}/g;
    }    #foreach
    return $templContent;
}

=pod

=head1 The procedure generates the html page from templates.

=over 4

=item fromHtmlTemplate($valuesRef, $templateName)

The procedure generates the html page from templates.

Parameters:
    $valuesRef    - HASH MAP ref with key value pairs to replace keys in the template 
    $templateName - template name (all templates have a default ext. html) 

Returns:
    html content
    
=back

=cut

sub fromHtmlTemplate(\%;$) {
    croak "$0:ERROR: fromHtmlTemplate\n" unless ( @_ == 2 );
    my ( $valuesRef, $templateName ) = @_;
    my $html = '';
    foreach my $key ( keys %defaultValues ) {
        unless ( defined $valuesRef->{$key} ) {
            $valuesRef->{$key} = $defaultValues{$key};
        }    #unless
    }    # foreach
    unless ( defined $valuesRef->{date} ) {
        $valuesRef->{date} = strftime '%d-%m-%y', localtime;
    }    #unless
    unless ( defined $valuesRef->{time} ) {
        $valuesRef->{time} = strftime '%H:%M:%S %d-%m-%y', localtime;
    }    #unless
    $html .= fromTemplate( %{$valuesRef}, 'header.html' );
    $html .= fromTemplate( %{$valuesRef}, $templateName . '.html' );
    $html .= fromTemplate( %{$valuesRef}, 'footer.html' );
    return $html;
}

=pod

=head1 The procedure generates the html page from templates to file.

=over 4

=item fromHtmlTemplate2HtmlFile($valuesRef, $templateName, $fileName)

The procedure generates the html page from templates.

Parameters:
    $valuesRef    - HASH MAP ref with key value pairs to replace keys in the template 
    $templateName - template name (all templates have a default ext. html) 
    $fileName     - file name were content of html will be stored 

Returns:
    html content
    
=back

=cut

sub fromHtmlTemplate2HtmlFile(\%;$;$) {
    croak "ERROR: fromHtmlTemplate\n" unless ( @_ == 3 );
    my ( $valuesRef, $templateName, $fileName ) = @_;
    my $html = fromHtmlTemplate( %{$valuesRef}, $templateName );
    my $file = IO::File->new( ">" . $fileName )
      || errorExit "ERROR: $0::Couldn't create/open file handle for $fileName $!\n";
    print $file $html;
    close $file;
    return $html;
}

=pod

=head1 The procedure generates the html table from structure.

=over 4

=item htmlTable($tableRef, $styleName)

The procedure generates the html table from structure.

Parameters:
    $tableRef     - ARRAYref of ARRAYS reference to table  
    $styleName    - style to be used (e.g.: $styleName(Th) $styleName(Td) $styleName(Table)) 

Returns:
    html content
    
=back

=cut

sub htmlTable(\@;$;) {
    croak "$0:ERROR: htmlTable\n" unless ( @_ == 2 );
    my ( $tableRef, $styleName ) = @_;

    my @totalsRow = ();
    my $html  = htmlTableWithTotals(@$tableRef, @totalsRow, $styleName);
    
    return $html;
}

=pod

=head1 The procedure generates the html table from structure.

=over 4

=item htmlTableWithTotals($tableRef, $totalsRowRef, $styleName)

The procedure generates the html table from structure.

Parameters:
    $tableRef     - ARRAYref of ARRAYS reference to table  
    $totalsRowRef - ARRAYref of values which are added to the bottom with header style
                    if referenced array is empty, totals are not added.
    $styleName    - style to be used (e.g.: $styleName(Th) $styleName(Td) $styleName(Table)) 

Returns:
    html content
    
=back

=cut

sub htmlTableWithTotals(\@;\@;$) {
    croak "$0:ERROR: htmlTableWithTotals\n" unless ( @_ == 3 );
    my ( $tableRef, $totalsRowRef, $styleName ) = @_;
    my $html  = "<TABLE class='" . $styleName . "Table'>\n";
    my $rowId = 0;

    foreach my $rowRef ( @{$tableRef} ) {
        $html .= " <TR>";
        if ( $rowId == 0 ) {
            foreach my $cell ( @{$rowRef} ) {
                $html .= "<TH class='" . $styleName . "Th'>$cell</TH>";
            }# foreach
        } else {
            foreach my $cell ( @{$rowRef} ) {
                $html .= "<Td class='" . $styleName . "Td'>" . (defined $cell ? $cell : 'undef') . "</Td>\n";
            }# foreach
        }
        $html .= "</TR>\n";
        $rowId++;
    } # foreach

    if (scalar @$totalsRowRef)
    {
        $html .= " <TR>";
        foreach my $cell ( @{$totalsRowRef} ) {
            $html .= "<TH class='" . $styleName . "Th'>$cell</TH>\n";
        }
        $html .= "</TR>\n";
    }

    $html .= "</TABLE>\n";
    #print $html;     
    return $html;
}

1;                                 # says use was ok
__END__

