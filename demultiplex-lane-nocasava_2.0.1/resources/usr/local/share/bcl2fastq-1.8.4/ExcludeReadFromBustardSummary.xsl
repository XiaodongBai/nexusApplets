<?xml version="1.0"?> 
<!--
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
-->
<xsl:stylesheet version="1.0" 
xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
xmlns:str="http://exslt.org/strings"
> 

<xsl:output method="xml" version="1.0" indent="yes"/>

<!--
    The variables below are used as strings containing single character per number.
    Make sure they don't contain any garbage such as tabs, spaces, etc
-->
<xsl:variable name="originalReadNumbers" select="translate(normalize-space($INCLUDED_ORIGINAL_READS_PARAM), ' ', '')"/>
<xsl:variable name="demuxReadNumbers" select="translate(normalize-space($DEMUX_READS_PARAM), ' ', '')"/>
<!--xsl:variable name="excludedReadNumbers" select="translate(normalize-space($EXCLUDED_READS_PARAM), ' ', '')"/-->

<xsl:template match="*"> 
    <xsl:copy>
        <xsl:apply-templates/>
    </xsl:copy>
</xsl:template>

<xsl:template match="Read">
    <xsl:choose>
        <xsl:when test="not(contains($originalReadNumbers, readNumber))">
        </xsl:when>
        <xsl:otherwise>
            <xsl:copy>
                <xsl:apply-templates/>
            </xsl:copy>
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<xsl:template match="readNumber">
   <xsl:element name="readNumber">
       <xsl:value-of select="translate(current(), $originalReadNumbers, $demuxReadNumbers)"/>
   </xsl:element>
</xsl:template>

</xsl:stylesheet>
