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

<xsl:include href="Utils.xsl"/>

<xsl:output method="xml" version="1.0" indent="yes"/>
<xsl:strip-space elements="*"/>

<!--
    The variables below are used as strings containing single character per number.
    Make sure they don't contain any garbage such as tabs, spaces, etc
-->
<xsl:variable name="demuxReadNumbers" select="translate(normalize-space($DEMUX_READS_PARAM), ' ', '')"/>
<xsl:variable name="demuxFirstCycles" select="translate(normalize-space($READ_DEMUX_FIRST_CYCLES_PARAM), ' ', '')"/>
<xsl:variable name="demuxLastCycles" select="translate(normalize-space($READ_DEMUX_LAST_CYCLES_PARAM), ' ', '')"/>

<xsl:template match="*"> 
    <xsl:copy>
        <xsl:copy-of select="@*"/>
        <xsl:apply-templates/>
    </xsl:copy>
</xsl:template>

<xsl:template match="RunParameters">
    <xsl:copy>
        <xsl:copy-of select="*[local-name() != 'Reads']"/>
        <xsl:variable name="demuxReadNumbersNodes" select="str:tokenize($DEMUX_READS_PARAM)"/>
        <xsl:variable name="demuxFirstCyclesNodes" select="str:tokenize($READ_DEMUX_FIRST_CYCLES_PARAM)"/>
        <xsl:variable name="demuxLastCyclesNodes" select="str:tokenize($READ_DEMUX_LAST_CYCLES_PARAM)"/>

        <xsl:variable name="counter">
            <xsl:call-template name="generateSequence">
                <xsl:with-param name="var" select="1"/>
                <xsl:with-param name="stop" select="count(str:tokenize($DEMUX_READS_PARAM)) + 1"/>
                <xsl:with-param name="step" select="1"/>
            </xsl:call-template>
        </xsl:variable>

        <xsl:for-each select="str:tokenize($counter)">
            <xsl:variable name="demuxReadIndex" select="number(.)"/>
            <!--xsl:message><xsl:value-of select ="$demuxReadIndex"/></xsl:message-->
            <xsl:element name="Reads">
                <xsl:attribute name="Index">
                    <xsl:value-of select="$demuxReadNumbersNodes[$demuxReadIndex]"/>
                </xsl:attribute>
                <xsl:element name="FirstCycle">
                    <xsl:value-of select="$demuxFirstCyclesNodes[$demuxReadIndex]"/>
                </xsl:element>
                <xsl:element name="LastCycle">
                    <xsl:value-of select="$demuxLastCyclesNodes[$demuxReadIndex]"/>
                </xsl:element>
            </xsl:element>
        </xsl:for-each>

    </xsl:copy>
</xsl:template>

</xsl:stylesheet>
