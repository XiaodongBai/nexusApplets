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
xmlns:math="http://exslt.org/math"
xmlns:casava="http://www.illumina.com/casava/alignment"
> 

<xsl:key name="samplesById" match="Sample" use="@SampleId"/>
<xsl:key name="projectsById" match="@ProjectId" use="."/>

<xsl:template name="extractSampleBustardSummaryLaneAvg">
    <xsl:param name="bustardSummaryNode"/>
    <xsl:param name="project"/>
    <xsl:param name="sampleSheetXml"/>
    <xsl:param name="read"/>
    <xsl:param name="sample"/>
    <xsl:param name="lane" select="'*'"/> <!-- '*' - all lanes, otherwise filter by lane-->
    <xsl:param name="variable"/>

    <xsl:variable name="resultsByLane">
        <xsl:for-each select="$sampleSheetXml/FlowcellInfo/Lane['*'=$lane or @Number=$lane]/Sample[@ProjectId=$project and @SampleId=$sample]">
            <xsl:variable name="laneNumber" select="../@Number"/>

            <xsl:value-of select="$bustardSummaryNode/ExpandedLaneSummary/Read[readNumber=$read]/Lane[laneNumber=$laneNumber]/*[local-name()=$variable]"/>
            <xsl:value-of select="' '"/>
        </xsl:for-each>
    </xsl:variable>

    <xsl:call-template name="meanAndStdev">
        <xsl:with-param name="nodes" select="str:tokenize($resultsByLane)"/>
    </xsl:call-template>
</xsl:template>

<xsl:template name="extractSampleBustardSummaryExpLaneAvgFromAvg">
    <xsl:param name="bustardSummaryNode"/>
    <xsl:param name="project"/>
    <xsl:param name="sampleSheetXml"/>
    <xsl:param name="read"/>
    <xsl:param name="sample"/>
    <xsl:param name="lane" select="'*'"/> <!-- '*' - all lanes, otherwise filter by lane-->
    <xsl:param name="variable"/>

    <xsl:variable name="resultsByLane">
        <xsl:for-each select="$sampleSheetXml/FlowcellInfo/Lane['*'=$lane or @Number=$lane]/Sample[@ProjectId=$project and @SampleId=$sample]">
            <xsl:variable name="laneNumber" select="../@Number"/>

            <xsl:value-of select="$bustardSummaryNode/ExpandedLaneSummary/Read[readNumber=$read]/Lane[laneNumber=$laneNumber]/*[local-name()=$variable]/mean"/>
            <xsl:value-of select="' '"/>
        </xsl:for-each>
    </xsl:variable>

    <xsl:call-template name="meanAndStdev">
        <xsl:with-param name="nodes" select="str:tokenize($resultsByLane)"/>
    </xsl:call-template>
</xsl:template>


<xsl:template name="extractSampleBustardSummaryExpLaneOriginalMeanAndStdev">
    <xsl:param name="bustardSummaryNode"/>
    <xsl:param name="read"/>
    <xsl:param name="lane"/>
    <xsl:param name="variable"/>

    <xsl:value-of select="$bustardSummaryNode/ExpandedLaneSummary/Read[readNumber=$read]/Lane[laneNumber=$lane]/*[local-name()=$variable]/mean"/>
    <xsl:value-of select="' +/-'"/>
    <xsl:value-of select="$bustardSummaryNode/ExpandedLaneSummary/Read[readNumber=$read]/Lane[laneNumber=$lane]/*[local-name()=$variable]/stdev"/>

</xsl:template>


<xsl:template name="extractSampleBustardSummaryLaneOriginalMeanAndStdev">
    <xsl:param name="bustardSummaryNode"/>
    <xsl:param name="read"/>
    <xsl:param name="lane"/>
    <xsl:param name="variable"/>

    <xsl:value-of select="$bustardSummaryNode/LaneResultsSummary/Read[readNumber=$read]/Lane[laneNumber=$lane]/*[local-name()=$variable]/mean"/>
    <xsl:value-of select="' +/-'"/>
    <xsl:value-of select="$bustardSummaryNode/LaneResultsSummary/Read[readNumber=$read]/Lane[laneNumber=$lane]/*[local-name()=$variable]/stdev"/>

</xsl:template>

<xsl:template name="extractSampleBustardSummaryLaneAvgFromAvg">
    <xsl:param name="bustardSummaryNode"/>
    <xsl:param name="project"/>
    <xsl:param name="sampleSheetXml"/>
    <xsl:param name="read"/>
    <xsl:param name="sample"/>
    <xsl:param name="lane" select="'*'"/> <!-- '*' - all lanes, otherwise filter by lane-->
    <xsl:param name="variable"/>
    <xsl:param name="meanFmt" select="'0.00'"/>
    <xsl:param name="devFmt" select="'0.00'"/>

    <xsl:variable name="resultsByLane">
        <xsl:for-each select="$sampleSheetXml/FlowcellInfo/Lane['*'=$lane or @Number=$lane]/Sample[@ProjectId=$project and @SampleId=$sample]">
            <xsl:variable name="laneNumber" select="../@Number"/>

            <xsl:value-of select="$bustardSummaryNode/LaneResultsSummary/Read[readNumber=$read]/Lane[laneNumber=$laneNumber]/*[local-name()=$variable]/mean"/>
            <xsl:value-of select="' '"/>
        </xsl:for-each>
    </xsl:variable>

    <xsl:call-template name="meanAndStdev">
        <xsl:with-param name="nodes" select="str:tokenize($resultsByLane)"/>
        <xsl:with-param name="meanFmt" select="$meanFmt"/>
        <xsl:with-param name="devFmt" select="$devFmt"/>
    </xsl:call-template>
</xsl:template>

</xsl:stylesheet>
