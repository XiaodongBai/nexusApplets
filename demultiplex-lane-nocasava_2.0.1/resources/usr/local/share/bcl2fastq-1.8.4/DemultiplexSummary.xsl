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
xmlns:casava="http://www.illumina.com/casava/alignment" 
xmlns:str="http://exslt.org/strings"
> 

<xsl:output method="html" version="4.0" indent="yes"/>

<xsl:include href="SampleSheet.xsl"/>
<!--xsl:include href="SampleSheet.xsl"/-->

<xsl:variable name="sampleSheetXml" select="document($DEMUX_CONFIG_XML_PATH_PARAM)/DemultiplexConfig"/>
<xsl:variable name="demuxSummaryXml" select="/"/>

<xsl:template match="/"> 

<html>
<link rel="stylesheet" href="css/Reports.css" type="text/css"/>
<body>

    <h1>Flowcell: <xsl:value-of select="$sampleSheetXml/FlowcellInfo/@ID"/></h1>

    <h2>Barcode lane statistics</h2>

    <div ID="ScrollableTableHeaderDiv">
    <table width="100%">
    <col width="4%"/><col width="5%"/><col width="19%"/><col width="8%"/><col width="7%"/><col width="5%"/><col width="12%"/><col width="7%"/><col width="4%"/><col width="5%"/><col width="4%"/><col width="5%"/><col width="6%"/><col width="5%"/><col/>
    <tr>
        <th>Lane</th>
        <th>Sample ID</th>
        <th>Sample Ref</th>
        <th>Index</th>
        <th>Description</th>
        <th>Control</th>
        <th>Project</th>
        <th>Yield (Mbases)</th>
        <th>% PF</th>
        <th># Reads</th>
        <th>% of raw clusters per lane</th>
        <th>% Perfect Index Reads</th>
        <th>% One Mismatch Reads (Index)</th>
        <th>% of >= Q30 Bases (PF)</th>
        <th>Mean Quality Score (PF)</th>
    </tr>
    </table>
    </div>

    <div ID="ScrollableTableBodyDiv">
    <table width="100%">
    <col width="4%"/><col width="5%"/><col width="19%"/><col width="8%"/><col width="7%"/><col width="5%"/><col width="12%"/><col width="7%"/><col width="4%"/><col width="5%"/><col width="4%"/><col width="5%"/><col width="6%"/><col width="5%"/><col/>
    <xsl:for-each select="$sampleSheetXml/FlowcellInfo/Lane/Sample">
    <xsl:sort select="concat(../@Number, @SampleId, '@Ref', @Index)"/>
        <xsl:variable name="lane" select="../@Number"/>
        <xsl:variable name="project" select="@ProjectId"/>
        <xsl:variable name="projectFolderName"><xsl:choose>
            <xsl:when test="'Undetermined_indices'=$project">Undetermined_indices</xsl:when>
            <xsl:otherwise><xsl:value-of select="concat('Project_', $project)"/></xsl:otherwise>
        </xsl:choose></xsl:variable>
        <xsl:variable name="reference" select="@Ref"/>
        <xsl:variable name="sample" select="@SampleId"/>
        <xsl:variable name="barcode" select="@Index"/>
        <xsl:variable name="description" select="@Desc"/>
        <xsl:variable name="control" select="@Control"/>

        <xsl:variable name="laneSample" select="$demuxSummaryXml/Summary/Lane[@index=$lane]/Sample"/>
        <xsl:variable name="laneBarcode" select="$laneSample/Barcode[@index=$barcode]/Tile/Read"/>
        <xsl:variable name="laneBarcodes" select="$laneSample/Barcode/Tile/Read"/>

        <xsl:variable name="clustersRaw" select="sum($laneBarcode/Raw/ClusterCount)"/>

        <xsl:variable name="clustersRawLane" select="sum($laneBarcodes/Raw/ClusterCount)"/>

        <xsl:variable name="yieldPf" select="sum($laneBarcode/Pf/Yield)"/>
        <xsl:variable name="clustersPf" select="sum($laneBarcode/Pf/ClusterCount)"/>
        <xsl:variable name="clusters0MismatchBarcodePf" select="sum($laneBarcode/Pf/ClusterCount0MismatchBarcode)"/>
        <xsl:variable name="clusters1MismatchBarcodePf" select="sum($laneBarcode/Pf/ClusterCount1MismatchBarcode)"/>
        <xsl:variable name="yieldQ30Pf" select="sum($laneBarcode/Pf/YieldQ30)"/>
        <xsl:variable name="qscoresSumPf" select="sum($laneBarcode/Pf/QualityScoreSum)"/>

        <tr>
            <td><xsl:value-of select="$lane"/></td>
            <td><xsl:value-of select="$sample"/></td>
            <td><xsl:value-of select="$reference"/></td>
            <td><xsl:value-of select="$barcode"/></td>
            <td><xsl:value-of select="$description"/></td>
            <td><xsl:value-of select="$control"/></td>
            <td><xsl:value-of select="$project"/></td>
            <td><xsl:value-of select="format-number(round($yieldPf div 1000000), '###,###,###,###,###')"/></td>
            <td><xsl:if test="0 != $clustersRaw"><xsl:value-of select="format-number($clustersPf div $clustersRaw * 100, '0.00')"/></xsl:if></td>
            <td><xsl:value-of select="format-number($clustersRaw, '###,###,###,###,###')"/></td>
            <td><xsl:if test="0 != $clustersRawLane"><xsl:value-of select="format-number($clustersRaw div $clustersRawLane * 100, '0.00')"/></xsl:if></td>
            <td><xsl:if test="0 != $clustersPf"><xsl:value-of select="format-number($clusters0MismatchBarcodePf div $clustersPf * 100, '0.00')"/></xsl:if></td>
            <td><xsl:if test="0 != $clustersPf"><xsl:value-of select="format-number($clusters1MismatchBarcodePf div $clustersPf * 100, '0.00')"/></xsl:if></td>
            <td><xsl:if test="0 != $yieldPf"><xsl:value-of select="format-number($yieldQ30Pf div $yieldPf * 100, '0.00')"/></xsl:if></td>
            <td><xsl:if test="0 != $yieldPf"><xsl:value-of select="format-number($qscoresSumPf div $yieldPf, '0.00')"/></xsl:if></td>
        </tr>
    </xsl:for-each>
    </table>
    </div>

    <p/>
    <h2>Sample information</h2>

    <div ID="ScrollableTableHeaderDiv">
    <table width="100%">
    <col width="10%"/><col width="10%"/><col width="7%"/><col />
    <tr>
        <th>Sample<p/>ID</th>
        <th>Recipe</th>
        <th>Operator</th>
        <th>Directory</th>
    </tr>
    </table>
    </div>

    <div ID="ScrollableTableBodyDiv">
    <table width="100%">
    <col width="10%"/><col width="10%"/><col width="7%"/><col />
    <xsl:for-each select="$sampleSheetXml/FlowcellInfo/Lane/Sample[generate-id()=generate-id(key('samplesById', @SampleId)[1])]">
    <xsl:sort select="concat(../@Number, @SampleId, '@Ref', @Index)"/>
        <xsl:variable name="project" select="@ProjectId"/>
        <xsl:variable name="projectFolderName"><xsl:choose>
            <xsl:when test="'Undetermined_indices'=$project">Undetermined_indices</xsl:when>
            <xsl:otherwise><xsl:value-of select="concat('Project_', $project)"/></xsl:otherwise>
        </xsl:choose></xsl:variable>
        <xsl:variable name="sample" select="@SampleId"/>

        <tr>
            <td><xsl:value-of select="$sample"/></td>
            <td><xsl:value-of select="$sampleSheetXml/FlowcellInfo/@Recipe"/></td>
            <td><xsl:value-of select="$sampleSheetXml/FlowcellInfo/@Operator"/></td>
            <td><xsl:value-of select="concat($PROJECTS_ROOT_PATH_PARAM, '/', $projectFolderName, '/', 'Sample_', $sample)"/></td>
        </tr>
    </xsl:for-each>
    </table>
    </div>


<p>bcl2fastq-1.8.4</p>
</body>
</html>

</xsl:template>


</xsl:stylesheet>
