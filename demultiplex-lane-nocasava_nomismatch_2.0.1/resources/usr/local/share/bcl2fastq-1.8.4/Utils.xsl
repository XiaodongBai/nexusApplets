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
xmlns:math="http://exslt.org/math"
xmlns:ckbk="http://www.oreily.com/xsltckbk"
> 

<xsl:template name="generateSequence">
    <xsl:param name="var"/>
    <xsl:param name="stop"/>
    <xsl:param name="step"/>

    <!--xsl:message><xsl:value-of select="concat($var, '!=', $stop, ':', $step)"/></xsl:message-->
    
    <xsl:choose>
        <xsl:when test="$var != $stop">
            <xsl:value-of select="concat($var, ' ')"/>
            <xsl:call-template name="generateSequence">
                <xsl:with-param name="var" select="number($var) + number($step)"/>
                <xsl:with-param name="stop" select="$stop"/>
                <xsl:with-param name="step" select="$step"/>
            </xsl:call-template>
        </xsl:when>
        <xsl:otherwise></xsl:otherwise>
    </xsl:choose>
    
</xsl:template>

<xsl:template name="ckbk:variance">
  <xsl:param name="nodes" select="/.."/>
  <xsl:param name="sum" select="0"/>
  <xsl:param name="sum-sq" select="0"/>
  <xsl:param name="count" select="0"/>
  <xsl:choose>
    <xsl:when test="not($nodes)">
      <xsl:value-of select="($sum-sq - ($sum * $sum) div $count) div ($count - 1)"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:variable name="value" select="$nodes[1]"/>
      <xsl:call-template name="ckbk:variance">
        <xsl:with-param name="nodes" select="$nodes[position( ) != 1]"/>
        <xsl:with-param name="sum" select="$sum + $value"/>
        <xsl:with-param name="sum-sq" select="$sum-sq + ($value * $value)"/>
        <xsl:with-param name="count" select="$count + 1"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>


<xsl:template name="mean">
  <xsl:param name="nodes" select="/.."/>
  <xsl:param name="sum" select="0"/>
  <xsl:param name="count" select="0"/>
  <xsl:choose>
    <xsl:when test="not($nodes)">
      <xsl:value-of select="$sum div $count"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:variable name="value" select="$nodes[1]"/>
      <xsl:call-template name="mean">
        <xsl:with-param name="nodes" select="$nodes[position( ) != 1]"/>
        <xsl:with-param name="sum" select="$sum + $value"/>
        <xsl:with-param name="count" select="$count + 1"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>


<xsl:template name="stdev">
  <xsl:param name="nodes"/>
  <xsl:param name="mean" select="0"/>
  <xsl:param name="sum" select="0"/>
  <xsl:param name="sum-sq" select="0"/>
  <xsl:param name="count" select="0"/>
  <!--xsl:message>tada <xsl:for-each select="$nodes"><xsl:value-of select="."/></xsl:for-each></xsl:message-->
  <xsl:choose>
    <xsl:when test="not($nodes)">
      <xsl:choose>
        <xsl:when test="3 &gt; $count">
            <xsl:value-of select="''"/> <!-- don't report standard deviation for sets where it makes no sense-->
        </xsl:when>
        <xsl:otherwise>
            <xsl:value-of select="math:sqrt($sum-sq div $count)"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:when>
    <xsl:otherwise>
      <xsl:variable name="value" select="$nodes[1]"/>
      <xsl:call-template name="stdev">
        <xsl:with-param name="nodes" select="$nodes[position( ) != 1]"/>
        <xsl:with-param name="mean" select="$mean"/>
        <xsl:with-param name="sum" select="$sum + $value"/>
        <xsl:with-param name="sum-sq" select="$sum-sq + (($mean - $value) * ($mean - $value))"/>
        <xsl:with-param name="count" select="$count + 1"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>


<xsl:template name="meanAndStdev">
  <xsl:param name="nodes" select="/.."/>
  <xsl:param name="meanFmt" select="'0.00'"/>
  <xsl:param name="devFmt" select="'0.00'"/>
  <xsl:choose>
    <xsl:when test="not($nodes)">
      <xsl:value-of select="0"/>
    </xsl:when>
    <xsl:otherwise>

      <xsl:variable name="m">
          <xsl:call-template name="mean">
            <xsl:with-param name="nodes" select="$nodes"/>
          </xsl:call-template>
      </xsl:variable>

      <xsl:variable name="std">
          <xsl:call-template name="stdev">
            <xsl:with-param name="nodes" select="$nodes"/>
            <xsl:with-param name="mean" select="$m"/>
          </xsl:call-template>
      </xsl:variable>

      <xsl:variable name="roundedMean" select="format-number($m, $meanFmt)"/>
      <xsl:variable name="roundedStdev" select="format-number($std, $devFmt)"/>
      <xsl:choose>
          <xsl:when test="string($std) = ''">
            <xsl:value-of select="$roundedMean"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="concat($roundedMean, ' +/- ', $roundedStdev)"/>
          </xsl:otherwise>
      </xsl:choose>

    </xsl:otherwise>
  </xsl:choose>
</xsl:template>


</xsl:stylesheet>
