<?xml version='1.0'  encoding="UTF-8" ?>
<xsl:stylesheet
	xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
	xmlns:x="http://exslt.org/dates-and-times"
	xmlns:foo="http://foo.com" exclude-result-prefixes="foo"
	version='1.0'
	>
<!--
############################
# Create a parameter file
#   from Illumina DemultiplexConfig (generated by CASAVA/configureBclToFastq.pl)
#   for NGS Workflow makefile
# Release: 0.9
# Date: 20/02/2014
# Author: Antony Le Béchec
############################

Options:
	count (int) : if set, only run on a subset of 'count' reads.
	outdir  (directory) : output directory
	indir (directory): input directory default: '.'
	reference (fasta): path to genome reference
Usage:
   xsltproc THIS.xsl DemultiplexConfig.xml

-->
<xsl:output method="text" encoding="UTF-8"/>

<xsl:key name="samples" match="Sample" use="@SampleId"/>

<!-- /(root) ==================================================================== -->
<xsl:template match="/">
<xsl:text>##################################
# PARAMETERS   
# Generated automatically
#   from a DemultiplexConfig.xml
##################################
</xsl:text>
	<xsl:apply-templates select="DemultiplexConfig"/>
</xsl:template>

<!-- DemultiplexConfig ==================================================================== -->
<xsl:template match="DemultiplexConfig">


<xsl:apply-templates select="Software"/>
<xsl:apply-templates select="FlowcellInfo"/>
</xsl:template>

<!-- Software ==================================================================== -->
<xsl:template match="Software">

	<xsl:variable name="Version" select="@Version"/>

	<xsl:text># Software Versions:
#    </xsl:text><xsl:value-of select="$Version"/><xsl:text>
</xsl:text>
	<xsl:for-each select="Software">
		<xsl:variable name="Name1" select="@Name"/>
		<xsl:variable name="Version1" select="@Version"/>
		<xsl:text>#    </xsl:text><xsl:value-of select="concat($Name1,'-',$Version1)"/><xsl:text> </xsl:text><xsl:text>
</xsl:text>
	</xsl:for-each>
	<xsl:for-each select="Software/Software">
		<xsl:variable name="Name2" select="@Name"/>
		<xsl:variable name="Version2" select="@Version"/>
		<xsl:text>#    </xsl:text><xsl:value-of select="concat($Name2,'-',$Version2)"/><xsl:text> </xsl:text><xsl:text>
</xsl:text>
	</xsl:for-each>
	<xsl:for-each select="Software/Software/Software">
		<xsl:variable name="Name3" select="@Name"/>
		<xsl:variable name="Version3" select="@Version"/>
		<xsl:text>#    </xsl:text><xsl:value-of select="concat($Name3,'-',$Version3)"/><xsl:text> </xsl:text><xsl:text>
</xsl:text>
	</xsl:for-each>


</xsl:template>

<!-- FlowcellInfo ==================================================================== -->
<xsl:template match="FlowcellInfo">

	<xsl:variable name="runID" select="@ID"/>


<xsl:text>
RUNS=</xsl:text><xsl:value-of select="$runID"/>

	<xsl:text>
RUNS_SAMPLES=</xsl:text>
	<xsl:for-each select="Lane/Sample[generate-id() = generate-id(key('samples',@SampleId)) and @Ref!='unknown']">
		<xsl:variable name="sampleID" select="@SampleId"/>
		<xsl:variable name="projectID" select="@ProjectId"/>
		<xsl:for-each select="../../Lane/Sample[@SampleId = $sampleID ]">
			<xsl:value-of select="concat($runID,':',$sampleID,':',$projectID)"/><xsl:text> </xsl:text>
		</xsl:for-each>
	</xsl:for-each>

	<xsl:text>

</xsl:text>

</xsl:template>

</xsl:stylesheet>







