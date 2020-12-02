<?xml version="1.0"?>
<!-- greeting.xsl -->
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
	<xsl:output method="html"/>
	<xsl:template match="surveillance_report">
		<html>
            <head>
                <style>
                    table, th, td {
                        border: 1px solid black;
                        border-collapse: collapse;
                        padding: 5px
                    }
                    .id_content {
                        max-width: 800px;
                        margin: auto;
                    }
                    .column {
                        float: left;
                        width: 20%;
                        padding: 5px;
                        margin-left:0 auto;
                        margin-right:0 auto;
                    }
                    .row:after {
                        content: "";
                        display: table;
                        clear: both;
                    }
                </style>
            </head>
			<body>
                <xsl:apply-templates select="title"/>
                <div class="id_content">
                    <div class="row">
                        <div class="column">
                            <b>Sample ID:</b>
                        </div>
                        <div class="column">
                            <xsl:apply-templates select="sample_id"/>
                        </div>
                        <div class="column">
                            <b>Date/Time:</b>
                        </div>
                        <div class="column">
                            <xsl:apply-templates select="date"/>
                        </div>
                    </div>
                    <div class="row">
                        <div class="column">
                            <b>Pipeline:</b>
                        </div>
                        <div class="column">
                            <xsl:apply-templates select="pipeline"/>
                        </div>
                        <div class="column">
                            <b>Lineage:</b>
                        </div>
                        <div class="column">
                            <xsl:apply-templates select="lineage"/>
                        </div>
                    </div>
                </div>
                <h2>Coverage</h2>
                <table style="100%; margin: 0px auto;">
                    <tr>
                        <th>Gene</th>
                        <th>Coverage Depth</th>
                        <th>Coverage Percent</th>
                    </tr>
                    <xsl:apply-templates select="coverage/gene"/>
                </table>
                <h2>Deletions</h2>
                <xsl:choose>
                    <xsl:when test="count(deletions/loci) > 0">
                        <table style="100%; margin: 0px auto;">
                            <tr>
                                <th>Gene</th>
                                <th>Deletion Type</th>
                            </tr>
                            <xsl:apply-templates select="deletions/loci"/>
                        </table>
                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:text>None Found</xsl:text>
                    </xsl:otherwise>
                </xsl:choose>
                <h2>SNPs</h2>
                    <xsl:apply-templates select="mutations/snp"/>
            </body>
        </html>
    </xsl:template>

    <xsl:template match="title">
        <h1 align="center">
            <xsl:value-of select="."/>
        </h1>
    </xsl:template>

    <xsl:template match="sample_id">
            <xsl:value-of select="."/>
    </xsl:template>

    <xsl:template match="date">
        <xsl:value-of select="substring-before(.,'.')"/>
    </xsl:template>

    <xsl:template match="pipeline">
        <xsl:value-of select="name"/>
        <xsl:text> </xsl:text>
        <xsl:value-of select="version"/>
    </xsl:template>

    <xsl:template match="lineage">
        <xsl:value-of select="code"/>
        <xsl:text>(</xsl:text>
        <xsl:value-of select="name"/>
        <xsl:text>)</xsl:text>
    </xsl:template>

    <xsl:template match="gene">
        <tr>
            <td><xsl:value-of select="@name"/></td>
            <td><xsl:value-of select="depth"/></td>
            <td><xsl:value-of select="percent"/></td>
        </tr>
    </xsl:template>

    <xsl:template match="loci">
        <tr>
            <td><xsl:value-of select="@name"/></td>
            <td><xsl:value-of select="type"/></td>
        </tr>
    </xsl:template>

    <xsl:template match="snp">
        <b><xsl:value-of select="@gene"/></b>
        <pre>&#xA;</pre>
        <xsl:apply-templates select="nuchange"/>
    </xsl:template>

    <xsl:template match="nuchange">
        <xsl:text>&#x9;</xsl:text>
        <xsl:text>Nucleotide Change: </xsl:text>
        <xsl:value-of select="@name"/>
        <xsl:text>&#x9;</xsl:text>
        <xsl:text>Amino Acid Change: </xsl:text>
        <xsl:value-of select="@aachange"/>
        <pre>&#xA;</pre>
        <xsl:apply-templates select="resistance"/>
    </xsl:template>

    <xsl:template match="resistance">
        <xsl:text>&#x9;&#9;Drug: </xsl:text>
        <xsl:value-of select="drug"/>
        <xsl:if test="string-length(drug) &lt; 10">
            <xsl:text>&#x9;</xsl:text>
        </xsl:if>
        <xsl:text>&#x9;</xsl:text>
        <xsl:text>Confidence: </xsl:text>
        <xsl:value-of select="confidence"/>
        <pre>&#xA;</pre>
    </xsl:template>

</xsl:stylesheet>