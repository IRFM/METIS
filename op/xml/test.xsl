<?xml version="1.0"?>
<?xml-stylesheet type="text/xsl" href="C:\Documents and Settings\Administrator\Desktop\dina2htm.xsl"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="html" version="1.0" encoding="UTF-8" indent="yes"/>

<!-- dina2htm.xls converts all xml files into html files at different stages -->

<xsl:template match="/*">
		<html>
		<head>
		</head>
		<body>
				
<xsl:apply-templates select="general" />
			
</body>
</html>
</xsl:template>
	
					
<xsl:template match="general">
<xsl:variable name="temp">t<xsl:value-of select="shortname"/>t</xsl:variable>
<xsl:element name="img"><xsl:attribute name="src"><xsl:value-of select="normalize-space(shortname)"/>2xml.gif</xsl:attribute><xsl:attribute name="align">center</xsl:attribute>
<xsl:attribute name="width">560</xsl:attribute>
<xsl:attribute name="height">420</xsl:attribute>
</xsl:element>
</xsl:template>

</xsl:stylesheet>
