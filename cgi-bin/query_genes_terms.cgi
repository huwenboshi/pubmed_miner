#!/usr/bin/python

import urllib2
import cgitb
import cgi
import sys

cgitb.enable()

# get form data
form = cgi.FieldStorage()

# get gene ids and construct gene ids query string
gene_ids = form['genes'].value
gene_ids_list = gene_ids.split()
gene_ids_qstr = 'id=' + '&id='.join(gene_ids_list)
gene_ids_summary_qstr = 'id=' + ','.join(gene_ids_list)

# get terms and construct terms query string
terms = form['terms'].value
terms_list = terms.split()
terms_qstr = 'term=' + '+OR+'.join(terms_list)

# eutils base url
eutils_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils'
esummary_url = eutils_url + '/esummary.fcgi'

# get gene summary information
genes_summary = urllib2.urlopen(esummary_url+'?db=gene&'+gene_ids_summary_qstr)
genes_summary_str = genes_summary.read()

##################################### HTML #####################################

# header
print 'Content-type:text/html\r\n\r\n'
print '<html>'

# head
print '<head>'
print '<title>Search Summary</title>'
print '</head>'

# body
print '<body>'
print '<h2>Search Summary</h2>'
print '<a>' + gene_ids_summary_qstr + '</a><br/>'
print '<a>' + gene_ids_qstr + '</a><br/>'
print '<a>' + terms_qstr + '</a><br/>'
print '</body>'

print '</html>'