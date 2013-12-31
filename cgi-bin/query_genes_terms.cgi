#!/usr/bin/python

import xml.etree.ElementTree as ET
import urllib2
import cgitb
import time
import cgi
import sys
import codecs

################################### Functions #################################
def xstr(s):
    if s is None:
        return ''
    return str(s)

def get_abstract_text(text_elements):
    str_list = []
    for element in text_elements:
        str_list.append(element.text)
    return '<br/>'.join(str_list)

def get_term_count(term_list, text):
    term_count = dict()
    for term in term_list:
        term_count[term] = text.count(term)
    return term_count

##################################### Handle ##################################
sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
cgitb.enable()

# get form data
form = cgi.FieldStorage()

# get gene ids and construct gene ids query string
gene_ids = form['genes'].value
gene_ids_list = gene_ids.split()
sorted_gene_ids_list = sorted(gene_ids_list, key=int)
gene_ids_qstr = 'id=' + '&id='.join(sorted_gene_ids_list)
gene_ids_summary_qstr = 'id=' + ','.join(sorted_gene_ids_list)

# get terms and construct terms query string
terms = form['terms'].value
terms_list = terms.split()
sorted_terms_list = sorted(terms_list)
terms_qstr = 'term=' + '+OR+'.join(terms_list)

# eutils base url
eutils_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils'
esummary_url = eutils_url + '/esummary.fcgi'
elink_url = eutils_url + '/elink.fcgi'
efetch_url = eutils_url + '/efetch.fcgi'

# pubmed url
pubmed_url = 'http://www.ncbi.nlm.nih.gov/pubmed'

# get gene summary information
genes_summary = urllib2.urlopen(esummary_url+'?db=gene&'+gene_ids_summary_qstr)
genes_summary_str = genes_summary.read()

# parse gene summary information xml
genes_summary_root = ET.fromstring(genes_summary_str)

# get abstract information
abstract_fixed_param = '?dbfrom=gene&db=pubmed&linkname=gene_pubmed&cmd=neighbor_history'
abstract_elink_qstr = elink_url+abstract_fixed_param+'&'+gene_ids_qstr+'&'+terms_qstr
abstract_elink = urllib2.urlopen(abstract_elink_qstr)
abstract_elink_str = abstract_elink.read()

# parse abstract elink information xml
abstract_elink_root = ET.fromstring(abstract_elink_str)
geneid_webenv_key = dict()
for i in xrange(len(abstract_elink_root)):
    gene_id = abstract_elink_root[i].find('.//Id').text
    query_key = abstract_elink_root[i].find('.//QueryKey').text
    webenv = abstract_elink_root[i].find('./WebEnv').text
    geneid_webenv_key[gene_id] = (webenv, query_key)

# fixed params for efetch
efetch_fixed_param = '?db=pubmed&rettype=abstract&retmode=xml'

##################################### HTML ####################################

# header
print 'Content-type:text/html\r\n\r\n'
print '<html>'

# head
print '<head>'
print '<title>Search Summary</title>'
print """<link rel="stylesheet" href="../tablesorter/themes/blue/style.css" """
print """ type="text/css" media="print, projection, screen" />"""
print """
        <style>
            table {
                border-bottom: 1px Solid Black;         
                border-right: 1px Solid Black;         
                border-collapse : collapse;  
            }
            table td, table th {    
                border-left: 1px Solid Black;         
                border-top: 1px Solid Black;              
                border-bottom:none;    
                border-right:none;  
            }
            table td {
                max-width:500px;
                word-wrap:break-word;
            }
            table.tablesorter {
                font-size: 14px; 
            }
        </style>
        """
print '<script src=\"../tablesorter/jquery-latest.js\"></script>'
print '<script src=\"../tablesorter/jquery.tablesorter.js\"></script>'
print '<script type="text/javascript">'
print '$(document).ready(function() {'
for gene_id in sorted_gene_ids_list:
    article_tb_id = '#articles_gene_id_'+str(gene_id)
    print '$(\"'+article_tb_id+'\").tablesorter({'
    print 'widgets: [\'blue\'], '
    print ' sortInitialOrder: \'desc\','
    print '});'  
print '});' 
print '</script>'
print '</head>'

# body
print '<body>'
print '<a id=\"top\">Search Summary</a><br/><br/>'

# create navigation
print '<a id=\"nav\">Navigation by Gene ID</a><br/>'
for gene_id in sorted_gene_ids_list:
    print '<a href=\"#gene_id_'+gene_id+'\">'+gene_id+'</a>'
print """<br/><a href="/">Make Another Search</a><br/><br/><hr/>"""

# iterate though gene ids
for i in xrange(len(sorted_gene_ids_list)):
    # parse gene summary xml
    gene_info = genes_summary_root[i]
    gene_id = gene_info.find('Id').text
    gene_name = gene_info.find("./Item[@Name='Name']").text
    gene_descp = gene_info.find("./Item[@Name='Description']").text
    gene_aliases = gene_info.find("./Item[@Name='OtherAliases']").text
    gene_desg = gene_info.find("./Item[@Name='OtherDesignations']").text
    gene_ns = gene_info.find("./Item[@Name='NomenclatureSymbol']").text
    gene_nn = gene_info.find("./Item[@Name='NomenclatureName']").text
    gene_summary = gene_info.find("./Item[@Name='Summary']").text
    gene_url = 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=Retrieve&dopt=full_report&list_uids='+gene_id
    rna_url = 'http://genomernai.de/GenomeRNAi/genedetails/'+gene_id
    biogps_url = 'http://biogps.org/#goto=genereport&id='+gene_id
    genecards_url = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene='+gene_name
    
    # create bookmark
    print '<a id=\"gene_id_'+gene_id+'\" href=\"#top\">'
    print 'Return to Top</a><br/><br/>'
    
    # print gene summary information table
    print 'Gene Summary</br>'
    print '<table>'
    print '<tr><td>Id</td><td>'+xstr(gene_id)+'</td></tr>'
    print '<tr><td>Name</td><td>'+xstr(gene_name)+'</td></tr>'
    print '<tr><td>Description</td><td>'+xstr(gene_descp)+'</td></tr>'
    print '<tr><td>Aliases</td><td>'+xstr(gene_aliases)+'</td></tr>'
    print '<tr><td>Other Designations</td><td>'+xstr(gene_desg)+'</td></tr>'
    print '<tr><td>Summary</td><td>'+xstr(gene_summary)+'</td></tr>'
    print '<tr><td>Gene URL</td><td><a href=\"'+gene_url+'\">click here</a></td></tr>'
    print '<tr><td>RNA URL</td><td><a href=\"'+rna_url+'\">click here</a></td></tr>'
    print '<tr><td>BioGPS URL</td><td><a href=\"'+biogps_url+'\">click here</a></td></tr>'
    print '<tr><td>GeneCards URL</td><td><a href=\"'+genecards_url+'\">click here</a></td></tr>'
    print '</table>'
    print '<br/>'
    
    # print article summary table
    print 'Article Summary</br>'
    print '<table id=\"articles_gene_id_'+xstr(gene_id)+'\"class=\"tablesorter\">'
    
    # print table header
    print '<thead>'
    num_terms = len(terms_list)
    print '<tr>'
    print '<td>Title & Abstract</td>'
    for term in sorted_terms_list:
        print '<th>'+term+'</th>'
    print '<tr/>'
    print '</thead>'
    
    # get article abstract from eutil server
    webenv = geneid_webenv_key[gene_id][0]
    query_key = geneid_webenv_key[gene_id][1]
    abstract_efetch_qstr = efetch_url+efetch_fixed_param+'&WebEnv='+webenv+'&query_key='+query_key
    abstract_efetch = urllib2.urlopen(abstract_efetch_qstr)
    abstract_efetch_str = abstract_efetch.read()
    abstract_efetch_root = ET.fromstring(abstract_efetch_str)
    num_articles = len(abstract_efetch_root)
    
    # print table body body
    print '<tbody>'
    for j in xrange(num_articles):
        pmid = abstract_efetch_root[j].find('.//PMID').text
        title = abstract_efetch_root[j].find('.//ArticleTitle').text
        text_elements = abstract_efetch_root[j].findall('.//AbstractText')
        text = get_abstract_text(text_elements)
        print '<tr>'
        print '<td><a href=\"'+pubmed_url+'/'+pmid+'\">'+title+'</a><br/>'+text+'</td>'
        term_count = get_term_count(terms_list, text)
        for term in sorted_terms_list:
            print '<td>'+str(term_count[term])+'</td>'
        print '</tr>'
    print '</tbody>'
    
    # end article summary table
    print '</table>'
    
    # print separator
    print '<hr/>'
    
    # sleep for 1 second
    time.sleep(1)
    
# end html
print '</body>'
print '</html>'