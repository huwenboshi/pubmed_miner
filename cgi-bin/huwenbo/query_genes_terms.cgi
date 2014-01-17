#!/usr/bin/python -u

import xml.etree.ElementTree as ET
import urllib2
import cgitb
import time
import cgi
import sys
import codecs

from utils import *
from consts import *

sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
cgitb.enable()

################################# Get User Input ###############################

# get form data
form = cgi.FieldStorage()

id_type = form['id_type'].value

# get entrez ids, convert to entrez gene id if input is gene symbol
gene_ids_list = []
if(id_type == 'gene_sym'):
    gene_symbols = form['genes'].value
    gene_symbols_list = gene_symbols.split()
    sym_id = symbol2entrez(gene_symbols_list)
    for key in sym_id:
        gene_ids_list.append(sym_id[key])
    gene_ids_list = sorted(gene_ids_list, key=int)
else:
    gene_ids = form['genes'].value
    gene_ids_list = sorted(gene_ids.split(), key=int)

terms = ''
tiab_only = True
if('terms' in form):
    terms = form['terms'].value
terms_list = sorted(terms.split())
if('tiab_only' not in form):
    tiab_only = False

############################## Download Gene Info ##############################

# download gene summary info
genes_info_list = get_gene_info(gene_ids_list)

# get all gene related article count
genes_pmids_cnt = get_genes_pmid_count(gene_ids_list)

# get WebEnv and QueryKey for searching within TIAB
gene_webenv_querykey = None
if(tiab_only):
    gene_webenv_querykey = get_webenv_querykey(gene_ids_list, terms_list)

# sleep for 1 second to obey the 3 queries/sec rule
time.sleep(1)

##################################### HTML #####################################

# print http and html header
print http_header
print html_header

# print body start, navigation bar
print """
<body>
    <h2><a id="top">Search Result</a></h2>
    <a href="../../huwenbo/index.html">Make Another Search</a><br/><br/>
"""
print '<a id="nav"><b>Navigation by Gene ID</b></a><br/>'
for gene_id in gene_ids_list:
    print '<a href="#gene_id_%s">%s</a>' % (gene_id, gene_id)
print '<br/><br/><hr/>'

sys.stdout.flush()

# create overview
print """
<div>
    <a href="#top">Return to Top</a><br/><br/>
    <b>Overview (choose terms to sort by the sum of number of abstracts
        containing these terms)
    </b>
    <button class="show_hide" type="button">hide</button><br/>
    <a>each cell shows the number of abstracts related to the gene 
        (left-most column) and term (top-most row)
    </a>
    <br/>
    <br/>
    <a id="loading">Loading...</a>
    <table id="overview_top">
    </table>
</div>
<br/>
<hr/>
"""

# initialize gene term count
gene_term_count = init_gene_term_cnt(gene_ids_list, terms_list)

# iterate though gene ids
for i in xrange(len(gene_ids_list)):

    # get gene id
    gene_id = gene_ids_list[i]

    # create "return to top" bookmark and title
    print '<a id="gene_id_%s" href="#top">Return to Top</a><br/><br/>'%gene_id
    print '<b>Info and Related Abstracts for Gene %s</b><br/><br/>'%gene_id
    
    # print gene info table
    print """
        <div>
            <a><b>Gene Info</b></a>
            <button class="show_hide" type="button">hide</button>
            <br/>
            
    """
    print_gene_info(genes_info_list[i], genes_pmids_cnt[i])
    print """
        </div>
        <br/>
    """
     
    # search result for tiab search
    if(tiab_only):
        gene_term_count = print_tiab_search_result(terms_list, gene_id,
            gene_webenv_querykey, gene_term_count)
        time.sleep(1)
    else:
        print_fulltext_search_result(terms_list, gene_id, gene_term_count)
    

# print overview bottom, this is a dummy and will not be displayed
print '<table id="overview_bottom">'
print '<tr>'
print '<td>Gene ID\Term</td>'
# print first row
for term in terms_list:
    print '<td><input type="checkbox" id="overview_%s" class="overview_opt">%s</td>' % (term, term)
print '</tr>'
# print following rows
for gene_id in gene_ids_list:
    print '<tr>'
    print '<td><a href="#gene_id_%s">%s</a></td>' % (gene_id, gene_id)
    for term in terms_list:
        print '<td><a class="abstract_count">%d</a></td>' % gene_term_count[gene_id][term]
    print '</tr>'
print '</table>'
    
# print body end
print '</body>'

# print html end
print '</html>'
