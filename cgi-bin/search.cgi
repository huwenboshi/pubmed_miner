#!/usr/bin/python -u

import xml.etree.ElementTree as ET
import urllib2
import cgitb
import time
import cgi
import sys
import codecs
import os

from utils import *
from search_utils import *
from consts import *

sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
cgitb.enable()

con = connect_db()

########################### HTTP HTML HEADER ###################################

# print http_header and html header
print http_header
print html_header

################################# Get User Input ###############################

# get form data
form = cgi.FieldStorage()

id_type = form['id_type'].value

# id2sym, sym2id conversion
sym_id = dict()
id_sym = dict()

# get entrez ids/symbols, convert between symbols and ids
gene_ids_list = []
if(id_type == 'gene_sym'):
    gene_symbols = form['genes'].value
    gene_symbols_list = gene_symbols.split()
    sym_id = symbol2entrez(gene_symbols_list)
    for key in sym_id:
        gene_ids_list.append(sym_id[key])
        id_sym[sym_id[key]] = key
    gene_ids_list = sorted(gene_ids_list, key=int)
else:
    gene_ids = form['genes'].value
    gene_ids_list = sorted(gene_ids.split(), key=int)

terms = ''
tiab_only = True
if('terms' in form):
    terms = form['terms'].value
terms_list = sorted(list(set([term.strip() for term in terms.split('\n')])))
if('tiab_only' not in form):
    tiab_only = False

############################## Download Gene Info ##############################

# download gene summary info
genes_info_list = get_gene_info(gene_ids_list)

# parse gene id gene symbol conversion
for gene_info in genes_info_list:
    id_sym[gene_info['id']] = gene_info['name']

# get all gene related article count
genes_pmids_cnt = get_genes_pmid_count(gene_ids_list)

# get WebEnv and QueryKey for searching within TIAB
gene_webenv_querykey = None
if(tiab_only):
    gene_webenv_querykey = get_webenv_querykey(gene_ids_list, terms_list)

# sleep for 1 second to obey the 3 queries/sec rule
time.sleep(1)

########################### DISPLAY CONTENT ####################################

# print body start, navigation bar
print """
<body>
    <h2><a id="top">Search Result</a></h2>
    <a href="../index.html">Make Another Search</a><br/>
    <hr/>
"""

sys.stdout.flush()

# create overview
print """
<div>
    <b>Overview</b>
    <button class="show_hide" type="button">hide</button><br/>
    <a id="loading">Loading...</a>
    <table id="overview_top" class="heat-map">
    </table>
</div>
<hr/>
"""

# create network
print """
<div>
    <b>Term Network</b>
    <button class="show_hide" type="button">hide</button><br/>
    <a id="loading_gene_term_network">Loading...</a>
    <table id="gene_term_network_top">
    </table>
</div>
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
    
    # print gene gwas info table
    gwas_info_list = search_nhgri_gwas_catalog(con, id_sym[gene_id])
    print_nhgri_gwas_info_list(gwas_info_list)
    
    # print clinical trial table
    clinical_trial_result = search_clinical_trial(con, gene_id)
    print_clinical_trial_result(clinical_trial_result)
    
    # print gene cluster talbe
    gene_cluser_result = search_hypervariable_count(con, gene_id)
    print_hypervariable_count(gene_cluser_result)
     
    # search result for tiab search
    if(len(terms_list) > 0):
        if(tiab_only):
            gene_term_count = print_tiab_search_result(terms_list, gene_id,
                gene_webenv_querykey, gene_term_count, id_sym)
        else:
            print_fulltext_search_result(terms_list, gene_id,
                gene_term_count, id_sym)
    else:
        time.sleep(1)
        
    print '<hr/>'

###########################OVERVIEW & NETWORK ##################################

# print overview bottom, this is a dummy and will not be displayed
print """<table id="overview_bottom" cellpadding="0"
          cellspacing="0" border="0">
         <thead>
         <tr>"""
print '<td class="first">Term\Gene</t>'
for i in xrange(len(gene_ids_list)):
    gene_id = gene_ids_list[i]
    gene_sym = id_sym[gene_id]
    if(i == len(gene_ids_list)-1):
        print """<td class="last">
                 <a href="#gene_id_%s">%s(%s)</a>
                 </td>"""%(gene_id, gene_sym, gene_id)
    else:
        print '<td><a href="#gene_id_%s">%s(%s)</a></td>'%(gene_id,
            gene_sym, gene_id)
print """</tr>
         </thead>
         <tbody>"""
for term in terms_list:
    print '<tr class="stats-row">'
    print '<td class="stats-title"><input type="radio" name="sort_choice"'
    print ' class="overview_opt" id="overview_%s">%s</td>' % (term, term)
    for gene_id in gene_ids_list:
        print '<td><a class="abstract_count"'
        if(gene_term_count[gene_id][term]):
            print 'href="#%s"' % (term+"_"+gene_id)
        print '>%d</a></td>' % gene_term_count[gene_id][term]
    print '</tr>'
print """</tbody>
         </table>"""

# create the term network
print '<table id="gene_term_network_bottom">'
print '<tr><td>'
create_gene_term_network(gene_term_count, gene_ids_list, terms_list, id_sym)
print '</td></tr>'
print '</table>'

#################################### HTML END ##################################

# print body end
print '</body>'
print '</html>'

con.close()
