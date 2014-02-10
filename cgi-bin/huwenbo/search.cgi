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
from consts import *
from genedb import *

sys.path.append('./httpagentparser-1.5.1')
import httpagentparser

sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
cgitb.enable()

con = connect_db('/home/huwenbo/pubmed_miner_db/pubmed_miner.db')

########################### HTTP HTML HEADER ###################################

# print http_header and html header
user_agent = os.environ.get("HTTP_USER_AGENT", "N/A")
os_browser = httpagentparser.simple_detect(user_agent)
browser = os_browser[1]
if(browser.lower().find('safari') >= 0):
    print http_header_mac
else:
    print http_header
print html_header

################################# Get User Input ###############################

# get form data
form = cgi.FieldStorage()

# get implication type and logic relation between selections
imp_types = form.getlist("imp_type")
imp_type_logic_sel = form.getvalue("imp_type_logic_sel")

max_distance = form.getvalue("max_distance")

assoc_logic_sel = form.getvalue("assoc_logic_sel")

gene_exp_pval = form.getvalue("gene_exp_pval")
protein_exp_pval = form.getvalue("protein_exp_pval")
trait_pval = form.getvalue("trait_pval")
trait_names = form.getlist("trait_names")

id_type = form.getvalue("id_type")
user_genes = form.getvalue("user_genes")
user_terms = form.getvalue("user_terms")
expand_term = form.getvalue("expand_term")
search_scope = form.getvalue("search_scope")

############################## Construct Query #################################
"""
query_gene = make_general_ewas_query('liver_expression_ewas',
    gene_exp_pval, max_distance)
print query_gene
print '<br/>'
query_prot = make_general_ewas_query('liver_proteomics_ewas',
    protein_exp_pval, max_distance)
print query_prot
print '<br/>'
query_trait = make_clinical_trait_ewas_query('clinical_metabolite_traits_ewas',
    trait_pval, max_distance, trait_names)
print query_trait
print '<br/>'
query_probe_to_entrez = make_to_human_gene_id_query('mouse_probe_mouse_mgi',
    'mouse_probe_id')
print query_probe_to_entrez
print '<br/>'
query_entrez_to_entrez = make_to_human_gene_id_query('mouse_entrez_mouse_mgi',
    'mouse_entrez_id')
print query_entrez_to_entrez
print '<br/>'
query_trans_to_entrez = make_to_human_gene_id_query('mouse_trans_mouse_mgi',
    'mouse_transcript_id')
print query_trans_to_entrez
print '<br/>'
print imp_type_logic_sel
print '<br/>'
print assoc_logic_sel
"""
ewas_query_result = get_ewas_query_result(gene_exp_pval, protein_exp_pval,
    trait_pval, max_distance, trait_names, con, assoc_logic_sel)
print ewas_query_result
