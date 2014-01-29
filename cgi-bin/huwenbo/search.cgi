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
print 'implicatoin type<br/>'
imp_types = form.getlist("imp_type")
imp_type_logic_sel = form.getvalue("imp_type_logic_sel")
print imp_types
print imp_type_logic_sel
print '<br/>'

print 'gene expressoin<br/>'
logic_gene_exp = form.getvalue("logic_gene_exp")
gene_exp_pval = form.getvalue("gene_exp_pval")
print logic_gene_exp
print gene_exp_pval
print '<br/>'

print 'protein expressoin<br/>'
logic_gene_exp = form.getvalue("logic_protein_exp")
gene_exp_pval = form.getvalue("protein_exp_pval")
print logic_gene_exp
print gene_exp_pval
print '<br/>'

print 'distance<br/>'
logic_distance = form.getvalue("logic_distance")
max_distance = form.getvalue("max_distance")
print logic_distance
print max_distance
print '<br/>'

print 'trait<br/>'
logic_trait = form.getvalue("logic_trait")
trait_pval = form.getvalue("trait_pval")
trait_logic_sel = form.getvalue("trait_logic_sel")
trait_names = form.getlist("trait_names")
print logic_trait
print trait_pval
print trait_logic_sel
print trait_names
print '<br/>'

print 'id type<br/>'
id_type = form.getvalue("id_type")
print id_type
print '<br/>'

print 'user genes<br/>'
user_genes = form.getvalue("user_genes")
print user_genes
print '<br/>'

print 'user terms<br/>'
user_terms = form.getvalue("user_terms")
print user_terms
print '<br/>'

print 'expand term<br/>'
expand_term = form.getvalue("expand_term")
print expand_term
print '<br/>'

print 'search scope<br/>'
search_scope = form.getvalue("search_scope")
print search_scope
print '<br/>'
