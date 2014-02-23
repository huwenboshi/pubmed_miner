#!/usr/bin/python -u

import xml.etree.ElementTree as ET
import urllib2
import cgitb
import time
import cgi
import sys
import codecs
import os

from consts import *
from db_utils import *
from filter_utils import *
from utils import *

sys.path.append('./httpagentparser-1.5.1')
import httpagentparser

sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
cgitb.enable()

con = connect_db()
con_lite = connect_db_lite(db_loc)

########################### HTTP HTML HEADER ###################################

# print http_header and html header
user_agent = os.environ.get("HTTP_USER_AGENT", "N/A")
os_browser = httpagentparser.simple_detect(user_agent)
browser = os_browser[1]
if(browser.lower().find('safari') >= 0):
    print http_header_mac
else:
    print http_header
print """
<!DOCTYPE html>
<html>
    <head>
        <title>
            PubMed Miner - Gene Filter Result
        </title>
        <style type="text/css">

            #user_input_form {
                display: none;
            }
        </style>
        <link rel="stylesheet" 
            href="../../huwenbo/table_sorter/themes/blue/style.css" />
        <script type="text/javascript"
                src="../../huwenbo/table_sorter/jquery-latest.js">
        </script>
        <script type="text/javascript"
                src="../../huwenbo/table_sorter/jquery.tablesorter.js">
        </script>
        <script type="text/javascript">
            $(document).ready(function(){
            
                ///////// table sorter for ewas gene /////////
                $.tablesorter.addParser({
                    id: 'ewas_gene_exp_pval',
                    is: function(s) {
                        return false;
                    },
                    format: function(text, table, cell) {
                        return parseFloat(text);
                    },
                    type: 'numeric'
                });
                $('#ewas_gene_exp_tbl').tablesorter({
                    headers: {
                        4: {
                            sorter: 'ewas_gene_exp_pval'
                        }
                    }
                });
                $('#gwas_gene_exp_tbl').tablesorter({
                    headers: {
                        4: {
                            sorter: 'ewas_gene_exp_pval'
                        }
                    }
                });
                
                ///////// table sorter for ewas prot /////////
                $.tablesorter.addParser({
                    id: 'ewas_prot_exp_pval',
                    is: function(s) {
                        return false;
                    },
                    format: function(text, table, cell) {
                        return parseFloat(text);
                    },
                    type: 'numeric'
                });
                $('#ewas_prot_exp_tbl').tablesorter({
                    headers: {
                        4: {
                            sorter: 'ewas_prot_exp_pval'
                        }
                    }
                });
                $('#gwas_prot_exp_tbl').tablesorter({
                    headers: {
                        4: {
                            sorter: 'ewas_prot_exp_pval'
                        }
                    }
                });
                
                ///////// table sorter for ewas trait /////////
                $.tablesorter.addParser({
                    id: 'ewas_trait_pval',
                    is: function(s) {
                        return false;
                    },
                    format: function(text, table, cell) {
                        return parseFloat(text);
                    },
                    type: 'numeric'
                });
                $('#ewas_trait_tbl').tablesorter({
                    headers: {
                        6: {
                            sorter: 'ewas_trait_pval'
                        }
                    }
                });
                
                ////////// table sorter for imp_summary_tbl /////////
                $('#imp_summary_tbl').tablesorter();
                
                ///////// submit button /////////
                $('#confirm_btn').live('click', function() {
                    $("#user_input_form").submit();
                });
                
                ///////// show hide tables ///////
                $('.show_hide').live('click', function() {
                if ($(this).html() == "hide") {
                    $(this).siblings("table").css({display:'none'});
                    $(this).html("show");
                }
                else {
                    $(this).siblings("table").css({display:'inline-table'});
                    $(this).html("hide");
                }
            });
            });
        </script>
    </head>
    <body>
        <h2>
            Filter Result
        </h2>
        <a href="../../huwenbo/index.html">Make Another Search</a>
        <br/>
        <hr/>
"""

################################# Get User Input ###############################

# get form data
form = cgi.FieldStorage()

# get selected implication types
imp_types = form.getlist('imp_type')

# get logic relation between implication types
imp_type_logic_sel = form.getvalue("imp_type_logic_sel")

# distance threshold between mCG site and implicated gene
ewas_gene_exp_max_distance = form.getvalue("ewas_gene_exp_max_distance")
ewas_prot_exp_max_distance = form.getvalue("ewas_prot_exp_max_distance")
ewas_trait_max_distance = form.getvalue("ewas_trait_max_distance")

# distance threshold between snp site and implicated gene
gwas_gene_exp_max_distance = form.getvalue("gwas_gene_exp_max_distance")
gwas_prot_exp_max_distance = form.getvalue("gwas_prot_exp_max_distance")
gwas_trait_max_distance = form.getvalue("gwas_trait_max_distance")

# intersect of union the genes found by different associations ewas
ewas_assoc_logic_sel = form.getvalue("ewas_assoc_logic_sel")

# intersect of union the genes found by different associations gwas
gwas_assoc_logic_sel = form.getvalue("gwas_assoc_logic_sel")

# get associations threshold ewas
ewas_gene_exp_pval = form.getvalue("ewas_gene_exp_pval")
ewas_prot_exp_pval = form.getvalue("ewas_prot_exp_pval")
ewas_trait_pval = form.getvalue("ewas_trait_pval")
ewas_trait_names = form.getlist("ewas_trait_names")
ewas_tables = form.getlist("ewas_tables")

gwas_gene_exp_pval = form.getvalue("gwas_gene_exp_pval")
gwas_prot_exp_pval = form.getvalue("gwas_prot_exp_pval")
gwas_trait_pval = form.getvalue("gwas_trait_pval")
gwas_trait_names = form.getlist("gwas_trait_names")
gwas_tables = form.getlist("gwas_tables")

# get user additional input
id_type = form.getvalue("id_type")
user_genes = form.getvalue("user_genes")
curated_terms = form.getlist("curated_terms")
user_terms = form.getvalue("user_terms")
expand_term = form.getvalue("expand_term")
search_scope = form.getvalue("search_scope")

############################## Cet Query Result ################################

# get query result from database
ewas_gwas_result = handle_query(con,
                            imp_types,
                            imp_type_logic_sel,
                            ewas_tables,
                            ewas_gene_exp_pval, 
                            ewas_prot_exp_pval,
                            ewas_trait_pval,
                            ewas_gene_exp_max_distance,
                            ewas_prot_exp_max_distance,
                            ewas_trait_max_distance,
                            ewas_trait_names,
                            ewas_assoc_logic_sel,
                            gwas_tables,
                            gwas_gene_exp_pval, 
                            gwas_prot_exp_pval,
                            gwas_trait_pval,
                            gwas_gene_exp_max_distance,
                            gwas_prot_exp_max_distance,
                            gwas_trait_max_distance,
                            gwas_trait_names,
                            gwas_assoc_logic_sel)

ewas_query_result = ewas_gwas_result['ewas_query_result']
gwas_query_result = ewas_gwas_result['gwas_query_result']
combined_entrez_id_set = ewas_gwas_result['human_gene_set']

# get gene supporting loci information ewas
ewas_gene_support_info = get_ewas_gene_supporting_info(ewas_query_result)

# get gene supporting loci information gwas
gwas_gene_support_info = get_gwas_gene_supporting_info(gwas_query_result)

# get user genes
gene_list = combined_entrez_id_set
if(user_genes != None):
    user_genes_list = user_genes.split()
    # conver user gene to entrez id if it's gene symbol
    if(id_type == 'GENE_SYM'):
        sym_id_map = symbol2entrez(user_genes_list)
        user_genes_list = list(sym_id_map.values())
    for gene in user_genes_list:
        if(gene.isdigit()):
            gene_list.add(int(gene))
gene_list = list(gene_list)

# get curated terms
all_terms_list = set()
if(curated_terms != None):
    for term in curated_terms:
        all_terms_list.add(term.strip())

# get user terms
if(user_terms != None):
    user_terms_list = user_terms.split('\n')
    for term in user_terms_list:
        term = term.strip()
        all_terms_list.add(term)
all_terms_list = list(all_terms_list)

# check user search scope
if(search_scope == 'TIAB'):
    search_scope = 'Title and abstract'
else:
    search_scope = 'Full-text'

############################### Print Result ###################################

# print query summary
print """
    <h3>
        Query summary
    </h3>
"""

print '<table id="summary_tbl">'
print '<tr>'
print '<td><b>Number of genes</b></td><td>%d</td>' % len(gene_list)
print '</tr>'
print '<tr>'
print '<td><b>Number of terms</b></td><td>%d</td>' % len(all_terms_list)
print '</tr>'
print '<tr>'
print '<td><b>Search scope</b></td><td>%s</td>' % search_scope
print '</tr>'
print '</table>'
print '<br/>'
print '<button id="confirm_btn">search</button>'
print '<hr/>'

# print gene implication summary
print """
        <div>
            <b>Gene implication summary</b>
            <button class="show_hide" type="button">hide</button>
            <br/>
            
      """
print '<table id="imp_summary_tbl" class="tablesorter">'
print '<thead>'
print """<tr>
            <th>Human gene Entrez ID</th>
            <th>Number of EWAS loci<br/>associated with gene expression</th>
            <th>Number of EWAS loci<br/>associated with protein expression</th>
            <th>Number of EWAS loci<br/>associated with phenotypes</th>
            <th>Number of GWAS loci<br/>associated with gene expression</th>
            <th>Number of GWAS loci<br/>associated with protein expression</th>
            <th>Total number of EWAS loci </th>
            <th>Total number of GWAS loci </th>
            <th>Total number of loci </th>
         </tr>"""
print '</thead>'
print '<tbody>'

for key in combined_entrez_id_set:
    key = key[0]
    ewas_assoc_pos = dict()
    if(key in ewas_gene_support_info):
        ewas_assoc_pos = ewas_gene_support_info[key]
    num_ewas_gene_exp = safe_len(safe_getval(ewas_assoc_pos, 'gene_exp'))
    num_ewas_prot_exp = safe_len(safe_getval(ewas_assoc_pos, 'prot_exp'))
    num_ewas_trait = safe_len(safe_getval(ewas_assoc_pos, 'trait'))
    num_ewas_tot = num_ewas_gene_exp + num_ewas_prot_exp + num_ewas_trait
    
    gwas_assoc_pos = dict()
    if(key in gwas_gene_support_info):
        gwas_assoc_pos = gwas_gene_support_info[key]
    num_gwas_gene_exp = safe_len(safe_getval(gwas_assoc_pos, 'gene_exp'))
    num_gwas_prot_exp = safe_len(safe_getval(gwas_assoc_pos, 'prot_exp'))
    num_gwas_tot = num_gwas_gene_exp + num_gwas_prot_exp
    
    num_tot = num_ewas_tot + num_gwas_tot
    
    print '<tr>'
    print '<td>%s</td>' % key
    print '<td>%d</td>' % num_ewas_gene_exp
    print '<td>%d</td>' % num_ewas_prot_exp
    print '<td>%d</td>' % num_ewas_trait
    print '<td>%d</td>' % num_gwas_gene_exp
    print '<td>%d</td>' % num_gwas_prot_exp
    print '<td>%d</td>' % num_ewas_tot
    print '<td>%d</td>' % num_gwas_tot
    print '<td>%d</td>' % num_tot
    
    print '</tr>'
    
print '</tbody>'
print '</table>'
print '</div>'
print '<hr/>'

print_ewas_query_result(ewas_query_result)
print '<hr/>'
print_gwas_query_result(gwas_query_result)

################################## Hidden Form #################################

genes_form_str = ''
for i in xrange(len(gene_list)):
    genes_form_str += str(gene_list[i])
    if(i < len(gene_list)-1):
        genes_form_str += '\n'

terms_form_str =''
for i in xrange(len(all_terms_list)):
    terms_form_str += all_terms_list[i].strip()
    if(i < len(all_terms_list)-1):
            terms_form_str += '\n'

print """
        <form id="user_input_form" name="user_input" target="_blank" 
            action="./search.cgi" method="post" style="display:none">
                
            <input type="radio" id="gene_entrez" name="id_type" 
                value="entrez_id" checked />Entrez ID
                   
            <input type="radio" id="gene_sym" name="id_type" 
                value="gene_sym" />Gene symbol
            
            <textarea id="gene_txt" name="genes" 
                rows="5" cols="20" >%s</textarea>
            
            <input type="checkbox" name="expand_term" />
            
            <textarea id="terms_txt" name="terms" 
                rows="5" cols="20">%s</textarea>
            
            <input type="checkbox" name="tiab_only" checked="checked"/> 
        </form>
""" % (genes_form_str, terms_form_str)


################################### Ending #####################################

print """
    </body>
</html>
"""
con.close()
con_lite.close()
