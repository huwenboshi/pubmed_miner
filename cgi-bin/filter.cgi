#!/usr/bin/python -u

import xml.etree.ElementTree as ET
import urllib2
import cgitb
import time
import cgi
import sys
import codecs
import os

"""
sys.path.insert(0,'/UCSC/Panel-Auxiliary/mysql_python/MySQL-python-1.2.3')
sys.path.insert(0,'/UCSC/Panel-Auxiliary/mysql_python/MySQL-python-1.2.3/build/lib.linux-x86_64-2.7')
sys.path.insert(0,'/UCSC/Panel-Auxiliary/mysql_python/MySQL-python-1.2.3/build/lib.linux-x86_64-2.7/_mysql.so')
"""

from consts import *
from filter_utils import *
from utils import *

sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
cgitb.enable()

con = connect_db()

########################### HTTP HTML HEADER ###################################

# print http_header and html header
print http_header
print '\n\n'

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
            href="../javascripts/table_sorter/themes/blue/style.css" />
        <link rel="stylesheet" type="text/css" href="../css/filter.css" />
        <script type="text/javascript"
                src="../javascripts/table_sorter/jquery-latest.js">
        </script>
        <script type="text/javascript"
                src="../javascripts/table_sorter/jquery.tablesorter.js">
        </script>
        <script type="text/javascript">
            $(document).ready(function(){
            
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
                
                $('#imp_summary_tbl').tablesorter();
                
                $('#confirm_btn').live('click', function() {
                    $("#user_input_form").submit();
                });
                
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
        <div align="left">
        <a href="../index.html">Make Another Search</a>
        </div>
        <div align="center">
        <h1>
            PubMed Miner - Gene Filter Result
        </h1>
        </div>
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

# get clinical trial threshold
gwas_gene_exp_cnt = form.getvalue("gwas_gene_exp_trial")
gwas_prot_exp_cnt = form.getvalue("gwas_prot_exp_trial")
gwas_trait_cnt = form.getvalue("gwas_trait_trial")
ewas_gene_exp_cnt = form.getvalue("ewas_gene_exp_trial")
ewas_prot_exp_cnt = form.getvalue("ewas_prot_exp_trial")
ewas_trait_cnt = form.getvalue("ewas_trait_trial")

# get gene trait assoc method
ewas_trait_dist_or_assoc = form.getvalue("ewas_trait_dist_or_assoc")

# get user additional input
id_type = form.getvalue("id_type")
user_genes = form.getvalue("user_genes")
curated_terms = form.getlist("curated_terms")
user_terms = form.getvalue("user_terms")
expand_term = form.getvalue("expand_term")
search_scope = form.getvalue("search_scope")

############################## Check User Input ################################

check_user_input(imp_types, imp_type_logic_sel, ewas_tables,
    ewas_gene_exp_pval, ewas_prot_exp_pval, ewas_trait_pval,
    ewas_gene_exp_max_distance, ewas_prot_exp_max_distance,
    ewas_trait_max_distance, ewas_gene_exp_cnt, ewas_prot_exp_cnt,
    ewas_trait_cnt, ewas_trait_names, ewas_assoc_logic_sel,
    gwas_tables, gwas_gene_exp_pval, gwas_prot_exp_pval, gwas_trait_pval,
    gwas_gene_exp_max_distance, gwas_prot_exp_max_distance,
    gwas_trait_max_distance, gwas_gene_exp_cnt, gwas_prot_exp_cnt,
    gwas_trait_cnt, gwas_trait_names, gwas_assoc_logic_sel)

############################## Cet Query Result ################################

# get query result from database
ewas_gwas_result = handle_query(con,imp_types, imp_type_logic_sel, ewas_tables,
    ewas_gene_exp_pval, ewas_prot_exp_pval, ewas_trait_pval,
    ewas_gene_exp_max_distance, ewas_prot_exp_max_distance,
    ewas_trait_max_distance, ewas_gene_exp_cnt, ewas_prot_exp_cnt,
    ewas_trait_cnt, ewas_trait_names, ewas_trait_dist_or_assoc,
    ewas_assoc_logic_sel,
    gwas_tables, gwas_gene_exp_pval, gwas_prot_exp_pval, gwas_trait_pval,
    gwas_gene_exp_max_distance, gwas_prot_exp_max_distance,
    gwas_trait_max_distance, gwas_gene_exp_cnt, gwas_prot_exp_cnt,
    gwas_trait_cnt, gwas_trait_names, gwas_assoc_logic_sel)
                            
ewas_query_result = ewas_gwas_result['ewas_query_result']
gwas_query_result = ewas_gwas_result['gwas_query_result']
combined_entrez_id_set = ewas_gwas_result['human_gene_set']

# count the number of implicating sites for each implicated gene
ewas_gene_support_info = count_implicating_sites_ewas(con, ewas_tables)
gwas_gene_support_info = count_implicating_sites_gwas(con, gwas_tables)
id_sym_trial = get_symbol_citeline_count(ewas_query_result, gwas_query_result)

############################ PROCESS USER ADD-ON ###############################

# get genes
gene_list = combined_entrez_id_set
if(user_genes != None):
    user_genes_list = user_genes.split()
    if(id_type == 'GENE_SYM'):
        sym_id_map = symbol2entrez(user_genes_list)
        user_genes_list = list(sym_id_map.values())
    for gene in user_genes_list:
        if(gene.isdigit()):
            gene_list.add(int(gene))
gene_list = list(gene_list)

# get terms
all_terms_list = set()
if('ewas_imp' in imp_types):
    ewas_trait_names = format_trait_names(ewas_trait_names)
    all_terms_list = all_terms_list.union(set(ewas_trait_names))
if('gwas_imp' in imp_types):
    gwas_trait_names = format_trait_names(gwas_trait_names)
    all_terms_list = all_terms_list.union(set(gwas_trait_names))
if(curated_terms != None):
    all_terms_list = all_terms_list.union(set(curated_terms))
if(user_terms != None):
    all_terms_list = all_terms_list.union(set(user_terms.split('\n')))
all_terms_list = list(all_terms_list)

# check user search scope
tiab_txt = ''
if(search_scope == 'TIAB'):
    search_scope = 'Title and abstract'
    tiab_txt = 'checked="checked"'
else:
    search_scope = 'Full-text'

############################### Print Result ###################################

# print query summary
print """
    <h2 align="center">
    Literature search query summary
    </h2>
"""

print '<table id="summary_tbl" align="center">'
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
print '<div align="center">'
print '<button id="confirm_btn">search in literature</button>'
print '</div>'
print '<hr/>'

# print gene implication summary
print """
    <div align="center">
    <b>Gene implication summary</b>
    <button class="show_hide" type="button">hide</button>
    <br/>
    <table id="imp_summary_tbl" class="tablesorter">
    <thead>
    <tr>
    <th>Human gene Entrez ID</th>
    <th>Human gene symbol</th>
    <th>Citeline count</th>
"""

if('tbl_ewas_gene_exp' in ewas_tables):
    print '<th>Number of EWAS loci<br/>associated with gene expression</th>'
if('tbl_ewas_prot_exp' in ewas_tables):
    print '<th>Number of EWAS loci<br/>associated with protein expression</th>'
if('tbl_ewas_trait' in ewas_tables):
    print '<th>Number of EWAS loci<br/>associated with phenotypes</th>'

if('tbl_gwas_gene_exp' in gwas_tables):
    print '<th>Number of GWAS loci<br/>associated with gene expression</th>'
if('tbl_gwas_prot_exp' in gwas_tables):
    print '<th>Number of GWAS loci<br/>associated with protein expression</th>'
if('tbl_gwas_trait' in gwas_tables):
    print '<th>Number of GWAS loci<br/>associated with phenotypes</th>'
    
if(len(ewas_tables) > 0):
    print '<th>Total number of EWAS loci </th>'
if(len(gwas_tables) > 0):
    print '<th>Total number of GWAS loci </th>'

print """
    <th>Total number of loci </th>
    </tr>
    </thead>
    <tbody>
"""

for key in combined_entrez_id_set:
    gene_sym = ''
    citeline_cnt = 0
    if(key in id_sym_trial):
        gene_sym = id_sym_trial[key]['gene_sym']
        citeline_cnt = id_sym_trial[key]['citeline_cnt']
    
    ewas_assoc_pos = dict()
    if(key in ewas_gene_support_info):
        ewas_assoc_pos = ewas_gene_support_info[key]
    num_ewas_gene_exp = safe_getval(ewas_assoc_pos, 'gene_exp')
    num_ewas_prot_exp = safe_getval(ewas_assoc_pos, 'prot_exp')
    num_ewas_trait = safe_getval(ewas_assoc_pos, 'trait')
    num_ewas_tot = num_ewas_gene_exp + num_ewas_prot_exp + num_ewas_trait
    
    gwas_assoc_pos = dict()
    if(key in gwas_gene_support_info):
        gwas_assoc_pos = gwas_gene_support_info[key]
    num_gwas_gene_exp = safe_getval(gwas_assoc_pos, 'gene_exp')
    num_gwas_prot_exp = safe_getval(gwas_assoc_pos, 'prot_exp')
    num_gwas_trait = safe_getval(gwas_assoc_pos, 'trait')
    num_gwas_tot = num_gwas_gene_exp + num_gwas_prot_exp + num_gwas_trait
    
    num_tot = num_ewas_tot + num_gwas_tot
    
    print '<tr>'
    print '<td>%s</td>' % key
    print '<td>%s</td>' % gene_sym
    print '<td>%d</td>' % citeline_cnt
    
    if('tbl_ewas_gene_exp' in ewas_tables):
        print '<td>%d</td>' % num_ewas_gene_exp
    if('tbl_ewas_prot_exp' in ewas_tables):
        print '<td>%d</td>' % num_ewas_prot_exp
    if('tbl_ewas_trait' in ewas_tables):
        print '<td>%d</td>' % num_ewas_trait
    
    if('tbl_gwas_gene_exp' in gwas_tables):
        print '<td>%d</td>' % num_gwas_gene_exp
    if('tbl_gwas_prot_exp' in gwas_tables):
        print '<td>%d</td>' % num_gwas_prot_exp
    if('tbl_gwas_trait' in gwas_tables):
        print '<td>%d</td>' % num_gwas_trait
    
    if(len(ewas_tables) > 0):
        print '<td>%d</td>' % num_ewas_tot
    if(len(gwas_tables) > 0):
        print '<td>%d</td>' % num_gwas_tot
    
    print '<td>%d</td>' % num_tot
    
    print '</tr>'
    
print """
    </tbody>
    </table>
    </div>
    <hr/>
"""

print_ewas_query_result(ewas_query_result, ewas_tables)
print_gwas_query_result(gwas_query_result, gwas_tables)

################################## Hidden Form #################################

genes_form_str = ''
for i in xrange(len(gene_list)):
    genes_form_str += str(gene_list[i][0]).strip()
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
            
            <input type="checkbox" name="tiab_only" %s/> 
        </form>
""" % (genes_form_str, terms_form_str, tiab_txt)


################################### Ending #####################################

print """
    </body>
    </html>
"""
con.close()
