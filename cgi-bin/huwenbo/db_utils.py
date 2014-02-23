#!/usr/bin/python

import math
import sys

from consts import *
from utils import *

############################### GENERAL QUERY ##################################

# constrcut query for searching the gene database - ewas
def make_ewas_pval_dist_query(table_name, pval, max_distance):
    query = 'select * from %s where ' % table_name
    
    # get gene expression query
    pval_str = str(math.pow(10.0, -1.0*float(pval)))
    query += ' pval < %s and ' % pval_str
    
    # get distance to gene query
    query += ' (dist_gene_start_methylation_site <> \'NULL\' and '
    query += ' dist_gene_end_methylation_site <> \'NULL\') '
    query += ' and ( '
    # distance threshold relative to gene start gene end position
    query += ' ( '
    query += ' abs(dist_gene_start_methylation_site) < %s ' % max_distance
    query += ' or '
    query += ' abs(dist_gene_end_methylation_site) < %s ' % max_distance
    query += ' ) '
    # intragenic case
    query += ' or ( '
    query += ' dist_gene_start_methylation_site < 0 '
    query += ' and dist_gene_end_methylation_site > 0 '
    query += ' ) '
    query += ' ) '
    
    return query

# constrcut query for searching the gene database - gwas
def make_gwas_pval_dist_query(table_name, pval, max_distance):
    query = 'select * from %s where ' % table_name
    
    # get gene expression query
    pval_str = str(math.pow(10.0, -1.0*float(pval)))
    query += ' pval < %s and ' % pval_str
    
    # get distance to gene query
    query += ' (dist_gene_start_snp_site <> \'NULL\' and '
    query += ' dist_gene_end_snp_site <> \'NULL\') '
    query += ' and ( '
    # distance threshold relative to gene start gene end position
    query += ' ( '
    query += ' abs(dist_gene_start_snp_site) < %s ' % max_distance
    query += ' or '
    query += ' abs(dist_gene_end_snp_site) < %s ' % max_distance
    query += ' ) '
    # intragenic case
    query += ' or ( '
    query += ' dist_gene_start_snp_site < 0 '
    query += ' and dist_gene_end_snp_site > 0 '
    query += ' ) '
    query += ' ) '
    
    return query

# construct query for mouse gene to human gene id conversion
def make_to_human_gene_id_query(from_table, from_id):
    query = 'select A.%s, B.human_entrez_id from %s' % (from_id, from_table)
    query += ' as A join mouse_mgi_human_entrez as B on '
    query += ' A.mouse_mgi_id = B.mouse_mgi_id'
    
    return query

############################ EWAS QUERY ########################################

# constrcut query for searching the gene database (clinical, metabolite trait)
def make_clinical_trait_ewas_query(table_name, pval, max_distance, trait_names):
        
    # add single quotes to trait names
    for i in xrange(len(trait_names)):
        trait_names[i] = '\''+trait_names[i]+'\''
    trait_names_sql_list = '('+','.join(trait_names)+')'
    
    query = make_ewas_pval_dist_query(table_name, pval, max_distance)
    query += ' and phenotype in %s ' % trait_names_sql_list
    
    return query

# make query for ewas
def get_ewas_query_result(gene_exp_pval, protein_exp_pval, trait_pval,
        ewas_gene_exp_max_distance, ewas_prot_exp_max_distance,
        ewas_trait_max_distance, trait_names, con, assoc_logic_sel, tables):
        
    if(con == None):
        return ([],[],[],set())
    c = con.cursor()
    
    # get gene expression result
    gene_exp_human_entrez_set = set()
    gene_exp_result_tmp = []
    # if the user wants to search in the ewas gene expression table
    if('tbl_ewas_gene_exp' in tables):
        gene_exp_result_tmp = []
        query_gene_exp = make_ewas_pval_dist_query('liver_expression_ewas',
            gene_exp_pval, ewas_gene_exp_max_distance)
        query = ' create temp table gene_exp_tmp as %s; ' % query_gene_exp
        query_result = c.execute(query)
        query = 'select * from gene_exp_tmp join mouse_probe_human_entrez on '
        query += ' gene_annot_gene_probe_id = mouse_probe_id;'
        query_result = c.execute(query)
        for result in query_result:
            gene_exp_result_tmp.append(result)
            gene_exp_human_entrez_set.add(result[len(result)-1])
    
    # get protein expression result 
    prot_exp_human_entrez_set = set()
    prot_exp_result_tmp = []
    # if the user wants to search in the protein expression table
    if('tbl_ewas_prot_exp' in tables):
        prot_exp_result_tmp = []
        query_prot_exp = make_ewas_pval_dist_query('liver_proteomics_ewas',
            protein_exp_pval, ewas_prot_exp_max_distance)
        query = ' create temp table prot_exp_tmp as %s; ' % query_prot_exp
        query_result = c.execute(query)
        query = 'select * from prot_exp_tmp join mouse_entrez_human_entrez on '
        query += ' gene_annot_gene_entrez_id = mouse_entrez_id;'
        query_result = c.execute(query)
        for result in query_result:
           prot_exp_result_tmp.append(result)
           prot_exp_human_entrez_set.add(result[len(result)-1])
    
    # get trait expression result
    trait_exp_human_entrez_set = set()
    trait_exp_result_tmp = []
    # if the user wants to search in the clinical/metabolite trait table
    if('tbl_ewas_trait' in tables):
        trait_exp_result_tmp = []
        query_trait = make_clinical_trait_ewas_query(
            'clinical_metabolite_traits_ewas', trait_pval,
            ewas_trait_max_distance, trait_names)
        query = ' create temp table trait_tmp as %s; ' % query_trait
        query_result = c.execute(query)
        query = 'select * from trait_tmp join mouse_trans_human_entrez on '
        query += ' gene_annot_trans_id = mouse_transcript_id;'
        query_result = c.execute(query)
        for result in query_result:
            trait_exp_result_tmp.append(result)
            trait_exp_human_entrez_set.add(result[len(result)-1])
    
    # map between table name and gene set
    tbl_nm_gene_set = dict()
    tbl_nm_gene_set['tbl_ewas_gene_exp'] = gene_exp_human_entrez_set
    tbl_nm_gene_set['tbl_ewas_prot_exp'] = prot_exp_human_entrez_set
    tbl_nm_gene_set['tbl_ewas_trait'] = trait_exp_human_entrez_set
    
    # intersect or union the gene ids based on the user input
    final_entrez_set = set()
    if(assoc_logic_sel == 'INTERSECTION'):
        tables_list = list(tables)
        if(len(tables_list) > 0):
            final_entrez_set = tbl_nm_gene_set[tables_list[0]]
        for i in xrange(1, len(tables_list)):
            final_entrez_set = final_entrez_set.intersection(
                tbl_nm_gene_set[tables_list[i]])
    elif((assoc_logic_sel == 'UNION')):
        for tbl in tables:
            final_entrez_set = final_entrez_set.union(
                tbl_nm_gene_set[tbl])
    
    # combine the result together
    # result[len(result)-1] gives the human gene entrez id
    gene_exp_result = []
    for result in gene_exp_result_tmp:
        if(result[len(result)-1] in final_entrez_set):
            gene_exp_result.append(result)
    prot_exp_result = []
    for result in prot_exp_result_tmp:
        if(result[len(result)-1] in final_entrez_set):
            prot_exp_result.append(result)
    trait_exp_result = []
    for result in trait_exp_result_tmp:
        if(result[len(result)-1] in final_entrez_set):
            trait_exp_result.append(result)
    
    return (gene_exp_result, prot_exp_result,
        trait_exp_result, final_entrez_set)

############################### GWAS QUERY #####################################

# remove query result if the human gene in the result is not in the gene set
def filter_query_result(query_result, gene_set):
    filtered_result = []
    for result in query_result:
        if(result[len(result)-1] in gene_set):
            filtered_result.append(result)
    return filtered_result

# extract human entrez gene ids 
def extract_human_entrez_id(query_result):
    gene_id_set = set()
    for result in query_result:
        gene_id_set.add(result[len(result)-1])
    return gene_id_set

# get genes from single gene expression table
def get_gwas_gene_exp_query_result_single(cursor, tbl_nm, pval, max_distance):
    
    # get result from cis eqtl
    gwas_gene_exp_result = []
    query = make_gwas_pval_dist_query(tbl_nm, pval, max_distance)
    query = ' create temp table %s_tmp as %s; ' % (tbl_nm, query)
    query_result = cursor.execute(query)
    query = 'select * from %s_tmp join mouse_probe_human_entrez ' % tbl_nm
    query += ' on probe_id = mouse_probe_id;'
    query_result = cursor.execute(query)
    for result in query_result:
        gwas_gene_exp_result.append(result)
    
    return gwas_gene_exp_result

# get gwas genes associated with gene expression
def get_gwas_gene_exp_query_result(con, pval, max_distance):
    
    if(con == None):
        return []
    c = con.cursor()
    
    # get result from cis and trans eqtl
    gwas_gene_exp_result_cis =  get_gwas_gene_exp_query_result_single(c,
        'gwas_liver_expression_cis_eqtl', pval, max_distance)
    gwas_gene_exp_result_trans =  get_gwas_gene_exp_query_result_single(c,
        'gwas_liver_expression_trans_eqtl', pval, max_distance)
    
    return gwas_gene_exp_result_cis + gwas_gene_exp_result_trans

# get gwas genes associated with protein expression
def get_gwas_prot_exp_query_result(con, pval, max_distance):

    if(con == None):
        return []
    c = con.cursor()
    
    gwas_prot_exp_result = []
    query_tmp = make_gwas_pval_dist_query('gwas_liver_protein_trans_eqtl',
        pval, max_distance)
    query = 'create temp table gwas_liver_protein_trans_eqtl_tmp '
    query += ' as %s; ' % query_tmp
    query_result = c.execute(query)
    query = 'select * from gwas_liver_protein_trans_eqtl_tmp join '
    query += ' mouse_sym_human_entrez '
    query += ' on gene_symbol = mouse_gene_sym;'
    query_result = c.execute(query)
    
    for result in query_result:
        gwas_prot_exp_result.append(result)
    
    return gwas_prot_exp_result

# get gwas query result from all gwas tables
def get_gwas_query_result(gene_exp_pval, prot_exp_pval,
        gene_exp_max_distance, prot_exp_max_distance, con,
        assoc_logic_sel, tables):
    
    # get result from database
    gene_exp_result_tmp = get_gwas_gene_exp_query_result(con, gene_exp_pval,
        gene_exp_max_distance)
    gene_exp_human_entrez_set = extract_human_entrez_id(gene_exp_result_tmp)
    prot_exp_result_tmp = get_gwas_prot_exp_query_result(con, prot_exp_pval,
        prot_exp_max_distance)
    prot_exp_human_entrez_set = extract_human_entrez_id(prot_exp_result_tmp)
    
    # map between table name and gene set
    tbl_nm_gene_set = dict()
    tbl_nm_gene_set['tbl_gwas_gene_exp'] = gene_exp_human_entrez_set
    tbl_nm_gene_set['tbl_gwas_prot_exp'] = prot_exp_human_entrez_set
    
    # intersect or union the gene ids based on the user input
    final_entrez_set = set()
    if(assoc_logic_sel == 'INTERSECTION'):
        tables_list = list(tables)
        if(len(tables_list) > 0 and tables_list[0] != 'tbl_gwas_trait'):
            final_entrez_set = tbl_nm_gene_set[tables_list[0]]
        for i in xrange(1, len(tables_list)):
            if(tables_list[i] == 'tbl_gwas_trait'):
                continue
            final_entrez_set = final_entrez_set.intersection(
                tbl_nm_gene_set[tables_list[i]])
    elif((assoc_logic_sel == 'UNION')):
        for tbl in tables:
            if(tbl == 'tbl_gwas_trait'):
                continue
            final_entrez_set = final_entrez_set.union(
                tbl_nm_gene_set[tbl])
    
    # get final result
    gene_exp_result = filter_query_result(gene_exp_result_tmp,final_entrez_set)
    prot_exp_result = filter_query_result(prot_exp_result_tmp,final_entrez_set)
    
    return (gene_exp_result, prot_exp_result, final_entrez_set)
