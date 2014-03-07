from consts import *
from utils import *

import sys
import fractions
import math

############################## CONSTANTS #######################################

# ewas user selection to table name
ewas_tbl_map = {'tbl_ewas_gene_exp': 'liver_expression_ewas',
                'tbl_ewas_prot_exp': 'liver_proteomics_ewas',
                'tbl_ewas_trait':    'clinical_metabolite_traits_ewas'}

ewas_simp_map = {'tbl_ewas_gene_exp': 'gene_exp',
                'tbl_ewas_prot_exp': 'prot_exp',
                'tbl_ewas_trait':    'trait'}

gwas_tbl_map = {'tbl_gwas_gene_exp': 'liver_expression_eqtl_gwas', 
                'tbl_gwas_prot_exp': 'liver_protein_trans_eqtl_gwas',
                'tbl_gwas_trait':    'clinical_metabolite_traits_gwas'}

gwas_simp_map = {'tbl_gwas_gene_exp': 'gene_exp',
                 'tbl_gwas_prot_exp': 'prot_exp',
                 'tbl_gwas_trait':    'trait'}

################################ INPUT #########################################

# return true if input is a number
def isnumber(s):
   try:
     float(s)
     return True
   except ValueError:
     try: 
       fractions.Fraction(s)
       return True
     except ValueError: 
       return False

# check input validity
def is_input_valid(pval, max_distance, cnt):
    if(pval == None or max_distance == None or cnt == None):
        return False
    return (isnumber(pval) and isnumber(max_distance) and isnumber(cnt))

# check if user input is valid
def check_user_input_ewas_gwas(tables, gene_exp_pval, prot_exp_pval,
    trait_pval, gene_exp_max_distance, prot_exp_max_distance,
    trait_max_distance, gene_exp_cnt, prot_exp_cnt, trait_cnt, trait_names,
    assoc_logic_sel):
    
    # check if any table is selected
    if(len(tables) == 0):
        print '<b>No association type selected</b>'
        sys.exit()
        
    # for each selected table check validity of input
    status = True
    for tbl in tables:
        if(tbl.find('gene_exp') >= 0):
            status = is_input_valid(gene_exp_pval,
                gene_exp_max_distance, gene_exp_cnt)
            if(status == False):
                break
        if(tbl.find('prot_exp') >= 0):
            status = is_input_valid(prot_exp_pval,
                prot_exp_max_distance, prot_exp_cnt)
            if(status == False):
                break
        if(tbl.find('trait') >= 0):
            status = is_input_valid(trait_pval,
                trait_max_distance, trait_cnt)
            # trait table selected, trait_names should not be empty
            if(len(trait_names) == 0):
                status = False
            if(status == False):
                break
    
    # exist if invalid input detected
    if(status == False):
        print '<b>Invalid input detected</b>'
        sys.exit()
    
    # if trait table not selected, trait_names should be empty
    if(len(trait_names) > 0
       and ('tbl_ewas_trait' not in tables)
       and ('tbl_gwas_trait' not in tables)):
       print """<b>Please check "Association with clinical and
           metabolite traits" if you want to search trait related genes</b>"""
       sys.exit()
    
    # assoc_logic_sel can take only INTERSECTION or UNION
    if(assoc_logic_sel.upper() != 'INTERSECTION' and
       assoc_logic_sel.upper() != 'UNION'):
        print "<b>Only intersection or union is supported</b>"
        sys.exit()

# check user input for ewas and gwas
def check_user_input(imp_types, imp_type_logic_sel, ewas_tables,
    ewas_gene_exp_pval, ewas_prot_exp_pval, ewas_trait_pval,
    ewas_gene_exp_max_distance, ewas_prot_exp_max_distance,
    ewas_trait_max_distance, ewas_gene_exp_cnt, ewas_prot_exp_cnt,
    ewas_trait_cnt, ewas_trait_names, ewas_assoc_logic_sel,
    gwas_tables, gwas_gene_exp_pval, gwas_prot_exp_pval, gwas_trait_pval,
    gwas_gene_exp_max_distance, gwas_prot_exp_max_distance,
    gwas_trait_max_distance, gwas_gene_exp_cnt, gwas_prot_exp_cnt,
    gwas_trait_cnt, gwas_trait_names, gwas_assoc_logic_sel):
    
    # assoc_logic_sel can take only INTERSECTION or UNION
    if(imp_type_logic_sel.upper() != 'INTERSECTION' and
       imp_type_logic_sel.upper() != 'UNION'):
        print "<b>Only intersection or union is supported</b>"
        sys.exit() 
    
    # if ewas_tables not empty ewas_imp must be in imp_types
    if(len(ewas_tables) > 0 and 'ewas_imp' not in imp_types):
        print """<b>Please check Implication by EWAS if you want to search
            by EWAS"""
        sys.exit()

    # check if there is any table for ewas implication
    if('ewas_imp' in imp_types):
        
        # if ewas implication is selected, there must be some ewas table
        cnt = 0
        for key in ewas_tbl_map:
            if(key in ewas_tables):
                cnt += 1
        if(cnt == 0):
            print """<b>At least one assocation type
                must be selected for EWAS</b>"""
            sys.exit()
            
        check_user_input_ewas_gwas(ewas_tables, ewas_gene_exp_pval, 
            ewas_prot_exp_pval, ewas_trait_pval, ewas_gene_exp_max_distance,
            ewas_prot_exp_max_distance, ewas_trait_max_distance,
            ewas_gene_exp_cnt, 
            ewas_prot_exp_cnt, ewas_trait_cnt, ewas_trait_names,
            ewas_assoc_logic_sel)
    
    # if ewas_tables not empty ewas_imp must be in imp_types
    if(len(gwas_tables) > 0 and 'gwas_imp' not in imp_types):
        print """<b>Please check Implication by GWAS if you want to search
            by GWAS"""
        sys.exit()

    # check if there is any table for gwas implication
    if('gwas_imp' in imp_types):
    
        # if gwas implication is selected, there must be some ewas table
        cnt = 0
        for key in gwas_tbl_map:
            if(key in gwas_tables):
                cnt += 1
        if(cnt == 0):
            print """<b>At least one assocation type
                must be selected for GWAS</b>"""
            sys.exit()
        
        check_user_input_ewas_gwas(gwas_tables, gwas_gene_exp_pval, 
            gwas_prot_exp_pval, gwas_trait_pval, gwas_gene_exp_max_distance,
            gwas_prot_exp_max_distance, gwas_trait_max_distance,
            gwas_gene_exp_cnt, 
            gwas_prot_exp_cnt, gwas_trait_cnt, gwas_trait_names,
            gwas_assoc_logic_sel)

############################# DATABASE #########################################

# create temporary tables for handling ewas database query
def handle_ewas_gwas_query_general(dbcon, ewas_or_gwas, tables, tbl_map,
    pval_map, dist_map, cnt_map, trait_names, assoc_logic_sel):

    # get cursor, mapping between user input and data
    cur = dbcon.cursor()
    
    # create temp tables for ewas tables
    for tbl in tables:
    
        # convert from user input to database table name
        db_table = tbl_map[tbl]
        
        # add additional constraint for searching in trait table
        trait_additional = ''
        if(db_table.find('clinical_metabolite') >= 0):
            trait_names_tmp = [];
            for i in xrange(len(trait_names)):
                clean_trait_name = trait_names[i].replace('\'','\'\'')
                trait_names_tmp.append('\''+clean_trait_name+'\'')
            trait_additional += ' and phenotype in '
            trait_additional += ' (%s) ' % (','.join(trait_names_tmp))
        
        # convert numers to strings
        pval = pval_map[tbl]
        max_distance = str(dist_map[tbl])
        citeline_cnt = str(cnt_map[tbl])
        pval_str = str(math.pow(10.0, -1.0*float(pval)))
        
        # construct query
        query = """
                create temporary table %s_tmp as
                select A.*, symbol, citeline_count from
                (select * from %s join mouse_sym_human_entrez on 
                    gene_annot_gene_sym = mouse_gene_sym
                    where (pval < %s) and 
                    (dist_gene_start_site <> 'NULL' and 
                     dist_gene_end_site <> 'NULL') and 
                    ((abs(dist_gene_start_site) < %s or
                      abs(dist_gene_end_site) < %s) or
                     (dist_gene_start_site < 0 and
                      dist_gene_end_site > 0))
                    %s) as A join clinical_trial on
                    human_entrez_id = gene_id
                    where citeline_count > %s
        """ % (db_table, db_table, pval_str, max_distance,
               max_distance, trait_additional, citeline_cnt)

        # execute query
        cur.execute(query)

    # create temp human gene id table, ewas intersection case
    if(assoc_logic_sel == 'INTERSECTION'):
    
        # create query to get intersection of human entrez ids
        sub_query = ''
        for i in xrange(len(tables)):
            tbl = tables[i]
            db_table = tbl_map[tbl]
            if(i == 0):
                sub_query = """
                    select distinct(%s_tmp.human_entrez_id) from %s_tmp
                """ % (db_table, db_table)
            if(i > 0):
                sub_query += """
                    join
                    (select distinct(human_entrez_id) from %s_tmp)
                    as %s_tmp_id
                    on
                    %s_tmp.human_entrez_id = %s_tmp_id.human_entrez_id
                """ % (db_table,db_table,tbl_map[tables[0]],db_table)
        query = """
            create temporary table human_entrez_id_%s_tmp as %s
        """ % (ewas_or_gwas, sub_query)
        
        # execute query
        cur.execute(query)

    # create temp human gene id table, ewas union case
    elif(assoc_logic_sel == 'UNION'):
    
        # create query to get union of human entrez ids
        sub_query = ''
        sub_query_list = []
        for tbl in tables:
            db_table = tbl_map[tbl]
            sub_query = """
                (select distinct(human_entrez_id) from %s_tmp 
                as %s_tmp_id)
            """ % (db_table, db_table)
            sub_query_list.append(sub_query)
        query = """
            create temporary table human_entrez_id_%s_tmp as
            select distinct(human_entrez_id) from ((%s) as A)
        """ % (ewas_or_gwas, ' union '.join(sub_query_list))
        
        # execute query
        cur.execute(query)
    
    return

# create temporary tables for handling ewas database query
def handle_ewas_query(dbcon,
                      ewas_tables,
                      ewas_gene_exp_pval, 
                      ewas_prot_exp_pval,
                      ewas_trait_pval,
                      ewas_gene_exp_max_distance,
                      ewas_prot_exp_max_distance,
                      ewas_trait_max_distance,
                      ewas_gene_exp_cnt, 
                      ewas_prot_exp_cnt,
                      ewas_trait_cnt,
                      ewas_trait_names,
                      ewas_assoc_logic_sel):

    ewas_pval_map = {'tbl_ewas_gene_exp': ewas_gene_exp_pval,
                     'tbl_ewas_prot_exp': ewas_prot_exp_pval,
                     'tbl_ewas_trait':    ewas_trait_pval}
    ewas_dist_map = {'tbl_ewas_gene_exp': ewas_gene_exp_max_distance,
                     'tbl_ewas_prot_exp': ewas_prot_exp_max_distance,
                     'tbl_ewas_trait':    ewas_trait_max_distance}
    ewas_cnt_map =  {'tbl_ewas_gene_exp': ewas_gene_exp_cnt,
                     'tbl_ewas_prot_exp': ewas_prot_exp_cnt,
                     'tbl_ewas_trait':    ewas_trait_cnt}
    
    handle_ewas_gwas_query_general(dbcon,
                                   'ewas',
                                   ewas_tables,
                                   ewas_tbl_map,
                                   ewas_pval_map,
                                   ewas_dist_map,
                                   ewas_cnt_map,
                                   ewas_trait_names,
                                   ewas_assoc_logic_sel)
    
    return

# create temporary tables for handling gwas database query
def handle_gwas_query(dbcon,
                      gwas_tables,
                      gwas_gene_exp_pval, 
                      gwas_prot_exp_pval,
                      gwas_trait_pval,
                      gwas_gene_exp_max_distance,
                      gwas_prot_exp_max_distance,
                      gwas_trait_max_distance,
                      gwas_gene_exp_cnt, 
                      gwas_prot_exp_cnt,
                      gwas_trait_cnt,
                      gwas_trait_names,
                      gwas_assoc_logic_sel):
                      
    gwas_pval_map = {'tbl_gwas_gene_exp': gwas_gene_exp_pval,
                     'tbl_gwas_prot_exp': gwas_prot_exp_pval,
                     'tbl_gwas_trait':    gwas_trait_pval}
    gwas_dist_map = {'tbl_gwas_gene_exp': gwas_gene_exp_max_distance,
                     'tbl_gwas_prot_exp': gwas_prot_exp_max_distance,
                     'tbl_gwas_trait':    gwas_trait_max_distance}
    gwas_cnt_map =  {'tbl_gwas_gene_exp': gwas_gene_exp_cnt,
                     'tbl_gwas_prot_exp': gwas_prot_exp_cnt,
                     'tbl_gwas_trait':    gwas_trait_cnt}
    
    handle_ewas_gwas_query_general(dbcon,
                                   'gwas',
                                   gwas_tables,
                                   gwas_tbl_map,
                                   gwas_pval_map,
                                   gwas_dist_map,
                                   gwas_cnt_map,
                                   gwas_trait_names,
                                   gwas_assoc_logic_sel)

    return

# handle ewas query and gwas query
def handle_query(dbcon,
                 imp_types,
                 imp_type_logic_sel,
                 ewas_tables,
                 ewas_gene_exp_pval, 
                 ewas_prot_exp_pval,
                 ewas_trait_pval,
                 ewas_gene_exp_max_distance,
                 ewas_prot_exp_max_distance,
                 ewas_trait_max_distance,
                 ewas_gene_exp_cnt, 
                 ewas_prot_exp_cnt,
                 ewas_trait_cnt,
                 ewas_trait_names,
                 ewas_assoc_logic_sel,
                 gwas_tables,
                 gwas_gene_exp_pval, 
                 gwas_prot_exp_pval,
                 gwas_trait_pval,
                 gwas_gene_exp_max_distance,
                 gwas_prot_exp_max_distance,
                 gwas_trait_max_distance,
                 gwas_gene_exp_cnt, 
                 gwas_prot_exp_cnt,
                 gwas_trait_cnt,
                 gwas_trait_names,
                 gwas_assoc_logic_sel):
    
    # initialization
    cur = dbcon.cursor()
    
    human_entrez_gene_id_set = []
    
    ewas_result = {'gene_exp': [],
                   'prot_exp': [],
                   'trait':    []}
                   
    gwas_result = {'gene_exp': [],
                   'prot_exp': [],
                   'trait':    []}

    # human entrez id tables
    human_entrez_id_tables = []

    # handle ewas query
    if('ewas_imp' in imp_types):
        human_entrez_id_tables.append('human_entrez_id_ewas_tmp')
        handle_ewas_query(dbcon,
                          ewas_tables,
                          ewas_gene_exp_pval, 
                          ewas_prot_exp_pval,
                          ewas_trait_pval,
                          ewas_gene_exp_max_distance,
                          ewas_prot_exp_max_distance,
                          ewas_trait_max_distance,
                          ewas_gene_exp_cnt, 
                          ewas_prot_exp_cnt,
                          ewas_trait_cnt,
                          ewas_trait_names,
                          ewas_assoc_logic_sel)
    
    # handle gwas query
    if('gwas_imp' in imp_types):
        human_entrez_id_tables.append('human_entrez_id_gwas_tmp')
        handle_gwas_query(dbcon,
                          gwas_tables,
                          gwas_gene_exp_pval, 
                          gwas_prot_exp_pval,
                          gwas_trait_pval,
                          gwas_gene_exp_max_distance,
                          gwas_prot_exp_max_distance,
                          gwas_trait_max_distance,
                          gwas_gene_exp_cnt, 
                          gwas_prot_exp_cnt,
                          gwas_trait_cnt,
                          gwas_trait_names,
                          gwas_assoc_logic_sel)
    # apply intersectoin
    if(imp_type_logic_sel == 'INTERSECTION'):
    
        # create query to get intersection of human entrez ids
        sub_query = ''
        for i in xrange(len(human_entrez_id_tables)):
            tbl = human_entrez_id_tables[i]
            if(i == 0):
                sub_query = """
                    select distinct(%s.human_entrez_id) from %s
                """ % (tbl, tbl)
            if(i > 0):
                sub_query += """
                    join (select distinct(human_entrez_id) from %s) as %s_sub
                    on
                    %s.human_entrez_id = %s_sub.human_entrez_id
                """ % (tbl,tbl,human_entrez_id_tables[0],tbl)
        query = """
            create temporary table human_entrez_id_ewas_gwas_tmp as %s
        """ % (sub_query)
        
        cur.execute(query)
    
    # apply union
    elif(imp_type_logic_sel == 'UNION'):

        # create query to get union of human entrez ids
        sub_query = ''
        sub_query_list = []
        for tbl in human_entrez_id_tables:
            sub_query = """
                (select distinct(human_entrez_id) from %s as %s_sub)
            """ % (tbl, tbl)
            sub_query_list.append(sub_query)
        query = """
            create temporary table human_entrez_id_ewas_gwas_tmp as
            select distinct(human_entrez_id) from ((%s) as A)
        """ % (' union '.join(sub_query_list))
        
        # execute query
        cur.execute(query)
    
    # parse out result
    query = """select * from human_entrez_id_ewas_gwas_tmp"""
    cur.execute(query)
    human_entrez_gene_id_set = fetch_from_db(cur)
    
    for table in ewas_tables:
        db_tbl = ewas_tbl_map[table]
        tmp_query = """select * from %s_tmp where human_entrez_id in
            (select * from human_entrez_id_ewas_gwas_tmp)""" % db_tbl
        query = """create temporary table %s_tmp_final 
            as %s""" % (db_tbl, tmp_query)
        cur.execute(query)
        query = """select * from %s_tmp_final""" % db_tbl
        cur.execute(query)
        ewas_result[ewas_simp_map[table]] = fetch_from_db(cur)
    
    for table in gwas_tables:
        db_tbl = gwas_tbl_map[table]
        tmp_query = """select * from %s_tmp where human_entrez_id in
            (select * from human_entrez_id_ewas_gwas_tmp)""" % db_tbl
        query = """create temporary table %s_tmp_final 
            as %s""" % (db_tbl, tmp_query)
        cur.execute(query)
        query = """select * from %s_tmp_final""" % db_tbl
        cur.execute(query)
        gwas_result[gwas_simp_map[table]] = fetch_from_db(cur)
    
    # combine result
    result = {'human_gene_set': human_entrez_gene_id_set,
              'ewas_query_result': ewas_result,
              'gwas_query_result': gwas_result}
    return result
    
############################## META ############################################

# count the number of implicating sites for each gene
def count_implicating_sites_ewas_gwas_general(dbcon, tables, tbl_map, simp_map):
    
    imp_count = {} 
    
    cur = dbcon.cursor()
    for table in tables:
        db_tbl = tbl_map[table]
        query = """
            select human_entrez_id, symbol, count(distinct site_chr, site_bp)
                from %s_tmp_final 
            group by human_entrez_id
        """ % db_tbl
        cur.execute(query)
        result = fetch_from_db(cur)
        for row in result:
            human_gene_id = row[0]
            symbol = row[1]
            count = row[2]
            assoc_type = simp_map[table]
            if(human_gene_id not in imp_count):
                imp_count[human_gene_id] = dict()
            if(assoc_type not in imp_count[human_gene_id]):
                imp_count[human_gene_id][assoc_type] = dict()
            imp_count[human_gene_id][assoc_type] = count
    
    return imp_count

# count the number of implicating sites for each gene
def count_implicating_sites_gwas(dbcon, gwas_tables):
    imp_count = count_implicating_sites_ewas_gwas_general(dbcon, gwas_tables,
        gwas_tbl_map, gwas_simp_map)
    
    return imp_count

# count the number of implicating sites for each gene
def count_implicating_sites_ewas(dbcon, ewas_tables):
    imp_count = count_implicating_sites_ewas_gwas_general(dbcon, ewas_tables,
        ewas_tbl_map, ewas_simp_map)    
    
    return imp_count

# get gene symbol citeline count data
def get_symbol_citeline_count(dbcon):
    cur = dbcon.cursor()
    id_sym_trial = dict()
    query = """select gene_id, symbol, citeline_count from clinical_trial where
               gene_id in (select * from human_entrez_id_ewas_gwas_tmp)"""
    cur.execute(query)
    result = fetch_from_db(cur)
    for row in result:
        gene_id = row[0]
        gene_sym = row[1]
        citeline_cnt = row[2]
        id_sym_trial[str(gene_id)] = dict()
        id_sym_trial[str(gene_id)]['gene_sym'] = gene_sym
        id_sym_trial[str(gene_id)]['citeline_cnt'] = citeline_cnt
    return id_sym_trial

############################# DISPLAY ##########################################

# display ewas query result
def print_ewas_query_result(ewas_query_result, ewas_tables):
    
#------------------------------ GENE TABLE ------------------------------------#

    if('tbl_ewas_gene_exp' in ewas_tables):
        print """
            <div>
                <b>Genes implicated by CG methylations associated 
                with gene expression</b>
                <button class="show_hide" type="button">hide</button>
                <br/>
        """
        gene_exp_result = ewas_query_result['gene_exp']
        print """
                <table id="ewas_gene_exp_tbl" class="tablesorter">
                <thead>
                <tr>
                    <th>mCG position</th>
                    <th>1Mb-window</th>
                    <th>Implicated mouse gene symbol</th>
                    <th>Gene position</th>
                    <th>p-value</th>
                    <th>Human homolog entrez ID</th>
                </tr>
                </thead>
                <tbody>
        """
        for result in gene_exp_result:
        
            print '<tr>'
            
            # mCG 
            print '<td>chr%s:%s</td>' % (result[5], result[6])
            
            # mCG window
            window_st = result[6]-500000 if(result[6]-500000 > 0) else 0
            window_ed = result[6]+500000
            window_str = "chr%s:%s-%s"%(result[5],str(window_st),str(window_ed))
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                window_str, window_str)
            
            # gene symbol
            print '<td>%s</td>' % result[16]
            
            # gene position
            gene_pos_str = 'chr%s:%s-%s' % (result[10], result[11], result[12])
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                gene_pos_str, gene_pos_str)
            
            # p-value
            print '<td>%s</td>' % result[0]
            
            # human entrez gene id 
            print '<td>%s</td>' % result[18]
            
            print '</tr>'
        
        print """
            </tbody>
            </table>
            </div>
        """
    
#------------------------------ PROT TABLE ------------------------------------#

    if('tbl_ewas_prot_exp' in ewas_tables):
        print """
            <div>
                <b>Genes implicated by CG methylations associated 
                with protein expression</b>
                <button class="show_hide" type="button">hide</button>
                <br/>
        """
        prot_exp_result = ewas_query_result['prot_exp']
        print """
                <table id="ewas_prot_exp_tbl" class="tablesorter">
                <thead>
                <tr>
                    <th>mCG position</th>
                    <th>1Mb-window</th>
                    <th>Implicated mouse<br/>gene symbol</th>
                    <th>Gene position</th>
                    <th>p-value</th>
                    <th>Human homolog<br/>entrez ID</th>
                </tr>
                </thead>
                <tbody>
        """
        for result in prot_exp_result:
        
            print '<tr>'
            
            # mCG 
            print '<td>chr%s:%s</td>' % (result[5],result[6])
            
            # mCG window
            window_st = result[6]-500000 if(result[6]-500000 > 0) else 0
            window_ed = result[6]+500000
            window_str = "chr%s:%s-%s"%(result[5],str(window_st),str(window_ed))
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                window_str, window_str)
            
            # gene symbol
            print '<td>%s</td>' % result[16]
            
            # gene position
            gene_pos_str = 'chr%s:%s-%s' % (result[10], result[11], result[12])
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                gene_pos_str, gene_pos_str)
            
            # p-value
            print '<td>%s</td>' % result[0]
            
            # human entrez gene id
            print '<td>%s</td>' % result[18]
            
            print '</tr>'
            
        print """
            </tbody>
            </table>
            </div>
        """
    
#----------------------------- TRAIT TABLE ------------------------------------#

    if('tbl_ewas_trait' in ewas_tables):
        print """
            <div>
                <b>Genes implicated by CG methylations associated 
                with clinical and metabolite trait</b>
                <button class="show_hide" type="button">hide</button>
                <br/>
        """
        trait_exp_result = ewas_query_result['trait']
        print """
                <table id="ewas_trait_tbl" class="tablesorter">
                <thead>
                <tr>
                    <th>mCG position</th>
                    <th>1-Mb window</th>
                    <th>Implicated mouse<br/>gene symbol</th>
                    <th>Gene position</th>
                    <th>Phenotype</th>
                    <th>Phenotype class</th>
                    <th>p-value</th>
                    <th>Human homolog<br/>entrez ID</th>
                </tr>
                </thead>
                <tbody>
        """
        for result in trait_exp_result:
            print '<tr>'
            
            # mCG position
            print '<td>chr%s:%s</td>' % (result[4],result[5])
            
            # mCG window
            window_st = result[5]-500000 if(result[5]-500000 > 0) else 0
            window_ed = result[5]+500000
            window_str = "chr%s:%s-%s" % (result[4], str(window_st),
                str(window_ed))
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                window_str, window_str)
            
            # gene symbol
            print '<td>%s</td>' % result[18]
            
            # gene position
            gene_pos_str = 'chr%s:%s-%s' % (result[12], result[13], result[14])
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                gene_pos_str, gene_pos_str)
            
            # phenotype
            print '<td>%s</td>' % result[2]
            
            # phenotype class
            print '<td>%s</td>' % result[3]
            
            # pvalue
            print '<td>%s</td>' % result[0]
            
            # human entrez id
            print '<td>%s</td>' % result[20]
            
            print '</tr>'
        
        print """
            </tbody>
            </table>
            </div>
        """

#------------------------------------------------------------------------------#

# display gwas query result
def print_gwas_query_result(gwas_query_result, gwas_tables):
    
#----------------------------- GENE TABLE -------------------------------------#

    if('tbl_gwas_gene_exp' in gwas_tables):
        print """
            <div>
                <b>Genes implicated by SNPs associated with gene expression</b>
                <button class="show_hide" type="button">hide</button>
                <br/>
        """
        gene_exp_result = gwas_query_result['gene_exp']
        print '<table id="gwas_gene_exp_tbl" class="tablesorter">'
        print """
                <thead>
                <tr>
                    <th>SNP name</th>
                    <th>SNP position</th>
                    <th>1Mb-window</th>
                    <th>Implicated mouse gene symbol</th>
                    <th>Gene position</th>
                    <th>p-value</th>
                    <th>Human homolog entrez ID</th>
                </tr>
                </thead>
                <tbody>
        """
        for result in gene_exp_result:
            print '<tr>'
            
            # snp name
            print '<td>%s</td>' % result[11]
            
            # snp pos
            print '<td>chr%s:%s</td>' % (result[8], result[9])
            
            # snp window
            snp_pos = int(result[9])
            window_st = snp_pos-500000 if(snp_pos-500000 > 0) else 0
            window_ed = snp_pos+500000
            window_str = "chr%s:%s-%s"%(result[8],str(window_st),str(window_ed))
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                window_str, window_str)
            
            # gene symbol
            print '<td>%s</td>' % result[6]
            
            # gene position
            gene_pos_str = 'chr%s:%s-%s' % (result[1], result[2], result[3])
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                gene_pos_str, gene_pos_str)
            
            # p-value
            print '<td>%s</td>' % result[7]
            
            # human entrez gene id
            print '<td>%s</td>' % result[18]
            
            print '</tr>'
        
        print """
            </tbody>
            </table>
            </div>
        """
    
#----------------------------- PROT TABLE -------------------------------------#

    if('tbl_gwas_prot_exp' in gwas_tables):
        print """
            <div>
              <b>Genes implicated by SNPs associated with protein expression</b>
              <button class="show_hide" type="button">hide</button>
              <br/>
        """
        prot_exp_result = gwas_query_result['prot_exp']
        print '<table id="gwas_prot_exp_tbl" class="tablesorter">'
        print """
                <thead>
                <tr>
                    <th>SNP name</th>
                    <th>SNP position</th>
                    <th>1Mb-window</th>
                    <th>Implicated mouse<br/>gene symbol</th>
                    <th>Gene position</th>
                    <th>p-value</th>
                    <th>Human homolog<br/>entrez ID</th>
                </tr>
                </thead>
                <tbody>
        """
        for result in prot_exp_result:
            print '<tr>'
            
            # snp name
            print '<td>%s</td>' % result[15]
            
            # snp position
            print '<td>chr%s:%s</td>' % (result[12], result[13])
            
            # snp window
            snp_pos = int(result[13])
            window_st = snp_pos-500000 if(snp_pos-500000 > 0) else 0
            window_ed = snp_pos+500000
            window_str = "chr%s:%s-%s" % (result[12], str(window_st),
                str(window_ed))
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                window_str, window_str)
            
            # gene symbol
            print '<td>%s</td>' % result[7]
            
            # gene position
            gene_pos_str = 'chr%s:%s-%s' % (result[2], result[3], result[4])
            print '<td><a target="_blank" href="%s%s">%s</a></td>' % (ucsc_url,
                gene_pos_str, gene_pos_str)
            
            # p-value
            print '<td>%s</td>' % result[11]
            
            # human entrez gene id
            print '<td>%s</td>' % result[19]
            
            print '</tr>'
        
        print """
            </tbody>
            </table>
            </div>
        """
            
#----------------------------- TRAIT TABLE -------------------------------------#

        if('tbl_gwas_trait' in gwas_tables):
            print """
                <div>
                  <b>Genes implicated by SNPs associated with phenotypes</b>
                  <button class="show_hide" type="button">hide</button>
                  <br/>
            """
            trait_result = gwas_query_result['trait']
            print '<table id="gwas_trait_tbl" class="tablesorter">'
            print """
                    <thead>
                    <tr>
                        <th>SNP name</th>
                        <th>SNP position</th>
                        <th>1Mb-window</th>
                        <th>Trait</th>
                        <th>Trait type</th>
                        <th>Implicated mouse<br/>gene symbol</th>
                        <th>Gene position</th>
                        <th>p-value</th>
                        <th>Human homolog<br/>entrez ID</th>
                    </tr>
                    </thead>
                    <tbody>
            """
            for result in trait_result:
                print '<tr>'
                
                # snp name
                print '<td>%s</td>' % result[6]
                
                # snp position
                print '<td>chr%s:%s</td>' % (result[3], result[4])
                
                # snp window
                snp_pos = int(result[4])
                window_st = snp_pos-500000 if(snp_pos-500000 > 0) else 0
                window_ed = snp_pos+500000
                window_str = "chr%s:%s-%s" % (result[3], str(window_st),
                    str(window_ed))
                print '<td><a target="_blank" href="%s%s">%s</a></td>' % (
                    ucsc_url, window_str, window_str)
                
                # trait
                print '<td>%s</td>' % result[0]
                
                # trait type
                print '<td>%s</td>' % result[1]
                
                # gene symbol
                print '<td>%s</td>' % result[20]
                
                # gene position
                gene_pos_str = 'chr%s:%s-%s' % (result[14], result[15],
                    result[16])
                print '<td><a target="_blank" href="%s%s">%s</a></td>' % (
                    ucsc_url, gene_pos_str, gene_pos_str)
                
                # p-value
                print '<td>%s</td>' % result[2]
                
                # human entrez gene id
                print '<td>%s</td>' % result[22]
                
                print '</tr>'
            
            print """
                </tbody>
                </table>
                </div>
            """
