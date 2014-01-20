#!/usr/bin/python

import sqlite3 as lite
import sys

################################# DB STUFF #####################################

# connect to database
def connect_db(db_name):
    con = None
    try:
        con = lite.connect(db_name)
        con.text_factory = str
        return con
    except lite.Error, e:
        return None
        
# search database by gene symbol
def search_nhgri_gwas_catalog(con, genesym):
    if(con == None):
        return []

    c = con.cursor()
    query = """select * from nhgri_gwas_catalog 
        where Reported_Genes like '%s' or 
        Mapped_gene like '%s'""" % (genesym, genesym)
    txt = c.execute(query)
    info_list = []
    for content in txt:
        info_dict = dict()
        info_dict['author'] = content[2]
        info_dict['date'] = content[3]
        info_dict['journal'] = content[4]
        info_dict['link'] = content[5]
        info_dict['study'] = content[6]
        info_dict['trait'] = content[7]
        info_dict['region'] = content[10]
        info_dict['reported'] = content[13]
        info_dict['mapped'] = content[14]
        info_dict['snp_allele'] = content[20]
        info_dict['pval'] = content[27]
        info_list.append(info_dict)
    return info_list
    
############################### DISPLAY STUFF ##################################

# print info list
def print_info_list(info_list):
    print """
        <div>
            <b>Gene GWAS Info</b>
            <button class="show_hide" type="button">hide</button>
            <br/>
    """

    # print table header
    print '<table class="info_table">'
    print """
        <tr>
            <td>Publish Date</td>
            <td>Article</td>
            <td>Trait</td>
            <td>Region</td>
            <td>Reported Gene</td>
            <td>Mapped Gene</td>
            <td>Strongest SNP-Risk Allele</td>
            <td>P-value</td>
        </tr>
    """
    print """
        </table>
        </div>
        <br/>
    """

# search and display in one function
def search_and_display(con, genesym):
    info_list = search_nhgri_gwas_catalog(con, genesym)
    print_info_list(info_list)
