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
        where Reported_Genes like '%"""+genesym+"""%'"""
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
            <td class="gwas_gene">Reported Gene</td>
            <td class="gwas_title">Publication</td>
            <td class="gwas_trait">Disease/Trait</td>
            <td>Region</td>
            <td class="gwas_snp">Strongest SNP-Risk Allele</td>
            <td>P-value</td>
        </tr>
    """
    for info_dict in info_list:
        date = info_dict['date']
        link = info_dict['link']
        study = info_dict['study']
        trait = info_dict['trait']
        region = info_dict['region']
        reported = info_dict['reported']
        mapped = info_dict['mapped']
        snp_allele = info_dict['snp_allele']
        pval = info_dict['pval']
        print '<tr>'
        print '<td class="gwas_gene">%s</td>' % reported
        print '<td class="gwas_title"><a href="%s">%s</a></td>' % (link, study)
        print '<td class="gwas_trait">%s</td>' % trait
        print '<td>%s</td>' % region
        print '<td class="gwas_snp">%s</td>' % snp_allele
        print '<td>%s</td>' % pval
        print '</tr>'
    print """
        </table>
        </div>
        <br/>
    """

# search and display in one function
def search_and_display(con, genesym):
    info_list = search_nhgri_gwas_catalog(con, genesym)
    print_info_list(info_list)
