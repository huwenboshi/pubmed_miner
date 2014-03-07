#!/UCSC/Panel-Auxiliary/python/Python-2.7.3/python -u

from utils import *
from consts import *
import urllib
import urllib2
import xml.etree.ElementTree as ET
import time
import sys
import os

os.umask(0o027)

import networkx as nx

os.environ['HOME'] = '/tmp/'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dir_prefix = '/UCSC/Apache-2.2.11/htdocs-UCLApanel'
dir_prefix = '..'

#################################### HELPERS ###################################

# check if string is None type
def xstr(s):
    if s is None:
        return ''
    return str(s)

# given XML text elements, parse out abstracts
def get_abstract_text(text_elements):
    str_list = []
    for element in text_elements:
        str_list.append(element.text)
    return '<br/>'.join(str_list)

# count number of occurence of terms in text
def get_term_count(term_list, text):
    text = text.lower()
    chars = """ '";:.>,<?/}]{[|\\+=_-!@#$%^&*()"""
    term_count = dict()
    for term in term_list:
        term_count_list = []
        for char in chars:
            term_count_list.append(text.count(term+char))
        term_count[term] = sum(term_count_list)
    return term_count

# initialize gene term count
def init_gene_term_cnt(gene_ids_list, terms_list):
    gene_term_count = dict()
    for gene_id in gene_ids_list:
        gene_term_count[gene_id] = dict()
        for term in terms_list:
            gene_term_count[gene_id][term] = 0
    return gene_term_count

# safe url open - catches exceptions if there is any
def safe_urlopen(url):
    data = None
    
    for i in xrange(5):
        try:
            data = urllib2.urlopen(url)
            break
        except urllib2.HTTPError as e:
            if(e.code == 414):
                print """<b>Query too long, please reduce
                    the number of genes or terms</b>"""
                sys.exit()
        except:
            continue
    
    if(data == None):
        print """Network error, please try again later"""
        sys.exit()
    
    return data
    
################################## GENE STUFF ##################################

# get gene info, return a dictionary
def get_gene_info(gene_ids_list):

    # download gene summary from eutils esummary
    gene_ids_esummary_qstr = 'id='+','.join(gene_ids_list)
    qstr_url = esummary_url+'&db=gene&'+gene_ids_esummary_qstr
    genes_info = None
    genes_info = safe_urlopen(qstr_url)
    genes_info_str = genes_info.read()
    genes_info_root = ET.fromstring(genes_info_str)
    
    # parse gene summary xml
    genes_info_list = []
    for i in xrange(len(gene_ids_list)):
        gene_info = genes_info_root[i]
        gene_id = gene_info.find('Id').text
        gene_name = gene_info.find("./Item[@Name='Name']").text
        tmp = dict()
        tmp['id'] = gene_id
        tmp['name'] = gene_name
        tmp['descp'] = gene_info.find("./Item[@Name='Description']").text
        tmp['aliases'] = gene_info.find("./Item[@Name='OtherAliases']").text
        tmp['desg'] = gene_info.find("./Item[@Name='OtherDesignations']").text
        tmp['ns'] = gene_info.find("./Item[@Name='NomenclatureSymbol']").text
        tmp['nn'] = gene_info.find("./Item[@Name='NomenclatureName']").text
        tmp['summary'] = gene_info.find("./Item[@Name='Summary']").text
        tmp['gene_url'] = gene_web_url + gene_id
        tmp['rna_url'] = genomernai_url + gene_id
        tmp['biogps_url'] = biogpsorg_url + gene_id
        tmp['genecards_url'] = genecards_url + gene_name
        genes_info_list.append(tmp)
    
    return genes_info_list

# get all gene related article count
def get_genes_pmid_count(gene_ids_list):

    # download lists of pmids for all genes
    gene_ids_elink_qstr = 'id='+'&id='.join(gene_ids_list)
    qstr_url = elink_url + dbfrom_db_lknm + gene_ids_elink_qstr
    
    genes_pmids = safe_urlopen(qstr_url)
    genes_pmids_str = genes_pmids.read()
    genes_pmids_root = ET.fromstring(genes_pmids_str)
    
    # parse the count
    genes_pmids_cnt = []
    for i in xrange(len(gene_ids_list)):
        gene_pmids = genes_pmids_root[i]
        cnt = len(gene_pmids.findall("./LinkSetDb/Link"))
        genes_pmids_cnt.append(cnt)
        
    return genes_pmids_cnt

# print gene info table
def print_gene_info(gene_info, gene_pmids_cnt):
    gene_id = gene_info['id']
    gene_name = gene_info['name']
    gene_descp = gene_info['descp']
    gene_aliases = gene_info['aliases']
    gene_desg = gene_info['desg']
    gene_ns = gene_info['ns']
    gene_nn = gene_info['nn']
    gene_summary = gene_info['summary']
    gene_url = gene_info['gene_url']
    rna_url = gene_info['rna_url']
    biogps_url = gene_info['biogps_url']
    genecards_url = gene_info['genecards_url']
    gene_pmids_cnt = gene_pmids_cnt
    print '<table class="info_table">'
    print '<tr><td>Id</td><td>%s</td></tr>'%gene_id
    print '<tr><td>Name</td><td>%s</td></tr>'%xstr(gene_name)
    print '<tr><td>Description</td><td>%s</td></tr>'%xstr(gene_descp)
    print '<tr><td>Aliases</td><td>%s</td></tr>'%xstr(gene_aliases)
    print '<tr><td>Other Designations</td><td>%s</td></tr>'%xstr(gene_desg)
    print '<tr><td>Summary</td><td>%s</td></tr>'%xstr(gene_summary)
    print '<tr><td>Gene URL</td><td>'
    print '<a href="%s" target="_blank">%s</a>'%(gene_url,gene_url)
    print '</td></tr>'
    print '<tr><td>RNA URL</td><td>'
    print '<a href="%s" target="_blank">%s</a>'%(rna_url,rna_url)
    print '</td></tr>'
    print '<tr><td>BioGPS URL</td><td>'
    print '<a href="%s" target="_blank">%s</a>'%(biogps_url,biogps_url)
    print '</td></tr>'
    print '<tr><td>GeneCards URL</td><td>'
    print '<a href="%s" target="_blank">%s</a>'%(genecards_url,genecards_url)
    print '</td></tr>'
    print '<tr><td>Total PMID Count</td><td>%d</td></tr>'%gene_pmids_cnt
    print '</table>'

############################### TIAB STUFF #####################################

# get webenv and querykey
def get_webenv_querykey(gene_ids_list, terms_list):

    # download WebEnv and QueryKey
    terms_qstr = 'term='+'+OR+'.join([s+"[TIAB]" for s in terms_list])
    terms_qstr = terms_qstr.replace(' ', '+')
    gene_ids_qstr = 'id='+'&id='.join(gene_ids_list)
    qstr_url = elink_url+dbfrom_db_lknm_cmd+'&'+gene_ids_qstr+'&'+terms_qstr
    abstract_elink = safe_urlopen(qstr_url)
    abstract_elink_str = abstract_elink.read()
    
    # parse out WebEnv and QueryKey
    abstract_elink_root = ET.fromstring(abstract_elink_str)
    webenv_key = dict()
    for i in xrange(len(abstract_elink_root)):
        gene_id = abstract_elink_root[i].find('.//Id').text
        query_key_element = abstract_elink_root[i].find('.//QueryKey')
        query_key = None
        if(query_key_element != None):
            query_key = abstract_elink_root[i].find('.//QueryKey').text
        webenv = abstract_elink_root[i].find('./WebEnv').text
        webenv_key[gene_id] = (webenv, query_key)
    
    return webenv_key

# fetch abstract for tiab search, returnx xml root
def fetch_abstract_tiab(webenv_querykey):
    
    # get article abstract from eutil server
    webenv = webenv_querykey[0]
    querykey = webenv_querykey[1]
    efetch_qstr = None
    efetch_root = []
    if(querykey != None):
        efetch_qstr = efetch_url+efetch_fixed_param
        efetch_qstr += '&WebEnv='+webenv+'&query_key='+querykey
        efetch = safe_urlopen(efetch_qstr)
        efetch_str = efetch.read()
        efetch_root = ET.fromstring(efetch_str)
        
    return efetch_root

# group abstract by term occurence
def group_abstract_by_term(abstract_efetch_root, terms_list):
    term_count = dict()
    term_abstract = dict()
    for term in terms_list:
        if(term not in term_abstract):
            term_abstract[term] = []
        for j in xrange(len(abstract_efetch_root)):
            try:
                pmid = abstract_efetch_root[j].find('.//PMID').text
                article_url = pubmed_url+'/'+pmid
                title = abstract_efetch_root[j].find('.//ArticleTitle').text
                text_elements=abstract_efetch_root[j].findall('.//AbstractText')
                text = get_abstract_text(text_elements)
                title_and_text = title + text
                if(title_and_text.lower().find(term) >= 0):
                    article_obj = {'url':article_url,'title':title,'text':text}
                    term_abstract[term].append(article_obj)
            except:
                continue
    return term_abstract

# print abstract table for tiab search
def print_abstract_by_term_tiab(term_abstract, terms_list, gene_id, id_sym):
    print '<table id="articles_by_term_gene_id_%s"' % gene_id
    print 'class="abstract_table">'
    print """
        <tr><td>Term</td><td colspan="2">Title & Abstract 
        Containing the Term</td></tr>
    """
    for term in terms_list:
        td_id = term+"_"+gene_id
        in_total = len(term_abstract[term])
        if(in_total > 0):
            print '<tr>'
            print """<td id="%s" valign="top" rowspan="%d">
                        %s<br/>
                        (%d in total)<br/>
                        <a href="#gene_id_%s">%s(%s)</a><br/>
                        <a href="#top">Return to Top</a>
                     </td>""" % (td_id, in_total, term, in_total,
                        gene_id, id_sym[gene_id], gene_id)

        for j in xrange(in_total):
            if(j != 0):
                print '<tr>'
            print '<td>%d<br/><a href="#%s">back</a></td>'%((j+1),
                term+"_"+gene_id)
            print '<td>'
            url = term_abstract[term][j]['url']
            title = term_abstract[term][j]['title']
            print '<a href="%s" target="_blank">%s</a>' % (url, title)
            """
            print '<div class="abstract_txt"><button class="show_more_less"'
            text = term_abstract[term][j]['text']
            print 'type="button">show more</button>%s</div>' % text
            print '</td>'
            """
            print '</tr>'
    print '</table>'

# print abstract with sort
def print_abstract_with_sort(abstract_efetch_root, terms_list, gene_id):
    num_articles = len(abstract_efetch_root)
    num_terms = len(terms_list)
    
    # print table header
    print '<table id="articles_gene_id_%s" class="abstract_table">' % xstr(gene_id)
    print '<tr><td></td><td>Title & Abstract'
    print ' (%d related articles in total)</td>' % num_articles
    for term in terms_list:
        print '<td><input type="checkbox" '
        print 'id="local_opt_gene_id_%s_%s" class="local_opt">'%(gene_id,term)
        print '%s</td>' % term
    print '</tr>'
    
    # print table body
    for j in xrange(num_articles):
        pmid = abstract_efetch_root[j].find('.//PMID').text
        title = abstract_efetch_root[j].find('.//ArticleTitle').text
        text_elements = abstract_efetch_root[j].findall('.//AbstractText')
        text = get_abstract_text(text_elements)
        print '<tr>'
        print '<td>%d</td>' % (j+1)
        print '<td><a href="%s/%s"' % (pubmed_url, pmid)
        print 'target="_blank">%s</a>' % title
        print '<div class="abstract_txt">'
        print '<button class="show_more_less" type="button">'
        print 'show more</button>%s</div></td>' % text
        term_count = get_term_count(terms_list, title+' '+text)
        for term in terms_list:
            print '<td><a class="term_count">%d</a>' % term_count[term]
            print '<br/>%s</td>' % term
        print '</tr>'
    
    # end article summary table
    print '</table>'

# print tiab search result
def print_tiab_search_result(terms_list, gene_id, 
    gene_webenv_querykey, gene_term_count, id_sym):
    abstract_efetch_root = fetch_abstract_tiab(gene_webenv_querykey[gene_id])
    term_abstract = group_abstract_by_term(abstract_efetch_root, terms_list)
    for term in terms_list:
        gene_term_count[gene_id][term] += len(term_abstract[term])
    print """
        <div>
            <b>Abstracts (grouped by terms)</b>
            <button class="show_hide" type="button">hide</button>
            <br/>
    """
    print_abstract_by_term_tiab(term_abstract, terms_list, gene_id, id_sym)
    #print """
    #    </div>
    #    <br/>
    #    <div>
    #    <a><b>Abstracts (choose terms to sort by the sum of their
    #    occurences in titles and abstracts)</b></a>
    #    <button class="show_hide" type="button">hide</button>
    #    <br/>
    #"""
    #print_abstract_with_sort(abstract_efetch_root, terms_list, gene_id)
    print """
        </div>
    """
    time.sleep(1)
    return gene_term_count

############################# FULL-TEXT STUFF ##################################

# get webenv and query key for full text search
def get_webenv_querykey_full_text(gene_id, term):
    terms_qstr = 'term=' + term
    gene_ids_qstr = 'id=' + gene_id
    qstr_url = elink_url+dbfrom_db_lknm_cmd+'&'+gene_ids_qstr+'&'+terms_qstr
    abstract_elink = safe_urlopen(qstr_url)
    abstract_elink_str = abstract_elink.read()
    abstract_elink_root = ET.fromstring(abstract_elink_str)
    query_key_element = abstract_elink_root[0].find('.//QueryKey')
    query_key = None
    if(query_key_element != None):
        query_key = abstract_elink_root[0].find('.//QueryKey').text
    webenv = abstract_elink_root[0].find('./WebEnv').text
    return (webenv, query_key)

# print full text search result
def print_fulltext_search_result(terms_list, gene_id, gene_term_count, id_sym):

    print """
        <div>
            <b>Abstracts (grouped by terms)</b>
            <button class="show_hide" type="button">hide</button>
            <br/>
    """

    # print table header
    print '<table id="articles_by_term_gene_id_%s"' % gene_id
    print 'class="abstract_table">'
    print """
        <tr><td>Term</td><td colspan="2">Title & Abstract 
        Containing the Term</td></tr>
    """
    
    # iterate through terms
    for i in xrange(len(terms_list)):
        term = terms_list[i]
        
        # get WebEnv and QueryKey
        webenv_querykey = get_webenv_querykey_full_text(gene_id, term)
        webenv = webenv_querykey[0]
        query_key = webenv_querykey[1]
        
        # sleep for 1 second to keep query under limit
        time.sleep(1)
        
        # fetch abstract
        efetch_root = []
        if(query_key != None):
            efetch_qstr = efetch_url+efetch_fixed_param
            efetch_qstr += '&WebEnv='+webenv+'&query_key='+query_key
            efetch = safe_urlopen(efetch_qstr)
            efetch_str = efetch.read()
            efetch_root = ET.fromstring(efetch_str)
        
        # display result in table
        in_total = len(efetch_root)
        gene_term_count[gene_id][term] += in_total
        td_id = term+"_"+gene_id
        if(in_total > 0):
            print """<td id="%s" valign="top" rowspan="%d">
                        %s<br/>
                        (%d in total)<br/>
                        <a href="#gene_id_%s">%s(%s)</a><br/>
                        <a href="#top">Return to Top</a>
                     </td>""" % (td_id, in_total, term, in_total,
                        gene_id, id_sym[gene_id], gene_id)
        
        # print table body
        for j in xrange(in_total):
            pmid = efetch_root[j].find('.//PMID').text
            url = pubmed_url + '/' + pmid
            title = efetch_root[j].find('.//ArticleTitle').text
            text_elements = efetch_root[j].findall('.//AbstractText')
            text = get_abstract_text(text_elements)
            if(j != 0):
                    print '<tr>'
            print '<td>%d<br/><a href="#%s">back</a></td>'%((j+1),
                term+"_"+gene_id)
            print '<td>'
            print '<a href="%s" target="_blank">%s</a>' % (url, title)
            print '<div class="abstract_txt"><button class="show_more_less"'
            print 'type="button">show more</button>%s</div>' % text
            print '</td>'
            print '</tr>'
            
        time.sleep(1)
        
    print """
        </table>
        </div>
    """
    
    return gene_term_count

########################### GWAS CATALOG STUFF #################################

# search database by gene symbol
def search_nhgri_gwas_catalog(con, genesym):
    if(con == None):
        return []

    cur = con.cursor()
    query = """select * from nhgri_gwas_catalog 
        where Reported_Genes like '%["""+genesym+"""]%'"""
    cur.execute(query)
    info_list = []
    for content in cur.fetchall():
        info_dict = dict()
        info_dict['author'] = content[2]
        info_dict['date'] = content[3]
        info_dict['journal'] = content[4]
        info_dict['link'] = content[5]
        info_dict['study'] = content[6]
        info_dict['trait'] = content[7]
        info_dict['region'] = content[10]
        info_dict['reported'] = content[13].replace('[','').replace(']','<br/>')
        info_dict['mapped'] = content[14]
        info_dict['snp_allele'] = content[20]
        info_dict['pval'] = content[27]
        info_list.append(info_dict)
    return info_list

# print info list
def print_nhgri_gwas_info_list(info_list):
    print """
        <div>
            <b>Gene NHGRI GWAS Catalog Info</b>
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
########################### CLINICAL TRIAL STUFF ###############################

# search database by gene id
def search_clinical_trial(con, gene_id):
    if(con == None):
        return []
    cur = con.cursor()
    query = """select * from clinical_trial 
        where gene_id = %s""" % gene_id
    cur.execute(query)
    return fetch_from_db(cur)

# print clinical trial result
def print_clinical_trial_result(result):
    print """
        <div>
            <b>Gene Clinical Trial Info</b>
            <button class="show_hide" type="button">hide</button>
            <br/>
    """

    # print table header
    if(len(result) > 0):
        row = result[0]
        print '<table class="info_table">'
        print '<tr><td>Gene ID</td><td>%d</td></tr>' % row[0]
        print '<tr><td>Gene symbol</td><td>%s</td></tr>' % row[1]
        print '<tr><td>Citeline count<td>%s</td></tr>' % row[2]
        print '<tr><td>Small molecule count<td>%s</td></tr>' % row[3]
        print '<tr><td>Small molecule launched count<td>%s</td></tr>' % row[4]
        print '<tr><td>Small molecule registered count<td>%s</td></tr>'%row[5]
        print """<tr><td>Small molecule 
                         pre-registered count<td>%s</td></tr>"""%row[6]
        print '<tr><td>Small molecule phase 3 count<td>%s</td></tr>'%row[7]
        print '<tr><td>Small molecule phase 2 count<td>%s</td></tr>'%row[8]
        print '<tr><td>Small molecule phase 1 count<td>%s</td></tr>'%row[9]
        print """<tr><td>Small molecule
                         pre-clinical count<td>%s</td></tr>"""%row[10]
        print """<tr><td>Small molecule no
                         development count<td>%s</td></tr>"""%row[11]
        print """<tr><td>Small molecule
                         discontinued count<td>%s</td></tr>"""%row[12]
        print '<tr><td>Small molecule suspended count<td>%s</td></tr>'%row[13]
        print '<tr><td>Small molecule withdrawn count<td>%s</td></tr>'%row[14]
        print '<tr><td>Antibody count<td>%s</td></tr>' % row[15]
        print '<tr><td>Antibody launched count<td>%s</td></tr>' % row[16]
        print '<tr><td>Antibody registered count<td>%s</td></tr>'%row[17]
        print '<tr><td>Antibody pre-registered count<td>%s</td></tr>'%row[18]
        print '<tr><td>Antibody phase 3 count<td>%s</td></tr>'%row[19]
        print '<tr><td>Antibody phase 2 count<td>%s</td></tr>'%row[20]
        print '<tr><td>Antibody phase 1 count<td>%s</td></tr>'%row[21]
        print '<tr><td>Antibody pre-clinical count<td>%s</td></tr>'%row[22]
        print '<tr><td>Antibody no development count<td>%s</td></tr>'%row[23]
        print '<tr><td>Antibody discontinued count<td>%s</td></tr>'%row[24]
        print '<tr><td>Antibody suspended count<td>%s</td></tr>'%row[25]
        print '<tr><td>Antibody withdrawn count<td>%s</td></tr>'%row[26]
        print """
            </table>
        """
    
    print """
            </div>
            <br/>
          """

########################### GENE CLUSTER INFO ##################################

# search database by gene id
def search_hypervariable_count(con, gene_id):
    if(con == None):
        return []
    cur = con.cursor()
    query = """select * from mouse_gene_hypervariable_count 
        where human_entrez_id = %s""" % gene_id
    cur.execute(query)
    return fetch_from_db(cur)

# print gene cluster info
def print_hypervariable_count(result):
    print """
        <div>
            <b>Number of Hyper/Hypervariable methylation
               sites around Gene on Mouse Genome</b>
            <button class="show_hide" type="button">hide</button>
            <br/>
    """
    
    print '<table class="info_table">'
    print '<tr><td>Human Entrez ID</td>'
    for row in result:
        print '<td>%d</td>' % row[0]
    print '</tr>'    
    
    print '<tr><td>Mouse homolog gene symbol</td>'
    for row in result:
        print '<td>%s</td>' % row[1]
    print '</tr>'
    
    print '<tr><td>Number of intragenic hypermethylation sites</td>'
    for row in result:
        print '<td>%d</td>' % row[2]
    print '</tr>'
    
    print '<tr><td>Number of hypermethylation sites 100 kb upstream</td>'
    for row in result:
        print '<td>%d</td>' % row[3]
    print '</tr>'
    
    print """<tr><td>Number of hypervariable methylation
                     between 100kb and 1mb upstream</td>"""
    for row in result:
        print '<td>%d</td>' % row[4]
    print '</tr>'
    
    print '<tr><td>Number of hypervariable CGs 100 kb downstream</td>'
    for row in result:
        print '<td>%d</td>' % row[5]
    print '</tr>'
    
    print '<tr><td>Number of hypervariable CGs 1mb downstream</td>'
    for row in result:
        print '<td>%d</td>' % row[6]
    print '</tr>'
        
    print '</table>'
    print """
            </div>
            <br/>
          """

############################## GENE TERM NETWORK ###############################

def draw_gene_term_graph(graph, filename, node_size=400, node_alpha=0.5,
               node_text_size=8, edge_color='blue', edge_alpha=0.3,
               edge_tickness=1, edge_text_pos=0.3,text_font='sans-serif'):

    # set figure size
    plt.figure(figsize=(30,30)) 
    
    # create networkx graph
    G=nx.Graph()

    labels = dict()
    
    # gene node
    gene_nodes = set()
    for edge in graph:
        gene_nodes.add(edge[0])
        labels[edge[0]] = edge[0]
    G.add_nodes_from(gene_nodes)
    
    # term node
    term_nodes = set()
    for edge in graph:
        term_nodes.add(edge[1])
        labels[edge[1]] = edge[1]
    G.add_nodes_from(term_nodes)

    # add edge
    G.add_edges_from(set(graph))

    # set style
    graph_pos=nx.graphviz_layout(G)

    # draw nodes
    nx.draw_networkx_nodes(G,graph_pos,nodelist=gene_nodes, 
        node_color='r', node_size=node_size, alpha=node_alpha)
    nx.draw_networkx_nodes(G,graph_pos,nodelist=term_nodes,
        node_color='b', node_size=node_size, alpha=node_alpha)
    
    # draw edges
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
        alpha=edge_alpha,edge_color=edge_color)
    
    # draw labels                       
    nx.draw_networkx_labels(G,graph_pos,labels,font_size=node_text_size,
        font_family=text_font)

    # save graph
    plt.axis('off')
    plt.savefig(filename)

    # change permission
    # os.chown(filename, 407, 507)

# display the network
def create_gene_term_network(gene_term_count, gene_ids_list,
    terms_list, id_sym):

    # get connected nodes
    graph = []
    
    for i in xrange(len(terms_list)):
        for j in xrange(len(gene_ids_list)):
            term = terms_list[i]
            gene_id = gene_ids_list[j]
            if(gene_term_count[gene_id][term] > 1):
                graph.append((id_sym[gene_id], term))
                    
    # draw the graph
    fnm = dir_prefix+'/tmp/gene_term_network_%d.png' % os.getpid()
    fnm_html = '../tmp/gene_term_network_%d.png' % os.getpid()
    draw_gene_term_graph(graph, filename=fnm)
 
    # display the graph
    print """<a target="_blank" href="%s">
            <img src="%s" height="200" width="200">
            </a>""" % (fnm_html, fnm_html)

############################## TERM TERM NETWORK ###############################

def draw_x_x_graph(graph,filename,fsr=30, fsc=30, node_size=400,
               node_alpha=0.5,ncolor='b',
               node_text_size=8, edge_color='blue', edge_alpha=0.3,
               edge_tickness=1, edge_text_pos=0.3,text_font='sans-serif'):

    # set figure size
    plt.figure(figsize=(fsr,fsc)) 
    
    # create networkx graph
    G=nx.Graph()

    labels = dict()
    
    # term node
    nodes = set()
    for edge in graph:
        nodes.add(edge[0])
        nodes.add(edge[1])
        labels[edge[0]] = edge[0]
        labels[edge[1]] = edge[1]
    G.add_nodes_from(nodes)

    # add edge
    G.add_edges_from(set(graph))

    # set style
    graph_pos=nx.graphviz_layout(G)

    # draw nodes
    nx.draw_networkx_nodes(G,graph_pos,nodelist=nodes,
        node_color=ncolor, node_size=node_size, alpha=node_alpha)
    
    # draw edges
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
        alpha=edge_alpha,edge_color=edge_color)
    
    # draw labels                       
    nx.draw_networkx_labels(G,graph_pos,labels,font_size=node_text_size,
        font_family=text_font)

    # save graph
    plt.axis('off')
    plt.savefig(filename)

    # change permission
    # os.chown(filename, 407, 507)

# create term term network
def create_term_term_network(gene_term_count, gene_ids_list,
    terms_list, id_sym):
    
    graph = []
    for i in xrange(len(terms_list)):
        for j in xrange(len(terms_list)):
            term_i = terms_list[i]
            term_j = terms_list[j]
            if(i != j):
                count_term_i = []
                count_term_j = []
                for k in xrange(len(gene_ids_list)):
                    gene_id = gene_ids_list[k]
                    count_term_i.append(gene_term_count[gene_id][term_i])
                    count_term_j.append(gene_term_count[gene_id][term_j])
                cor = pearson(count_term_i, count_term_j)
                if(cor > 0.8 and cor < 1.0):
                    graph.append((term_i, term_j))
    # draw the graph
    fnm = dir_prefix+'/tmp/term_term_network_%d.png' % os.getpid()
    fnm_html = '../tmp/term_term_network_%d.png' % os.getpid()
    draw_x_x_graph(graph, filename=fnm)
    
    # display the graph
    print """<a target="_blank" href="%s">
            <img src="%s" height="200" width="200">
            </a>""" % (fnm_html, fnm_html)

# create gene gene network
def create_gene_gene_network(gene_term_count, gene_ids_list,
    terms_list, id_sym):
    
    graph = []
    for i in xrange(len(gene_ids_list)):
        for j in xrange(len(gene_ids_list)):
            gene_i = gene_ids_list[i]
            gene_j = gene_ids_list[j]
            if(i != j):
                count_gene_i = []
                count_gene_j = []
                for k in xrange(len(terms_list)):
                    term = terms_list[k]
                    count_gene_i.append(gene_term_count[gene_i][term])
                    count_gene_j.append(gene_term_count[gene_j][term])
                cor = pearson(count_gene_i, count_gene_j)
                if(cor > 0.8 and cor < 1.0):
                    graph.append((id_sym[gene_i], id_sym[gene_j]))
    # draw the graph
    fnm = dir_prefix+'/tmp/gene_gene_network_%d.png' % os.getpid()
    fnm_html = '../tmp/gene_gene_network_%d.png' % os.getpid()
    draw_x_x_graph(graph, fsr=10, fsc=10, ncolor='r', filename=fnm)
    
    # display the graph
    print """<a target="_blank" href="%s">
            <img src="%s" height="200" width="200">
            </a>""" % (fnm_html, fnm_html)
