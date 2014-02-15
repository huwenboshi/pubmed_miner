from consts import *
import urllib
import urllib2
import xml.etree.ElementTree as ET
import time
import sys

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
    return '<br/><br/>'.join(str_list)

# count number of occurence of terms in text
def get_term_count(term_list, text):
    text = text.lower()
    term_count = dict()
    for term in term_list:
        term_count[term] = text.count(term)
    return term_count

# initialize gene term count
def init_gene_term_cnt(gene_ids_list, terms_list):
    gene_term_count = dict()
    for gene_id in gene_ids_list:
        gene_term_count[gene_id] = dict()
        for term in terms_list:
            gene_term_count[gene_id][term] = 0
    return gene_term_count
    
################################## GENE STUFF ##################################

# get gene info, return a dictionary
def get_gene_info(gene_ids_list):

    # download gene summary from eutils esummary
    gene_ids_esummary_qstr = 'id='+','.join(gene_ids_list)
    qstr_url = esummary_url+'&db=gene&'+gene_ids_esummary_qstr
    genes_info = urllib2.urlopen(qstr_url)
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
    genes_pmids = urllib2.urlopen(qstr_url)
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
    abstract_elink = urllib2.urlopen(qstr_url)
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
        efetch = urllib2.urlopen(efetch_qstr)
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
            pmid = abstract_efetch_root[j].find('.//PMID').text
            article_url = pubmed_url+'/'+pmid
            title = abstract_efetch_root[j].find('.//ArticleTitle').text
            text_elements = abstract_efetch_root[j].findall('.//AbstractText')
            text = get_abstract_text(text_elements)
            title_and_text = title + text
            if(title_and_text.lower().find(term) >= 0):
                article_obj = {'url':article_url, 'title':title, 'text':text}
                term_abstract[term].append(article_obj)
    return term_abstract

# print abstract table for tiab search
def print_abstract_by_term_tiab(term_abstract, terms_list, gene_id):
    print '<table id="articles_by_term_gene_id_%s"' % gene_id
    print 'class="info_table">'
    print """
        <tr><td>Term</td><td colspan="2">Title & Abstract 
        Containing the Term</td></tr>
    """
    for term in terms_list:
        in_total = len(term_abstract[term])
        if(in_total > 0):
            print '<tr>'
            print '<td valign="top" rowspan="%d">%s<br/>' % (in_total, term)
            print '(%d in total)</td>' % in_total
        else:
            print '<tr><td>%s<br/>(%d in total)</td>' % (term, in_total)
            print '<td>0</td><td>N/A</td></tr>'
        for j in xrange(in_total):
            if(j != 0):
                print '<tr>'
            print '<td>%d</td>' % (j+1)
            print '<td>'
            url = term_abstract[term][j]['url']
            title = term_abstract[term][j]['title']
            print '<a href="%s" target="_blank">%s</a>' % (url, title)
            print '<div class="abstract_txt"><button class="show_more_less"'
            text = term_abstract[term][j]['text']
            print 'type="button">show more</button>%s</div>' % text
            print '</td>'
            print '</tr>'
    print '</table>'

# print abstract with sort
def print_abstract_with_sort(abstract_efetch_root, terms_list, gene_id):
    num_articles = len(abstract_efetch_root)
    num_terms = len(terms_list)
    
    # print table header
    print '<table id="articles_gene_id_%s" class="info_table">' % xstr(gene_id)
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
gene_webenv_querykey, gene_term_count):
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
    print_abstract_by_term_tiab(term_abstract, terms_list, gene_id)
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
        <hr/>
    """
    return gene_term_count

############################# FULL-TEXT STUFF ##################################

# get webenv and query key for full text search
def get_webenv_querykey_full_text(gene_id, term):
    terms_qstr = 'term=' + term
    gene_ids_qstr = 'id=' + gene_id
    qstr_url = elink_url+dbfrom_db_lknm_cmd+'&'+gene_ids_qstr+'&'+terms_qstr
    abstract_elink = urllib2.urlopen(qstr_url)
    abstract_elink_str = abstract_elink.read()
    abstract_elink_root = ET.fromstring(abstract_elink_str)
    query_key_element = abstract_elink_root[0].find('.//QueryKey')
    query_key = None
    if(query_key_element != None):
        query_key = abstract_elink_root[0].find('.//QueryKey').text
    webenv = abstract_elink_root[0].find('./WebEnv').text
    return (webenv, query_key)

# print full text search result
def print_fulltext_search_result(terms_list, gene_id, gene_term_count):

    print """
        <div>
            <b>Abstracts (grouped by terms)</b>
            <button class="show_hide" type="button">hide</button>
            <br/>
    """

    # print table header
    print '<table id="articles_by_term_gene_id_%s"' % gene_id
    print 'class="info_table">'
    print """
        <tr><td>Term</td><td colspan="2">Title & Abstract 
        Containing the Term</td></tr>
    """
    
    # iterate through terms
    for term in terms_list:
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
            efetch = urllib2.urlopen(efetch_qstr)
            efetch_str = efetch.read()
            efetch_root = ET.fromstring(efetch_str)
        
        # display result in table
        in_total = len(efetch_root)
        gene_term_count[gene_id][term] += in_total
        
        if(in_total > 0):
            print '<tr>'
            print '<td valign="top" rowspan="%d">%s<br/>' % (in_total, term)
            print '(%d in total)</td>' % in_total
        else:
            print '<tr><td>%s<br/>(%d in total)</td>' % (term, in_total)
            print '<td>0</td><td>N/A</td></tr>'
        
        # print table body
        for j in xrange(in_total):
            pmid = efetch_root[j].find('.//PMID').text
            url = pubmed_url + '/' + pmid
            title = efetch_root[j].find('.//ArticleTitle').text
            text_elements = efetch_root[j].findall('.//AbstractText')
            text = get_abstract_text(text_elements)
            if(j != 0):
                    print '<tr>'
            print '<td>%d</td>' % (j+1)
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
        <hr/>
    """
    
    return gene_term_count

########################### GWAS CATALOG STUFF #################################

# print info list
def print_nhgri_gwas_info_list(info_list):
    print """
        <div>
            <b>Gene GWAS Catalog Info</b>
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
