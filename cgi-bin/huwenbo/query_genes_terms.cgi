#!/usr/bin/python -u

import xml.etree.ElementTree as ET
import urllib2
import cgitb
import time
import cgi
import sys
import codecs

################################### Functions #################################
def xstr(s):
    if s is None:
        return ''
    return str(s)

def get_abstract_text(text_elements):
    str_list = []
    for element in text_elements:
        str_list.append(element.text)
    return '<br/><br/>'.join(str_list)

def get_term_count(term_list, text):
    text = text.lower()
    term_count = dict()
    for term in term_list:
        term_count[term] = text.count(term)
    return term_count

##################################### Handle ##################################
sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
cgitb.enable()

# get form data
form = cgi.FieldStorage()

# get gene ids and construct gene ids query string
gene_ids = form['genes'].value
gene_ids_list = gene_ids.split()
sorted_gene_ids_list = sorted(gene_ids_list, key=int)
gene_ids_qstr = 'id=' + '&id='.join(sorted_gene_ids_list)
gene_ids_summary_qstr = 'id=' + ','.join(sorted_gene_ids_list)

# get terms and construct terms query string
terms = ''
if('terms' in form):
    terms = form['terms'].value
terms_list = terms.split()
sorted_terms_list = sorted(terms_list)
terms_qstr = ''
if('tiab_only' not in form):
    # search in full-text
    terms_qstr = 'term=' + '+OR+'.join(sorted_terms_list)
else:
    # search in title and abstract only
    sorted_terms_abstract_list = [s + "[TIAB]" for s in sorted_terms_list]
    terms_qstr = 'term=' + '+OR+'.join(sorted_terms_abstract_list)

# eutils base url
tool_email = 'tool=pubmed_article_rec&email=pubmed_article_rec@yahoo.com'
eutils_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils'
esummary_url = eutils_url + '/esummary.fcgi?' + tool_email
elink_url = eutils_url + '/elink.fcgi?' + tool_email
efetch_url = eutils_url + '/efetch.fcgi?' + tool_email

# pubmed url
pubmed_url = 'http://www.ncbi.nlm.nih.gov/pubmed'

# get gene summary information
genes_summary = urllib2.urlopen(esummary_url+'&db=gene&'+gene_ids_summary_qstr)
genes_summary_str = genes_summary.read()
genes_all_articles = urllib2.urlopen(elink_url+'&dbfrom=gene&db=pubmed&linkname=gene_pubmed&'+gene_ids_qstr)
genes_all_articles_str = genes_all_articles.read()

# parse gene summary information xml
genes_summary_root = ET.fromstring(genes_summary_str)
genes_all_articles_root = ET.fromstring(genes_all_articles_str)

# get abstract information
abstract_fixed_param = '&dbfrom=gene&db=pubmed&linkname=gene_pubmed&cmd=neighbor_history'
abstract_elink_qstr = elink_url+abstract_fixed_param+'&'+gene_ids_qstr+'&'+terms_qstr
abstract_elink = urllib2.urlopen(abstract_elink_qstr)
abstract_elink_str = abstract_elink.read()
time.sleep(1)

# parse abstract elink information xml
abstract_elink_root = ET.fromstring(abstract_elink_str)
geneid_webenv_key = dict()
for i in xrange(len(abstract_elink_root)):
    gene_id = abstract_elink_root[i].find('.//Id').text
    query_key_element = abstract_elink_root[i].find('.//QueryKey')
    query_key = None
    if(query_key_element != None):
        query_key = abstract_elink_root[i].find('.//QueryKey').text
    webenv = abstract_elink_root[i].find('./WebEnv').text
    geneid_webenv_key[gene_id] = (webenv, query_key)

# fixed params for efetch
efetch_fixed_param = '&db=pubmed&rettype=abstract&retmode=xml'

##################################### HTML ####################################

# header
print 'Content-Type: text/html Connection: keep-alive;charset=utf-8\r\n\r\n'
print '<html>'

# head
print '<head>'
print '<title>Search Summary</title>'
print """
        <script src="http://code.jquery.com/jquery-1.10.1.min.js"></script>
        <script src="http://code.jquery.com/jquery-migrate-1.2.1.min.js"></script>
        <script src="../javascripts/sorttable.js"></script>
      """
print """
<script type="text/javascript">
    $(document).ready(function(){
        $('.global_opt').change(function() {
            var id_str = $(this).attr('id');
            var id_arr = id_str.split('_');
            var term = id_arr[2];
            alert(term);
        });
        $('.local_opt').change(function() {
            var id_str = $(this).attr('id');
            var id_arr = id_str.split('_');
            var gene_id = id_arr[4];
            var table_id = "articles_gene_id_"+gene_id;
            sorttable(table_id);
        });
        $('.show_more_less').live('click', function() {
            $div = $(this).parent();
            if ($div.data('open')) {
                $div.css({height:'75px', overflow:'hidden', 'line-height':'25px'});
                $div.data('open', 0);
                $(this).html("show more");
            }
            else {
                $div.css({height:'100%'});
                $div.data('open', 1);
                $(this).html("show less");
            }
        });
        $('.show_hide').live('click', function() {
            if ($(this).html() == "hide") {
                $(this).siblings("table").css({display:'none'});
                $(this).html("show");
            }
            else {
                $(this).siblings("table").css({display:'inherit'});
                $(this).html("hide");
            }
        });
    });
</script>
"""
print """
        <style>
            table {
                border-bottom: 1px Solid Black;         
                border-right: 1px Solid Black;         
                border-collapse : collapse;  
            }
            table td, table th {    
                border-left: 1px Solid Black;         
                border-top: 1px Solid Black;              
                border-bottom: 1px Solid Black;    
                border-right:none;  
            }
            table td {
                max-width:700px;
                word-wrap:break-word;
            }
            .abstract_txt {
                line-height:25px;
                height:75px;
                overflow: hidden;
            }
        </style>
        """
print """<script src="../tablesorter/jquery-latest.js"></script>"""
print '</head>'

# body
print '<body>'
print '<h2><a id="top">Search Result</a></h2>'
print '<a href="/">Make Another Search</a><br/><br/>'

# create navigation
print '<a style="font-weight:bold" id="nav">Navigation by Gene ID</a><br/>'
for gene_id in sorted_gene_ids_list:
    print '<a href="#gene_id_%s">%s</a>' % (gene_id, gene_id)
print '<br/><br/><hr/>'

# create sort options
"""
print '<a id="global_sort">Sort Articles by Term Occurence</a><br/>'
print '<a>applies to all Article Summary tables</a><br/>'
print '<a>choose multiple terms to sort by the sum of their occurence in abstracts</a><br/>'
for term in sorted_terms_list:
    print '<input type="checkbox" id="global_opt_%s" class="global_opt">%s' % (term, term)

print '<br/><br/><hr/>'
"""

# flush the stdout buffer
sys.stdout.flush()

# iterate though gene ids
for i in xrange(len(sorted_gene_ids_list)):
    # parse gene summary xml
    gene_info = genes_summary_root[i]
    gene_all_articles_info = genes_all_articles_root[i]
    gene_id = gene_info.find('Id').text
    gene_name = gene_info.find("./Item[@Name='Name']").text
    gene_descp = gene_info.find("./Item[@Name='Description']").text
    gene_aliases = gene_info.find("./Item[@Name='OtherAliases']").text
    gene_desg = gene_info.find("./Item[@Name='OtherDesignations']").text
    gene_ns = gene_info.find("./Item[@Name='NomenclatureSymbol']").text
    gene_nn = gene_info.find("./Item[@Name='NomenclatureName']").text
    gene_summary = gene_info.find("./Item[@Name='Summary']").text
    gene_url = 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=Retrieve&dopt=full_report&list_uids='+gene_id
    gene_article_cnt = len(gene_all_articles_info.findall("./LinkSetDb/Link"))
    rna_url = 'http://genomernai.de/GenomeRNAi/genedetails/'+gene_id
    biogps_url = 'http://biogps.org/#goto=genereport&id='+gene_id
    genecards_url = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene='+gene_name
    
    # create bookmark
    print '<a id="gene_id_%s" href="#top">Return to Top</a><br/><br/>' % gene_id
    
    # title
    print '<b>Info and Related Abstracts for Gene %s</b><br/><br/>' % gene_id
    
    # print gene summary information table
    print '<div>'
    print '<a style="font-weight:bold">Gene Info</a>'
    print '<button class="show_hide" type="button">hide</button>'
    print '<br/>'
    print '<table class="info_table">'
    print '<tr><td>Id</td><td>%s</td></tr>' % gene_id
    print '<tr><td>Name</td><td>%s</td></tr>' % xstr(gene_name)
    print '<tr><td>Description</td><td>%s</td></tr>' % xstr(gene_descp)
    print '<tr><td>Aliases</td><td>%s</td></tr>' % xstr(gene_aliases)
    print '<tr><td>Other Designations</td><td>%s</td></tr>' % xstr(gene_desg)
    print '<tr><td>Summary</td><td>%s</td></tr>' % xstr(gene_summary)
    print '<tr><td>Gene URL</td><td><a href="%s" target="_blank">%s</a></td></tr>' % (gene_url, gene_url)
    print '<tr><td>RNA URL</td><td><a href="%s" target="_blank">%s</a></td></tr>' % (rna_url, rna_url)
    print '<tr><td>BioGPS URL</td><td><a href="%s" target="_blank">%s</a></td></tr>' % (biogps_url, biogps_url)
    print '<tr><td>GeneCards URL</td><td><a href="%s" target="_blank">%s</a></td></tr>' % (genecards_url, genecards_url)
    print '<tr><td>Total PMID Count</td><td>%d</td></tr>' % gene_article_cnt
    print '</table>'
    print '</div>'
    print '<br/>'
     
    # flush the stdout buffer
    sys.stdout.flush() 
     
    # get article abstract from eutil server
    webenv = geneid_webenv_key[gene_id][0]
    query_key = geneid_webenv_key[gene_id][1]
    abstract_efetch_qstr = None
    abstract_efetch_root = []
    if(query_key != None):
        abstract_efetch_qstr = efetch_url+efetch_fixed_param+'&WebEnv='+webenv+'&query_key='+query_key
        abstract_efetch = urllib2.urlopen(abstract_efetch_qstr)
        abstract_efetch_str = abstract_efetch.read()
        abstract_efetch_root = ET.fromstring(abstract_efetch_str)
    num_articles = len(abstract_efetch_root)
    
    # group abstract by term
    term_abstract = dict()
    for term in sorted_terms_list:
        if(term not in term_abstract):
            term_abstract[term] = []
        for j in xrange(num_articles):
            pmid = abstract_efetch_root[j].find('.//PMID').text
            article_url = pubmed_url+'/'+pmid
            title = abstract_efetch_root[j].find('.//ArticleTitle').text
            text_elements = abstract_efetch_root[j].findall('.//AbstractText')
            text = get_abstract_text(text_elements)
            if(text.lower().find(term) >= 0):
                article_obj = {'url': article_url, 'title': title, 'text': text}
                term_abstract[term].append(article_obj)
    
    # print abstract table grouped by term
    print '<div>'
    print '<b>Abstracts</b> (grouped by terms)'
    print '<button class="show_hide" type="button">hide</button>'
    print '<br/>'
    print '<table id="articles_by_term_gene_id_%s" class="info_table">' % xstr(gene_id)
    print '<tr><td>Term</td><td colspan="2">Title & Abstract</td></tr>'
    for term in sorted_terms_list:
        num_articles_for_term = len(term_abstract[term])
        if(num_articles_for_term > 0):
            print '<tr><td rowspan="%d">%s<br/>(%d in total)</td>' % (num_articles_for_term,term,num_articles_for_term)
        else:
            print '<tr><td>%s<br/>(%d in total)</td><td>0</td><td>N/A</td></tr>' % (term, num_articles_for_term)
        for j in xrange(num_articles_for_term):
            if(j != 0):
                print '<tr>'
            print '<td>%d</td>' % (j+1)
            print '<td>'
            print '<a href="%s" target="_blank">%s</a>' % (term_abstract[term][j]['url'],term_abstract[term][j]['title'])
            print '<div class="abstract_txt"><button class="show_more_less" type="button">show more</button>%s</div>' % term_abstract[term][j]['text']
            print '</td>'
            print '</tr>'
    print '</table>'
    print '</div>'
    
    print '<br/>'
    
    # print abstracts table
    print '<div>'
    print '<a style="font-weight:bold">Abstracts</a>'
    print '<a>(choose terms to sort by the sum of their occurences in titles and abstracts)'
    print '<button class="show_hide" type="button">hide</button>'
    print '<br/>'
    print '<table id="articles_gene_id_%s" class="info_table">' % xstr(gene_id)
    
    # print table header
    num_terms = len(terms_list)
    print '<tr>'
    print '<td></td><td>Title & Abstract (%d related articles in total)</td>' % num_articles
    for term in sorted_terms_list:
        print '<td><input type="checkbox" id="local_opt_gene_id_%s_%s" class="local_opt">%s</td>' % (gene_id,term,term)
    print '</tr>'
    
    # print table body
    for j in xrange(num_articles):
        pmid = abstract_efetch_root[j].find('.//PMID').text
        title = abstract_efetch_root[j].find('.//ArticleTitle').text
        text_elements = abstract_efetch_root[j].findall('.//AbstractText')
        text = get_abstract_text(text_elements)
        print '<tr>'
        print '<td>%d</td>' % (j+1)
        print '<td><a href="%s/%s" target="_blank">%s</a>' % (pubmed_url,pmid,title)
        print '<div class="abstract_txt"><button class="show_more_less" type="button">show more</button>%s</div></td>' % text
        term_count = get_term_count(terms_list, title+' '+text)
        for term in sorted_terms_list:
            print '<td><a class="term_count">%d</a><br/>%s</td>' % (term_count[term], term)
        print '</tr>'
    
    # end article summary table
    print '</table>'
    print '</div>'
    
    # print separator
    print '<hr/>'
    
    # flush the stdout buffer
    sys.stdout.flush()
    
    # sleep for 1 second
    time.sleep(1)
    
# end html
print '</body>'
print '</html>'