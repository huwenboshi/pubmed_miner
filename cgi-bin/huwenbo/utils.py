# eutils base url
tool_email = 'tool=pubmed_article_rec&email=pubmed_article_rec@yahoo.com'
eutils_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils'
esummary_url = eutils_url + '/esummary.fcgi?' + tool_email
elink_url = eutils_url + '/elink.fcgi?' + tool_email
efetch_url = eutils_url + '/efetch.fcgi?' + tool_email

# pubmed url
pubmed_url = 'http://www.ncbi.nlm.nih.gov/pubmed'

# fixed param for abstract
abstract_fixed_param = '&dbfrom=gene&db=pubmed&linkname=gene_pubmed&cmd=neighbor_history'

# fixed params for efetch
efetch_fixed_param = '&db=pubmed&rettype=abstract&retmode=xml'

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