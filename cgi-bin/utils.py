import urllib
import urllib2
import xml.etree.ElementTree as ET
import MySQLdb
from math import sqrt

################################# DB STUFF #####################################

# connect to database
def connect_db():
    con = None
    try:
        """
            con = MySQLdb.connect(host="127.0.0.1", user="UCLApanelRW",
                                  passwd="PanelWrite", db="UCLApanel")
        """
        con = MySQLdb.connect(host="localhost", user="huwenbo",
                              passwd="goldandblue", db="pubmed_miner")
        return con
    except:
        return None

# get rows from database
def fetch_from_db(cur):
    table = []
    for row in cur.fetchall():
        col = []
        for i in xrange(len(row)):
            col.append(row[i])
        table.append(col)
    return table

########################### CONTAINER STUFF ####################################

# get length safely
def safe_len(container):
    if(container == None):
        return 0
    return len(container)

# get dictionary value safely
def safe_getval(container, key):
    if(key in container):
        return container[key]
    return 0

########################### ID CONVERSION STUFF ################################

# convert gene symbol to entrez id
def symbol2entrez(gene_symbols_list):
    
    # get conversion data from biomart
    xml_query = """<?xml version="1.0" encoding="UTF-8"?> \
        <!DOCTYPE Query> \
        <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" \
         uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >	\
	        <Dataset name = "hsapiens_gene_ensembl" interface = "default" > \
		        <Filter name = "hgnc_symbol" value = "%s"/> \
		        <Attribute name = "hgnc_symbol" /> \
		        <Attribute name = "entrezgene" /> \
	        </Dataset> \
        </Query>""" % ','.join(gene_symbols_list)
    url = 'http://www.ensembl.org/biomart/martservice?'
    values = {'query' : xml_query}
    data = urllib.urlencode(values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    the_page = response.read()
    
    # parse gene id mapping
    sym_id = dict()
    lines = the_page.split('\n')
    for line in lines:
        if(len(line) > 0):
            cols = line.split('\t')
            if(len(cols) > 1):
                gene_sym = cols[0]
                entrez_id = cols[1]
                sym_id[gene_sym] = entrez_id
    
    return sym_id

# convert entrez to gene symbol
def entrez2symbol(gene_ids_list):
    
    # get conversion data from biomart
    xml_query = """<?xml version="1.0" encoding="UTF-8"?> \
        <!DOCTYPE Query> \
        <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" \
         uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >	\
	        <Dataset name = "hsapiens_gene_ensembl" interface = "default" > \
		        <Filter name = "entrezgene" value = "%s"/> \
		        <Attribute name = "entrezgene" /> \
		        <Attribute name = "hgnc_symbol" /> \
	        </Dataset> \
        </Query>""" % ','.join(gene_ids_list)
    url = 'http://www.ensembl.org/biomart/martservice?'
    values = {'query' : xml_query}
    data = urllib.urlencode(values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    the_page = response.read()
    
    # parse gene id mapping
    id_sym = dict()
    lines = the_page.split('\n')
    for line in lines:
        if(len(line) > 0):
            cols = line.split('\t')
            if(len(cols) > 1):
                entrez_id = cols[0]
                gene_sym = cols[1]
                id_sym[entrez_id] = gene_sym
    
    return id_sym

######################### CORRELATION ##########################################

def pearson(dataseries1, dataseries2):
    result = 0.0
    sum_sq_x = 0.0
    sum_sq_y = 0.0
    sum_coproduct = 0.0
    mean_x = dataseries1[0]
    mean_y = dataseries2[0]

    for i in range(2,len(dataseries1)+1):
        sweep = (i-1)/float(i)
        delta_x = dataseries1[i-1]-mean_x
        delta_y = dataseries2[i-1]-mean_y
        sum_sq_x += delta_x * delta_x * sweep
        sum_sq_y += delta_y * delta_y * sweep
        sum_coproduct += delta_x * delta_y * sweep
        mean_x += delta_x / float(i)
        mean_y += delta_y / float(i)

    pop_sd_x = (sum_sq_x / float(len(dataseries1)))**0.5
    pop_sd_y = (sum_sq_y / float(len(dataseries1)))**0.5
    cov_x_y = sum_coproduct / float(len(dataseries1))

    if(pop_sd_x*pop_sd_y == 0):
        return 0
        
    result = cov_x_y / (pop_sd_x*pop_sd_y)

    return result
