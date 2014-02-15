import urllib
import urllib2
import xml.etree.ElementTree as ET

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
            entrez_id = cols[0]
            gene_sym = cols[1]
            id_sym[entrez_id] = gene_sym
    
    return id_sym
