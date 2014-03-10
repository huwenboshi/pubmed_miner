import urllib2
import cgitb
import time
import cgi
import sys
import codecs
import os

from utils import *

# get doids
def get_term_doid(con, term):
    if(con == None):
        return []
    cur = con.cursor()
    query = """select * from term_doid 
        where term like '%["""+term+"""]%'"""
    cur.execute(query)
    return fetch_from_db(cur)

# expand term
def expand_terms(con, terms):
    expanded_terms = dict()
    for term in terms:
        expanded_terms[term] = set()
        result = get_term_doid(con, term)
        for row in result:
            term_exp = row[0].replace(']',' ').replace('[','').lower()
            expanded_terms[term].add(term_exp)
    return expanded_terms
