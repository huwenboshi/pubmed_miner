#!/usr/bin/python -u

import xml.etree.ElementTree as ET
import urllib2
import cgitb
import time
import cgi
import sys
import codecs
import os

"""
sys.path.insert(0,'/UCSC/Panel-Auxiliary/mysql_python/MySQL-python-1.2.3')
sys.path.insert(0,'/UCSC/Panel-Auxiliary/mysql_python/MySQL-python-1.2.3/build/lib.linux-x86_64-2.7')
sys.path.insert(0,'/UCSC/Panel-Auxiliary/mysql_python/MySQL-python-1.2.3/build/lib.linux-x86_64-2.7/_mysql.so')
"""

from utils import *
from expand_utils import *
from consts import *

sys.stdout = codecs.getwriter('utf-8')(sys.stdout)
cgitb.enable()

con = connect_db()

########################### HTTP HTML HEADER ###################################

# print http_header and html header
print http_header
print

form = cgi.FieldStorage()
terms = form['terms'].value
terms_list = terms.split()
expanded_terms = expand_terms(con, terms_list)

for key in expanded_terms:
    val = expanded_terms[key]
    for t in val:
        print key+'\t'+t
