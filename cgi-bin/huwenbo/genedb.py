#!/usr/bin/python

import sqlite3 as lite
import sys


# connect to database
def connect_db(db_name):
    con = None
    try:
        con = lite.connect(db_name)
        con.text_factory = str
        return con
    except lite.Error, e:
        return None
