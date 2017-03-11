#!/usr/bin/env python2


import re
import os
import os.path
import sys
import json
#from pprint import pprint

compilation_database_pattern = re.compile('(?<=\s)-[DIOUWfgs][^=\s]+(?:=\\"[^"]+\\"|=[^"]\S+)?')

dbfile_search_patterns = [r'/build/compile_commands.json',
                          r'/debug/compile_commands.json',
                          r'/release/compile_commands.json']

def find_dffile(path):

    for dbfile_pat in dbfile_search_patterns:
        dbfile = path + dbfile_pat
        if os.path.isfile(dbfile):
            return dbfile

    # not found
    newpath = os.path.dirname(path)
    if path == newpath:
        raise Exception(r'Unable to find ' + dbfile_pat)
    return find_dffile(newpath)

def load_db(filename):
    compilation_database = {}
    with open(filename) as compilation_database_file:
        compilation_database_entries = json.load(compilation_database_file)

    entry = 0
    for compilation_entry in compilation_database_entries:
        entry = entry + 1
        compilation_database[compilation_entry["file"]] = [ p.strip() for p in compilation_database_pattern.findall(compilation_entry["command"]) ]
    return compilation_database

def guess_option(filename, db):
    guessed = []

    name, ext = os.path.splitext(filename)
    guesses = [name + '.cpp', name + '.c']

    for key in db.iterkeys():
        if key.endswith('main.cpp'):
            guesses.append(key)
            break

    for g in guesses:
        if os.path.isfile(g):
            # print ('guessed = ', g)
            for opt in db[g]:
                guessed.append(opt)
            break

    return guessed



srcfile = unicode(sys.argv[1],"utf-8")
dbfile = find_dffile(os.path.dirname(srcfile))
dbpath = os.path.dirname(dbfile)

db = load_db(dbfile)
srcfile = unicode(sys.argv[2],"utf-8")
if db:
    if srcfile in db:
        for option in db[srcfile]:
            print(option)
    else:
        for option in guess_option(srcfile, db):
            print(option)