#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import argparse
import subprocess
import os
import re
from datetime import datetime

from .git import getgitlastmoddate, getgittag
from .acronyms import Acronyms

def main(argv):
    parser = argparse.ArgumentParser(description='Update clinicalgenomics.se with gene list names and versions')
    parser.add_argument('infiles', nargs="+", type=argparse.FileType('r'), help='')
    args = parser.parse_args(argv)

    acronyms = Acronyms(os.path.dirname(os.path.dirname(os.path.realpath(args.infiles[0].name))))
    re_gl_name = re.compile(r'cust...-(.*).txt')

    versions = {} # Database => { 'Version': Version, 'Datum': Date, 'Beskrivning': description, 'Databas': database acronym, 'Panels': panels}
    for infile in args.infiles:

        # get the database name from the filename
        match = re.search(re_gl_name, infile.name)
        database = match.group(1)

        # fill versions dict
        mod_date = getgitlastmoddate(infile.name, '%Y-%m-%d %H:%M:%S')
        version = getgittag(infile.name, date=mod_date) # get version on that date
        full_name = acronyms[database]
        panels = acronyms.get_panels_of(database)
        versions[database] = {
            'Version': version,
            'Datum': mod_date.partition(' ')[0],
            'Beskrivning': full_name,
            'Databas': database,
            'Panels': panels
        }

    print("""---
layout: base
permalink: /namnpagenlistor/
title: Namn på genlistor
---

# Namn på genlistor
""")

    keys = ['Databas', 'Beskrivning', 'Version', 'Datum']
    for database in sorted(versions.keys()):
        print('|%s|' % '|'.join(keys))
        print('|%s' % ('---|' * len(keys)))
        print('|%s|' % '|'.join([ versions[database][key] for key in keys ]))

        if len(versions[database]['Panels']) > 1:
            print()
            print('Panels:')
            for panel in versions[database]['Panels']:
                descr = acronyms[panel]
                panel_descr = panel if len(descr) == 0 else ': '.join([panel, descr])
                print('- %s' % panel_descr)
        print()

if __name__ == '__main__':
    main(sys.argv[1:])
