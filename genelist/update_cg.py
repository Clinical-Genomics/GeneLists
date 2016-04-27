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

def panels_2_html(panels, acronyms):
    """Converts a list of panel acronyms to a HTML list.
    Will add the description of the panel to the list as well.

    Args:
        panels (list): a list of panel acronyms

    Returns: a HTML formatted string.

    """
    if len(panels) == 0:
        return ''

    li_panels = []
    i=0
    for panel in panels:
        i += 1
        descr = acronyms[panel]
        panel_descr = panel if len(descr) == 0 else ': '.join([panel, descr])
        li_panels.append('%d. %s<br />' % (i, panel_descr))

    return ''.join(li_panels)

def main(argv):
    parser = argparse.ArgumentParser(description='Update clinicalgenomics.se with gene list names and versions')
    parser.add_argument('infiles', nargs="+", type=argparse.FileType('r'), help='')
    args = parser.parse_args(argv)

    acronyms = Acronyms(os.path.dirname(os.path.dirname(os.path.realpath(args.infiles[0].name))))
    re_gl_name = re.compile(r'cust...-(.*).txt')

    versions = {} # Database => { 'Version': Version, 'Datum': Date, 'Beskrivning': description, 'Databas': database acronym}
    for infile in args.infiles:

        # get the database name from the filename
        match = re.search(re_gl_name, infile.name)
        database = match.group(1)

        # fill versions dict
        mod_date = getgitlastmoddate(infile.name, '%Y-%m-%d %H:%M:%S')
        if not mod_date:
            print('WARNING: {} not committed!'.format(infile.name))
            continue
        version = getgittag(infile.name, date=mod_date) # get version on that date
        if '-' in version: # this is a new tag with the GL name included
            version = version.split('-')[1]
        full_name = acronyms[database]
        panels = acronyms.get_panels_of(database)
        panels = [ panel for panel in panels if panel != 'FullList' ]
        panels = full_name if len(panels) == 1 else panels_2_html(panels, acronyms)
        versions[database] = {
            'Version': version,
            'Datum': mod_date.partition(' ')[0],
            'Beskrivning': panels,
            'Databas': database,
        }

    print("""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <title>Clinical Genomics customer portal</title>
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta name="description" content="Public customer portal at SciLifeLab Clinical Genomics.
">
  <meta name="theme-color" content="#ef7c00">
  <link rel="canonical" href="http://www.clinicalgenomics.se/namnpagenlistor/">

  <!-- Fonts -->
  <link href="http://fonts.googleapis.com/css?family=Source+Sans+Pro:400,700,400italic" rel="stylesheet" type="text/css">

  <!-- Custom CSS -->
  <link rel="stylesheet" href="/assets/css/style.css">
</head>


<body>

  <aside class="cg-sidebar">
    <div class="cg-sidebar-header">
      <div class="cg-logo">
        <div class="cg-logo-icon">
          <img src="/assets/img/cg-logo.svg"
               alt="Clinical Genomics logo">
        </div>
        <div class="cg-logo-label">
          <a href="/" class="cg-logo-label-title">Clinical Genomics</a>
          <div class="cg-logo-label-subtitle">a SciLifeLab facility</div>
        </div>
      </div>

      <div class="cg-sheet ticket-sheet">
        <a href="https://clinical-scilifelab.supportsystem.com/open.php"
           class="cg-sheet-button">
          Öppna ny ticket
        </a>
        <a href="https://clinical-scilifelab.supportsystem.com/view.php"
           class="cg-sheet-button">
          Kontrollera ticketstatus
        </a>
      </div>
    </div>
    <div class="cg-sidebar-body">
      <div class="cg-list">
          <a href="/analyser/" class="cg-list-item ">
            Analyser
          </a>
          <a href="/bestallningar/" class="cg-list-item ">
            Beställningar
          </a>
          <a href="/dataleverans/" class="cg-list-item ">
            Dataleveranser
          </a>
          <a href="/kontakt/" class="cg-list-item ">
            Kontakt
          </a>
          <a href="/namnpagenlistor/" class="cg-list-item current">
            Namn på genlistor
          </a>
          <a href="/provkrav/" class="cg-list-item ">
            Krav på prover
          </a>
      </div>
    </div>
    <div class="cg-sidebar-footer">
      &copy; SciLifeLab 2015
    </div>
  </aside>

  <main class="cg-content">
    <div class="cg-main">
      <div class="cg-sidebar-toogle">Meny</div>

      <h1 id="namn-p-genlistor">Namn på genlistor</h1>
    """)

    print("""<table>
    <thead>
        <tr>
    """)

    keys = ['Databas', 'Beskrivning', 'Version', 'Datum']
    print('<th>%s</th>' % '</th><th>'.join(keys))

    print("""
        </tr>
    </thead>
    <tbody>
    """)
    for database in sorted(versions.keys()):
        print('<tr>')
        print('<th>%s</th>' % '</th><th>'.join([ versions[database][key] for key in keys ]))
        print('</tr>')

    print("""
  </tbody>
</table>

    </div>
  </main>

  <script src="//ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js"></script>
  <script src="/assets/js/main.js"></script>

</body>
</html>
    """)

if __name__ == '__main__':
    main(sys.argv[1:])
