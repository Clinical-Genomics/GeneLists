#!/usr/bin/env python
# encoding: utf-8

# TODO load the acronyms from a folder provided through a config file

from __future__ import print_function
import re
import sys
import glob

class Acronyms(object):

    def __init__(self, base_dir):
        """Reads in the LISTS files from the git repo's to get gene list and panel descriptions.
        Reads in all gene lists to gather all the panels for each gene list. Will fail silently
        when no panels are found or on faulty lines in the gene list.

        Args:
            base_dir (str): full path to the basedir of the repo's
        """

        self.acronyms = {} # acro: long text
        self.panels   = {} # acro: list of panels (acro)

        for list_file in glob.glob(base_dir + '/*/LISTS'):
            f = open(list_file)
            for line in f.readlines():
                acro, sep, full_name = line.strip().partition(': ')
                self.acronyms[acro] = full_name

        re_gl_name = re.compile(r'cust...-(.*).txt')
        for gl in glob.glob(base_dir + '/*/cust*txt'):
            gl_name = re.search(re_gl_name, gl).group(1)

            with open(gl, 'r') as f:
                lines = (line.strip() for line in f.readlines())
                panels = {}

                for line in lines:
                    # skip comments
                    if line.startswith('#'):
                        continue

                    l = line.split('\t')

                    if len(l) < 18: # skip faulty lines silently
                        continue

                    cur_panels = l[17].split(',')
                    cur_panels = ( panel.strip() for panel in cur_panels )
                    for cur_panel in cur_panels:
                        panels[cur_panel] = 1

                self.panels[gl_name] = panels.keys()

    def __getitem__(self,  acronym):
        """Retrieve the acronym description

        Args:
            acronym (str): a gene list/panel acronym

        Returns: acronym description

        """
        return self.acronyms.get(acronym, '')

    def get_panels_of(self, acronym):
        """Retrieve the panels of a gene list

        Args:
            acronym (str): a gene list acronym.

        Returns: a list of panels (acronyms)
        """
        return self.panels.get(acronym, [])

def main(args):
    acrs = Acronyms(args[0])
    print(acrs['NMD'])

if __name__ == '__main__':
    main(sys.argv[1:])
