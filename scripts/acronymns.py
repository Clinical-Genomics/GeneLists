#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import glob

class Acronymns(object):

    def __init__(self, base_dir):
        """Reads in the LISTS files from the git repo's

        Args:
            base_dir (str): full path to the basedir of the repo's
        """
        self.acronymns = {} # acro: long text
        for list_file in glob.glob(base_dir + '/*/LISTS'):
            f = open(list_file)
            for line in f.readlines():
                acro, sep, full_name = line.strip().partition(': ')
                self.acronymns[acro] = full_name

    def __getitem__(self,  acronymn):
        """TODO: Docstring for __item__.

        Args:
            acronymn (TODO): TODO

        Returns: TODO

        """
        return self.acronymns.get(acronymn, None)

def main(args):
    acrs = Acronymns(args[0])
    print(acrs['NMD'])

if __name__ == '__main__':
    main(sys.argv[1:])
