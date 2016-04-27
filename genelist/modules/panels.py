""" module to list panels in a gene list """
# encoding: utf-8

from __future__ import print_function

def get_panels(genelist):
    """ List panels in a gene list

    Args:
        genelist (str): path to the gene list.

    Returns (list): list of panel names.
    """

    f = open(genelist)

    lines = (line.strip() for line in f.readlines())

    panels = {}

    for line in lines:
        # skip comments
        if line.startswith('#'):
            continue

        l = line.split('\t')

        cur_panels = l[17].split(',')
        cur_panels = (panel.strip() for panel in cur_panels)
        for cur_panel in cur_panels:
            panels[cur_panel] = 1

    return panels.keys()
