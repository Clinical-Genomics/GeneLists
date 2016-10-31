#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import yaml

from genelist import api
from genelist.utils import git

def get_display_name(target_list, target_panel):
    """ """

    list_dir = os.path.dirname(target_list)
    panel_names = dict((line.split(': ') for line in open(os.path.join(list_dir, 'LISTS'), 'r')))

    return panel_names[target_panel].rstrip()


def main(cml, target_list, target_panel):
    comments, cml_lines = api.readlist(open(cml, 'r'))
    comments, target_lines = api.readlist(open(target_list, 'r'))

    cml_hgnc_symbols = set()
    for line in cml_lines:
        panels = line['Clinical_db_gene_annotation'].split(',')
        if target_panel in panels:
            cml_hgnc_symbols.update({ line['HGNC_symbol'] })

    target_hgnc_symbols = set()
    for line in target_lines:
        panels = line['Clinical_db_gene_annotation'].split(',')
        if target_panel in panels:
            target_hgnc_symbols.update({ line['HGNC_symbol'] })

    missing_hgnc_symbols = target_hgnc_symbols - cml_hgnc_symbols

    output = {
        target_panel: {
            'display': get_display_name(target_list, target_panel),
            'version': git.getgittag(target_list),
            'date': git.getgitlastmoddate(target_list),
            'genes': list(missing_hgnc_symbols)
        }
    }

    print(yaml.dump(output, default_flow_style=False))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
