#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys

def main(argv):
    f = open(argv[0])

    lines = (line.strip() for line in f.readlines())

    panels = {}

    for line in lines:
        # skip comments
        if line.startswith('#'):
            continue

        l = line.split('\t')

        cur_panels = l[17].split(',')
        cur_panels = ( panel.strip() for panel in cur_panels )
        for cur_panel in cur_panels:
            panels[cur_panel] = 1

    print('\n'.join(panels.keys()))

if __name__ == '__main__':
    main(sys.argv[1:])
