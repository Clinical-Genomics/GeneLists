#!/usr/bin/env python
# encoding: utf-8

def readlist(lines):
    """ """

    # slurp and make a line
    raw_data = (line.strip() for line in lines) # sluuuurp
    parsable_data = (line.split("\t") for line in raw_data)

    # skip parsing of leading comments
    comments = []
    line = next(parsable_data)
    while line[0].startswith('##'):
        if not line[0].startswith('##contig'):
            comments.append(line) # skip all contig comments
        line = next(parsable_data)

    # list to dict
    header = line # get the header
    if header[0].startswith('#'):
        header[0] = header[0].lstrip('#')

    dict_data = ( dict(zip(header, line)) for line in parsable_data )

    return comments, dict_data
