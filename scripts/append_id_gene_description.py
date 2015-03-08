#!/usr/bin/env python
# encoding: utf-8

import sys

def main(argv):
    f = open(argv[0])

    lines = (line.strip() for line in f.readlines())

    for line in lines:
        # skip comments
        if line.startswith('#'):
            print(line)
            continue

        l = line.split('\t')

        ensembl_gene_id = l[14]
        HGNC_symbol = l[3]
#        gene_description = l[-1]
        transcripts = l[-2]

#        if len(gene_description) > 0 and gene_description != '#NA':
#            gene_description = '%s:%s' % (HGNC_symbol, gene_description)
#            l[-1] = gene_description

        if len(transcripts) > 0 and transcripts != '#NA':
            transcripts = transcripts.replace(ensembl_gene_id, HGNC_symbol, 1)
            l[-2] = transcripts

        print('\t'.join(l))

if __name__ == '__main__':
    main(sys.argv[1:])
