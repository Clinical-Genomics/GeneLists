#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import click

from .modules.genelist import Genelist
from .modules.sanity import Sanity 
from .modules.panels import get_panels
from .modules.merge import merge_panels

#logger = logging.getLogger(__name__)

__version__ = '1.20.0'

@click.group()
def run():
    """Handy dandy package to automatically annotate and validate gene lists"""
    pass

@run.command()
@click.argument('infile', nargs=1, type=click.Path(exists=True))
@click.argument('outfile', nargs=1, type=click.Path())
@click.option('--zero', is_flag=True, default=False, show_default=True,
              help='Will convert 0-based coordinates to 1-based.')
@click.option('--verbose', is_flag=True, default=False, show_default=True,
              help='Will show conflict messages from EnsEMBLdb inbetween the gene list lines.')
@click.option('--errors', is_flag=True, default=False, show_default=True,
              help='Will not output the gene list, but only error messages.')
@click.option('--download-mim2gene', is_flag=True, default=False, show_default=True,
              help='Will download a new version of mim2gene.txt, used to check the OMIM type.')
@click.option('--mim2gene', is_flag=True, default=False, show_default=True,
              help='Will try to resolve an HGNC symbol with the help of mim2gene.txt.')
def fetch(infile, outfile, zero, verbose, errors, download_mim2gene, mim2gene):
    """Fetch all annotations."""

    genelist = Genelist()
    genelist.annotate(infile, outfile, zero=zero, verbose=verbose, errors=
                      errors,download_mim2gene=download_mim2gene, mim2gene=mim2gene)

@run.command()
@click.argument('genelist', nargs=1, type=click.Path(exists=True))
def validate(genelist):
    """Validate a genelist. Will print out messages as to what is wrong."""

    sanity = Sanity()
    sanity.check(genelist)

@run.command()
@click.argument('genelist', nargs=1, type=click.Path(exists=True))
def panels(genelist):
    """List panels in a gene list """

    print('\n'.join(get_panels(genelist)))

@run.command()
@click.argument('infiles', nargs=-1, required=True, type=click.File('r'))
@click.option('--database', '-d', multiple=True, help='only take HGNC_symbols from this database.')
def merge(infiles, database):
    """ Merge gene lists. Will only output HGNC_symbol, EnsEMBL_gene_id and Database columns.

    Args:
        infiles: paths to gene lists.
    """
    merge_panels(infiles, database)

def setup_logging(level='INFO'):
    """Setup the loggin for this package

    Args:
        level (str): log level.

    """

    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # customize formatter, align each column
    template = "[%(asctime)s] %(name)-25s %(levelname)-8s %(message)s"
    formatter = logging.Formatter(template)

    # add a basic STDERR handler to the logger
    console = logging.StreamHandler()
    console.setLevel(level)
    console.setFormatter(formatter)

    root_logger.addHandler(console)
    return root_logger

if __name__ == "__main__":
    setup_logging(level='DEBUG')
    #logger.info('Version: %s %s', __file__, __version__)
    run()
