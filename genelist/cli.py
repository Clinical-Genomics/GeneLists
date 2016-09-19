#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import click

from .modules.fetch import Fetch
from .modules.mans import Mans
from .modules.sanity import Sanity
from .modules.panels import get_panels
from .modules.merge import merge_panels

#logger = logging.getLogger(__name__)

__version__ = '3.3.3'

@click.group()
def run():
    """Handy dandy package to automatically annotate and validate gene lists"""
    pass

@run.command()
@click.argument('infile', nargs=1, type=click.File('r'))
@click.argument('outfile', nargs=1, type=click.File('w'))
@click.option('--remove-non-genes', is_flag=True, default=False, show_default=True,
              help='Removes non-genes. Based on mim2gene.')
@click.option('--warn', is_flag=True, default=False, show_default=True,
              help='Print minor conflicts.')
@click.option('--error', is_flag=True, default=False, show_default=True,
              help='Print severe conflicts.')
@click.option('--info', is_flag=True, default=False, show_default=True,
              help='Be more verbose.')
@click.option('--report-empty', is_flag=True, default=False, show_default=True,
              help='Report warnings from empty fields.')
@click.option('--download-mim2gene', is_flag=True, default=False, show_default=True,
              help='Will download a new version of mim2gene.txt, used to check the OMIM type.')
@click.option('--config', '-c', required=True, type=click.File('r'),
              help='YAML config file.')
def fetch(infile, outfile, remove_non_genes, warn, error, info, download_mim2gene, report_empty, config):
    """Fetch all annotations."""

    fetch = Fetch(config, download_mim2gene=download_mim2gene)

    for line in fetch.annotate(lines=infile, remove_non_genes=remove_non_genes, info=info, error=error, warn=warn, report_empty=report_empty):
        outfile.write(line + '\n')

@run.command()
@click.argument('infile', nargs=1, type=click.File('r'))
@click.argument('outfile', nargs=1, type=click.File('w'))
@click.option('--warn', is_flag=True, default=False, show_default=True,
              help='Print minor conflicts.')
@click.option('--error', is_flag=True, default=False, show_default=True,
              help='Print severe conflicts.')
@click.option('--info', is_flag=True, default=False, show_default=True,
              help='Be more verbose.')
@click.option('--report-empty', is_flag=True, default=False, show_default=True,
              help='Report warnings from empty fields.')
@click.option('--config', '-c', required=True, type=click.File('r'),
              help='YAML config file.')
def mans(infile, outfile, warn, error, info, report_empty, config):
    """Fetch omim annotations.
    Rerquires following headers:
    #Description\tGene Symbols\tomim_morbid\tCyto Location
    """

    mans = Mans(config)

    for line in mans.annotate(lines=infile, info=info, error=error, warn=warn, report_empty=report_empty):
        outfile.write(line + '\n')

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
