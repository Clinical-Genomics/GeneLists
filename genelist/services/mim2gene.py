"""Basic interface to mim2gene. Requires a mim2gene file which can be optionally downloaded.
"""

import tempfile
from os.path import exists
from urllib.request import urlretrieve

class Mim2gene(object):
    """Basic interface to mim2gene. Requires a mim2gene file which can be optionally downloaded.

    Args:
        filename (str, None): Location of the mim2gene.txt file. If none given, it will be
                             downloaded from omim.org to a temporary file.
        download (bool, True): Force download, even if file exists.
    """

    def __init__(self, filename=None, download=False):

        if download or \
          filename is None or \
          not exists(filename):
            filename = self.download(filename)

        self.filename = filename

        # init
        self.omim_of = {} # HGNC_symbol: OMIM_id
        self.symbol_of = {} # OMIM_id: HGNC_symbol
        self.ensembl_gene_id_of = {} # OMIM_id: ensEMBL_gene_id
        self.type_of = {} # HGNC_symbol: OMIM_type

        # read in the mim2gene file
        if filename:
            self.read(self.filename)

    def download(self, filename=None):

        if filename is None:
            filename = tempfile.NamedTemporaryFile().name

        urlretrieve('http://omim.org/static/omim/data/mim2gene.txt', filename)

        return filename

    def read(self, filename):
        """Read in the mim2gene file and store it as a dict of OMIM id: HGNC_symbol.
        Only gene and gene/phenotype types will be saved.

        Kwargs:
                filename (str): the aboslute path to the mim2gene.txt file

        Returns: None
        """

        mim2gene_fh = open(filename, 'r')
        lines = (line.rstrip('\n') for line in mim2gene_fh)
        for line in lines:
            if line.startswith('#'):
                continue
            (file_omim_id, omim_type, gene_id, hgnc_symbol, ensembl) = line.split("\t")
            ensembl_gene_id = ensembl.split(',')[0]
            self.omim_of[hgnc_symbol] = file_omim_id
            self.symbol_of[file_omim_id] = hgnc_symbol if hgnc_symbol else False
            self.ensembl_gene_id_of[file_omim_id] = ensembl_gene_id if ensembl_gene_id else False
            self.type_of[hgnc_symbol] = omim_type
            # this works because omim_id != hgnc_symbol, for all hgnc_symbol
            self.type_of[file_omim_id] = omim_type

    def get_hgnc(self, omim_id, only_gene=True):
        """Looks up the omim_id in the mim2gene.txt file.
        If found and the omim type is 'gene', return the official HGNC symbol

        Args:
            omim_id (int): the omim id

        Returns:
            on omim morbid match, official HGNC symbol if type of gene or
            gene/phenotype otherwise False
        """

        if omim_id in self.symbol_of:
            hgnc_type = self.type_of[omim_id]
            if not only_gene or hgnc_type in ('gene', 'gene/phenotype'):
                return self.symbol_of[omim_id]
        return False

    def get_ensembl(self, omim_id):
        """Looks up the EnsEMBL gene id in the mim2gene.txt file.
        If found and the omim type is 'gene', return it

        Args:
            omim_id (int): the omim id

        Returns: on omim id match, EnsEMBL gene id if type of gene or gene/phenotype otherwise False
        """

        if omim_id in self.ensembl_gene_id_of and \
           self.ensembl_gene_id_of[omim_id] != None and \
           self.ensembl_gene_id_of[omim_id] != '-':
            return self.ensembl_gene_id_of[omim_id]
        return False

    def is_gene(self, hgnc_symbol):
        if hgnc_symbol in self.type_of and \
           self.type_of[hgnc_symbol] in ('gene', 'gene/phenotype'):
               return True
        return False
