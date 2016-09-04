""" Provide all annotation functionality in one place."""
# encoding: utf-8

from __future__ import print_function
import pymysql
import re
import os
import yaml
import logging
from io import StringIO
from collections import Counter

from ..services.omim import OMIM
from ..services.ensembl import Ensembl
from ..services.genenames import Genenames
from ..services.uniprot import Uniprot
from ..services.mim2gene import Mim2gene

def there(line, key):
    """ Checks if the key is in the line and has a value that doesn't resolve to False.

    Args:
        line (dict): a dict.
        key (str): a possible key in dict.

    return: Value if key exists and has a value that doesn't resolve to False otherwise False.
    """

    if key not in line:
        return False

    if not line[key]:
        return False

    return line[key]

class Fetch(object):

    """Provide all annotation functionality in one class. """

    def __init__(self, config, download_mim2gene):
        self.header = ['Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_symbol', 'Protein_name',
                          'Symptoms', 'Biochemistry', 'Imaging', 'Disease_trivial_name',
                          'Trivial_name_short',
                          'Phenotypic_disease_model', 'OMIM_morbid', 'Gene_locus', 'UniProt_id',
                          'Ensembl_gene_id', 'Ensemble_transcript_ID', 'Reduced_penetrance',
                          'Clinical_db_gene_annotation', 'Disease_associated_transcript',
                          'Ensembl_transcript_to_refseq_transcript', 'Gene_description',
                          'Genetic_disease_model', 'HGNC_RefSeq_NM', 'Uniprot_protein_name',
                          'Database_entry_version', 'Curator', 'Alias', 'Group_or_Pathway',
                          'Mosaicism', 'Comments']

        # columns that need a HGNC prefix
        self.prefix_header = ['Phenotypic_disease_model', 'OMIM_morbid', 'UniProt_id',
                              'Ensemble_transcript_ID', 'Disease_associated_transcript',
                              'Ensembl_transcript_to_refseq_transcript', 'Gene_description',
                              'HGNC_RefSeq_NM', 'Uniprot_protein_name']

        self.config = yaml.load(config)

        # EnsEMBL connection
        self.conn = pymysql.connect(
            host=self.config['ensembl']['host'],
            port=self.config['ensembl']['port'],
            user=self.config['ensembl']['user'],
            db=self.config['ensembl']['db'])

        self.logger = logging.getLogger(__name__)
        self.setup_logging(level='DEBUG')

        # check mem2gene.txt for HGNC symbol resolution
        self.mim2gene = self.init_mim2gene(download_mim2gene)

        self.reset()

    def reset(self):
        """ Reset state for a next genelist to annotate """

        self.delimiter = '|' # join elements of a field
        self.current_hgnc_id = ''
        self.current_line = {}
        self.line_nr = 0
        self.contigs = set()
        self.print_info = False
        self.print_warn = False
        self.print_error = False
        self.report_empty = False

        # reset the StringIO
        self.log_buffer.truncate(0)
        self.log_buffer.seek(0)

    def print_header(self, header=None):
        """
        Prints the contigs meta data headers.
        Prints the header.

        Args:
                header (list, optional): a list of strings
        Note:
                Prints to STDOUT
        """
        if not header:
            header = self.header
        print(self.get_header(header=header))

    def get_header(self, header=None):
        """
        Returns the header.

        Args:
                header (list, optional): a list of strings
        Return (str): the header
        """
        if not header:
            header = self.header

        return '#' + "\t".join(self.header)

    def get_contigs(self):
        """ Returns the contigs meta data headers.

        Yields (list): one item per header
        """

        for contig in sorted(self.contigs, key=lambda item:
                             (int(item) if item.isdigit() else float('inf'), item)):
            yield '##contig=<ID={}'.format(contig)

    def format_line(self, line):
        """Formats a line based on the order of the headers

        Args:
            line (dict): dict with values for one line of a gene list.
                         All values of headers should be present as keys.

        Returns (str): a tab-delim line
        """
        ordered_line = list()
        for column_name in self.header:
            ordered_line.append(str(line[column_name]))
        return "\t".join(ordered_line)

    def print_line(self, line):
        """Prints a line based on the order of the headers

        Args:
                line (dict): dict with values for one line of a gene list.
                             All values of headers should be present as keys

        Returns:
                pass

        Note:
                Will print the STDOUT
        """
        print(self.format_line(line))

    def warn(self, line, key=None):
        """print only if the verbose switch has been set

        Args:
                line (str): line to print to STDOUT
                key (str, optional): the key in current data.
                                     Used to check if the value behind the key is empty.

        Returns:
                pass
        """
        if self.print_warn:
            if key is None: # no key is given, report all warnings
                self.logger.warning(line, extra={'line_nr': self.line_nr,
                                    'hgnc_id': self.current_hgnc_id})
            elif self.report_empty or (key and
                key in self.current_line and self.current_line[key]):
                self.logger.warning(line, extra={'line_nr': self.line_nr,
                                    'hgnc_id': self.current_hgnc_id})

    def info(self, line):
        """print only if the verbose switch has been set

        Args:
                line (str): line to print to STDOUT

        Returns:
                pass
        """
        if self.print_info:
            self.logger.info(line, extra={'line_nr': self.line_nr, 'hgnc_id': self.current_hgnc_id})

    def error(self, line):
        """print only if the verbose switch has been set

        Args:
                line (str): line to print to STDOUT

        Returns:
                pass
        """
        if self.print_error:
            self.logger.error(line, extra={'line_nr': self.line_nr, 'hgnc_id': self.current_hgnc_id})

    def get_context(self, data):
        """Increments the global line_nr for each passing line.
        Gets the HGNC identifier so we can make more sensical warning messages.

        Args:
            lines (list of dicts): each dict contains a dict with header as keys
        Yields: line

        """
        for line in data:
            self.line_nr += 1
            self.current_hgnc_id = line['HGNC_symbol']
            self.current_line = line

            yield line

    def remove_hgnc_prefix(self, line):
        """ Removes the prefixed HGNC symbol from all fields

        Args
            line (dict): representing the line and columns in a gene list.

        Yields (dict):
            a line with the prefixed HGNC identifiers removed
        """

        for column in self.header:
            if column in line:
                line[column] = re.sub(r'^.*?:', '', str(line[column]))
        return line

    def prepend_hgnc(self, data):
        for line in data:
            for column in self.prefix_header:
                if column in line and line[column]:
                    line[column] = '%s:%s' % (line['HGNC_symbol'], line[column])

            yield line

    def gather_contig(self, data):
        """Aggregates the contigs so we can print them as headers

        Args:
            data (list of dicts): Inner dict represents a row in a gene list

        Yields:
            dict: with the added HGNC symbol prepended to the HGNC_symbol column.
        """
        for line in data:
            contig = line['Chromosome']
            self.contigs.add(contig)

            yield line

    def cleanup(self, data):
        """Will clean the data according to rules set by MIP
                # replace white space seperated comma's with just a comma
                # replace ; with comma
                # remove leading and trailing white space
                # remove trailing commas
                # collapse multiple commas
                # remove prefixed HGNC symbol

        Args:
                data (list of lists): Inner list represents a row in a gene list

        Yield:
                list: a cleaned up row
        """
        for line in data:
            for key, value in line.items():
                if isinstance(value, str):
                    if line[key] == '#NA':
                        line[key] = ''
                    else:
                        line[key] = re.sub(r'\s*[,;]\s*', ',', value.strip()) # rm whitespace
                        line[key] = re.sub(r',+', ',', line[key]) # collapse commas
                        line[key] = line[key].rstrip(',') # rm trailing commas
                elif value is False:
                    line[key] = ''
            line = self.remove_hgnc_prefix(line)
            yield line

    def list2dict(self, header, data):
        """Will convert each row in the data from a list to dict using the header list as keys.

        Args:
                header (list): A list containing the keys for the dict generation
                data (list of lists): Inner list represents a row in a gene list

        Yields:
                dict: the next dictified line
        """
        for line in data:
            yield dict(zip(header, line))

    def merge_line(self, line, client):
        """Will merge line with client.
           ens will take precedence over client. Changes will be reported.

        Args:
                line (dict): dict with new values.
                client (dict): with values coming from the gene list.

        Yields:
                dict: merged ens and client dict
        """

        for key, value in client.items():
            has_new_value = there(line, key)
            if key in ('Gene_start', 'Gene_stop'):
                # don't report start/stop mismatches
                continue
            elif not has_new_value:
                # skip if no new value
                continue
            else:
                if str(line[key]) != str(value):
                    self.warn("{}: line '{}' differs from client '{}'".\
                              format(key, line[key], value), key)

        merged = client.copy()
        merged.update(line)
        return merged

    def add_mim2gene_alias(self, data):
        """Looks up the most recent HGNC symbol for an HGNC alias in mim2gene.txt and
        prepends it to the HGNC_symbol column.
        Normally, the most recent symbol will have more chance to have a hit in EnsEMBLdb.
        Only use this function when mim2gene switch is active

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: with the added HGNC symbol prepended to the HGNC_symbol column.
        """
        for line in data:
            hgnc_symbols = line['HGNC_symbol'].split(',')
            if 'OMIM_morbid' in line:
                omim_id = line['OMIM_morbid']
                hgnc_symbol = self.mim2gene.resolve_gene(omim_id)
                ensembl_gene_id = self.mim2gene.resolve_ensembl_gene_id(omim_id)
                if hgnc_symbol != False and hgnc_symbol not in hgnc_symbols:
                    self.info("Add mim2gene HGNC symbol %s" % hgnc_symbol)
                    hgnc_symbols.insert(0, hgnc_symbol)
                if ensembl_gene_id != False and 'Ensembl_gene_id' in line.keys() \
                   and line['Ensembl_gene_id'] != ensembl_gene_id:
                    self.warn("morbidmap '{}' differs from local '{}'".\
                           format(line['Ensembl_gene_id'], ensembl_gene_id), 'Ensembl_gene_id')
            line['HGNC_symbol'] = ','.join(hgnc_symbols)
            yield line

    def init_mim2gene(self, download_mim2gene):
        mim2gene_filename = os.path.join(os.path.dirname(__file__), 'mim2gene.txt')
        if download_mim2gene:
            self.info('Downloading {} ... '.format(mim2gene_filename))
            return Mim2gene(filename=mim2gene_filename, download=True)
        else:
            return Mim2gene(filename=mim2gene_filename)

    def fill_from_mim2gene(self, data):
        """ Fill in HGNC symbol, OMIM id and ensembl_gene_id.

        Order of precedence, meaning first or following will determine the value of the others.
            - OMIM id
            - HGNC symbol
            - ensembl gene id

        Method will warn when any value is overwritten.

        Args:
            data (list of dicts): Inner dict represents a row in a gene list

        Yields:
            dict: with the hgnc symbol, omim id and ensembl_gene_id filled in.

        """
        for line in data:
            omim = there(line, 'OMIM_morbid')
            hgnc = there(line, 'HGNC_symbol')
            ensembl = there(line, 'Ensembl_gene_id')

            if omim:
                yield self.merge_line(
                    {
                        'HGNC_symbol': self.mim2gene.get_hgnc(omim),
                        'Ensembl_gene_id': self.mim2gene.get_ensembl(omim)
                    },
                    line
                )
                continue

            if hgnc:
                yield self.merge_line(
                    {
                        'OMIM_morbid': self.mim2gene.get_omim(hgnc),
                        'Ensembl_gene_id': self.mim2gene.get_ensembl(hgnc)
                    },
                    line
                )
                continue

            # let's skip EnsEMBL checking for now
            yield line

    def fill_from_ensembl(self, data):
        """ Fill in Gene_start, Gene_stop, Chromosome and HGNC_symbol

        Method will warn when any value is overwritten.

        Args:
            data (list of dicts): Inner dict represents a row in a gene list

        Yields:
            dict: with the Gene_start, Gene_stop, Chromome and HGNC_symbol filled in.

        """
        with Ensembl() as ensembldb:
            for line in data:
                ensembl_gene_id = there(line, 'Ensembl_gene_id')
                if ensembl_gene_id:
                    for ensembl_line in ensembldb.query(ensembl_gene_id):
                        yield self.merge_line(ensembl_line, line)

    def fill(self, data):
        """ Removes #NA's and fills in '' for missing values.

        Args:
                data (list of dicts): representing the lines and columns in a gene list

        Yields:
                dict: with all missing columns filled in with ''
        """
        defaults = dict((column_name, '') for column_name in self.header)
        for line in data:
            d = defaults.copy()
            d.update(line)
            yield d

    def get_log_messages(self):
        self.log_buffer_handler.flush()
        self.log_buffer.flush()

        return self.log_buffer.getvalue()

    def setup_logging(self, level='INFO'):
        """ Set up logging """
        root_logger = logging.getLogger(__name__) # only log this package
        root_logger.setLevel(level)
        
        # customize formatter, align each column
        template = "#%(line_nr)s [%(hgnc_id)s] %(message)s"
        formatter = logging.Formatter(template)
        fancy_formatter = logging.Formatter('\033[93m' + template + '\033[0m')
        
        # add a basic STDERR handler to the logger
        console = logging.StreamHandler()
        console.setLevel(level)
        console.setFormatter(fancy_formatter)
        root_logger.addHandler(console)

	# add a basic Memory handler so we can prepend them to the genelist
        self.log_buffer = StringIO()
        self.log_buffer_handler = logging.StreamHandler(self.log_buffer)
        self.log_buffer_handler.setLevel(level)
        self.log_buffer_handler.setFormatter(formatter)
        root_logger.addHandler(self.log_buffer_handler)

        return root_logger

    def annotate(self, lines, warn=False, error=False, info=False, report_empty=False, remove_non_genes=False):
        """ Annotate a gene list """

        self.reset()

        # make sure we print if we are asked to
        verbose = False
        if info:
            self.print_info = True
            verbose = True
        if warn:
            self.print_warn = True
            verbose = True
        if error:
            self.print_error = True
            verbose = True

        if report_empty:
            self.report_empty = True
            verbose = True

        # slurp and make a line
        raw_data = (line.strip() for line in lines) # sluuuurp
        parsable_data = (line.split("\t") for line in raw_data)

        # skip parsing of leading comments
        comments = []
        line = next(parsable_data)
        self.line_nr += 1
        while line[0].startswith('##'):
            if not line[0].startswith('##contig'):
                comments.append(line) # skip all contig comments
            self.line_nr += 1
            line = next(parsable_data)

        # list to dict
        header = line # get the header
        if header[0].startswith('#'):
            header[0] = header[0].lstrip('#')
        dict_data = self.list2dict(header, parsable_data)

        # get some context for error messages
        context_data = self.get_context(dict_data)

        # clean up the input
        clean_data = self.cleanup(context_data)

        # Get OMIM morbid number
        # Get E!
        # all from mim2gene, the magical file
        mim2gene_data = self.fill_from_mim2gene(clean_data)

        # fill in info from ensembl
        ensembl_data = self.fill_from_ensembl(mim2gene_data)

        # fill in missing values with ''
        completed_data = self.fill(ensembl_data)

        # at last, clean up the output
        cleaner_data = self.cleanup(completed_data)

        # prepend the HGNC symbol to some fields
        prefixed_data = self.prepend_hgnc(cleaner_data)

        # get all contigs
        final_data = self.gather_contig(prefixed_data)
        print_data = []
        for line in final_data:
            print(self.format_line(line))
            print_data.append(line)

        # print the errors and warnings
        if verbose:
            # split the lines for easier unit testing
            for line in self.get_log_messages().split('\n'):
                yield line

        # print the gene list
        for comment in comments:
            yield '\t'.join(comment)
        for line in self.get_contigs():
            yield line
        yield self.get_header()
        for line in print_data:
            yield self.format_line(line)
