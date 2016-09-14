""" Provide Mans' OMIM functionality """
# encoding: utf-8

from __future__ import print_function
import sys
import re
import os
import yaml
import logging
from io import StringIO

from ..services.omim import OMIM
from ..services.genenames import Genenames

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

class Mans(object):
    """Provide omim annotation functionality in one class. """

    def __init__(self, config):
        self.header = ['omim_morbid', 'omim_phenotype', 'description', 'inheritance_models']

        self.config = yaml.load(config)
        self.logger = logging.getLogger(__name__)
        self.setup_logging(level='DEBUG')

        # check mem2gene.txt for HGNC symbol resolution

        self.reset()

    def reset(self):
        """ Reset state for a next genelist to annotate """

        self.delimiter = '|' # join elements of a field
        self.current_line = {}
        self.current_omim_morbid = ''
        self.line_nr = 0
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
        """print only if the verbose switch has been set.
        Warn is used when a value in the genelist will be overwritten.

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
                                    'omim_morbid': self.current_omim_morbid})
            elif self.report_empty or (key and
                key in self.current_line and self.current_line[key]):
                self.logger.warning(line, extra={'line_nr': self.line_nr,
                                    'omim_morbid': self.current_omim_morbid})

    def info(self, line):
        """print only if the verbose switch has been set

        Args:
                line (str): line to print to STDOUT

        Returns:
                pass
        """
        if self.print_info:
            self.logger.info(line, extra={'line_nr': self.line_nr, 'omim_morbid': self.current_omim_morbid})

    def error(self, line):
        """print only if the verbose switch has been set
        Error is used when a mandatory value in the genelist cannot be retrieved.
        e.g. Ensembl_gene_id cannot be filled in.

        Args:
                line (str): line to print to STDOUT

        Returns:
                pass
        """
        if self.print_error:
            line = '\033[31m [ERROR]\033[93m ' + line # add some color
            self.logger.error(line, extra={'line_nr': self.line_nr, 'omim_morbid': self.current_omim_morbid})

    def get_context(self, data):
        """Increments the global line_nr for each passing line.

        Args:
            lines (list of dicts): each dict contains a dict with header as keys
        Yields: line

        """
        for line in data:
            self.current_omim_morbid = line['omim_morbid']
            self.line_nr += 1
            self.current_line = line

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
           line will take precedence over client. Changes will be reported.

           Following will not be reported:
           - Gene_start, Gene_stop changes
           - Updating a missing value
           - Updating HGNC_symbol when symbol is present in client's HGNC_symbols

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
                caller = sys._getframe(1).f_code.co_name
                if str(line[key]) != str(value):
                    # don't report HGNC mismatches if multiple given
                    #if key == 'HGNC_symbol' and line[key] in client['HGNC_symbols']:
                    #    continue
                    self.warn("[{}] {}: line '{}' differs from client '{}'".\
                              format(caller, key, line[key], value), key)

        merged = client.copy()
        merged.update(line)
        return merged

    def remove_comments(self, data):
        for line in data:
            if not line[0].startswith('#'):
                yield line

    def get_chromosome(self, data):
        for line in data:
            cyto_location = there(line, 'Cyto Location')
            new_line = {}
            if cyto_location:
                if any([x for x in cyto_location if x in ('q', 'p')]):
                    chromosome = re.compile('p|q').split(cyto_location)[0]
                    new_line['Chromosome'] = str(chromosome)
                else:
                    chromosome = re.sub(r'^Chr.', '', cyto_location)
                    new_line['Chromosome'] = str(chromosome)

            yield self.merge_line(new_line, line)

    def get_phenotype_number(self, data):
        for line in data:
            new_line = {}
            description = there(line, 'Description')

            if description:
                p = re.compile('.*, (\d+).*')
                m = p.search(description)
                if m:
                    new_line['phenotype_number'] = m.group(1)

            yield self.merge_line(new_line, line)

    def query_omim(self, data):
        """Queries OMIM to fill in the inheritance models

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: with the added HGNC symbol prepended to the HGNC_symbol column.
        """
        omim = OMIM(api_key=self.config['OMIM']['api_key'])
        for line in data:
            omim_morbid = there(line, 'omim_morbid')
            chromosome = there(line, 'Chromosome')
            phenotype_number = there(line, 'phenotype_number')

            if omim_morbid and 'Chromosome' in line:
                entry = omim.gene(mim_number=omim_morbid)
            else:
                func_name = sys._getframe().f_code.co_name
                self.warn('[{}] No entry in OMIM!')
                yield line
                continue

            new_line = {}

            new_line['omim_morbid'] = entry['mim_number']
            new_line['Cyto Location'] = entry['gene_location']
            if entry['gene_location']:
                if any([x for x in entry['gene_location'] if x in ('q', 'p')]):
                    chromosome = re.compile('p|q').split(entry['gene_location'])[0]
                    new_line['Chromosome'] = str(chromosome)

            # extract the inheritance model
            phenotypic_disease_models = omim.\
                parse_phenotypic_disease_models_ext(entry['phenotypes'],
                                                    chromosome=line['Chromosome'],
                                                    phenotype_number=phenotype_number)

            line_phenotypic_disease_models = []
            # if any inheritance models and omim numbers are present, use them!
            for omim_number, models_descriptions in phenotypic_disease_models.items():
                if omim_number is not None:
                    for models_description in models_descriptions:
                        new_line['inheritance_models'] = ','.join(models_description['models'])
                        new_line['omim_phenotype'] = omim_number
                        new_line['description'] = models_description['description']
                        yield self.merge_line(new_line, line)

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
        template = "#%(line_nr)s [%(omim_morbid)s] %(message)s"
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

    def annotate(self, lines, warn=False, error=False, info=False, report_empty=False):
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

        line = next(parsable_data)
        header = line # get the header
        if header[0].startswith('#'):
            header[0] = header[0].lstrip('#')

        # remove the comments
        comment_data = self.remove_comments(parsable_data)

        dict_data = self.list2dict(header, comment_data)

        # get some context for error messages
        context_data = self.get_context(dict_data)

        # fill in the chromosome
        chromosome_data = self.get_chromosome(context_data)

        # fill in the phenotype number
        phenotype_data = self.get_phenotype_number(chromosome_data)

        # fill in the inheritance models, chromosome
        omim_data = self.query_omim(phenotype_data)

        # fill in missing values with ''
        completed_data = self.fill(omim_data)

        print_data = []
        for line in completed_data:
            print(self.format_line(line))
            print_data.append(line)

        # print the errors and warnings
        if verbose:
            # split the lines for easier unit testing
            for line in self.get_log_messages().split('\n'):
                yield line

        # print the gene list
        yield self.get_header()
        for line in print_data:
            yield self.format_line(line)
