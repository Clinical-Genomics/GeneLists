""" Provide all annotation functionality in one place."""
# encoding: utf-8

from __future__ import print_function
import sys
import re
import os
import yaml
import logging
from io import StringIO

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
        self.logger = logging.getLogger(__name__)
        self.setup_logging(level='DEBUG')

        self.reset()

        # check mem2gene.txt for HGNC symbol resolution
        self.mim2gene = self.init_mim2gene(download_mim2gene)
        self.ensembldb = Ensembl(
            host=self.config['ensembl']['host'],
            port=self.config['ensembl']['port'],
            user=self.config['ensembl']['user'],
            db=self.config['ensembl']['db']
        )

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
        self.remove_non_genes = False

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
        Error is used when a mandatory value in the genelist cannot be retrieved.
        e.g. Ensembl_gene_id cannot be filled in.

        Args:
                line (str): line to print to STDOUT

        Returns:
                pass
        """
        if self.print_error:
            line = '\033[31m [ERROR]\033[93m ' + line # add some color
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
                    if line[key] == 'NA':
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

    def pick_hgnc_symbol(self, data):
        """ If multiple HGNC symbols, pick first one as main symbol.
        Rest is stored in HGNC_symbols.
        """
        for line in data:
            hgnc_symbol = there(line, 'HGNC_symbol')
            hgnc_symbols = hgnc_symbol.split(',')
            line['HGNC_symbols'] = hgnc_symbols # all symbols as a list
            line['HGNC_symbol'] = hgnc_symbols[0] # the proper HGNC smbol
            line['HGNC_symbol_start'] = hgnc_symbols[0] # in case the above changes we still have the oriinal one

            yield line


    def init_mim2gene(self, download_mim2gene):
        mim2gene_filename = os.path.join(os.path.dirname(__file__), 'mim2gene.txt')
        if download_mim2gene:
            self.info('Downloading {} ... '.format(mim2gene_filename))
            return Mim2gene(filename=mim2gene_filename, download=True)
        else:
            return Mim2gene(filename=mim2gene_filename)

    def remove_from_mim2gene(self, data):
        """Based on mim2gene.txt, you can remove all non genes.

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: with all-non genes removed
        """
        for line in data:
            if not self.remove_non_genes or self.mim2gene.is_gene(line['OMIM_morbid']):
                yield line
            else:
                self.warn('Removed non gene: {}'.format(line['HGNC_symbol']))

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
            omim_morbid = there(line, 'OMIM_morbid')
            hgnc_symbol = there(line, 'HGNC_symbol')

            #ensembl_gene_id = there(line, 'Ensembl_gene_id')
            func_name = sys._getframe().f_code.co_name

            if omim_morbid:
                hgnc_symbol = self.mim2gene.get_hgnc(omim_morbid)
                if not hgnc_symbol:
                    self.error('[{}] HGNC_symbol NOT FOUND!'.format(func_name))
                    yield line
                else:
                    yield self.merge_line(
                        {
                            'HGNC_symbol': hgnc_symbol,
                            'Ensembl_gene_id': self.mim2gene.get_ensembl(omim_morbid)
                        },
                        line,
                    )
                continue

            if hgnc_symbol:
                omim_morbid = self.mim2gene.get_omim(hgnc_symbol)
                if not omim_morbid:
                    self.error('[{}] OMIM_morbid not found!'.format(func_name))
                    yield line
                else:
                    yield self.merge_line(
                        {
                            'OMIM_morbid': omim_morbid,
                            'Ensembl_gene_id': self.mim2gene.get_ensembl(hgnc_symbol)
                        },
                        line,
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
            dict: with the Gene_start, Gene_stop, Chromosome and HGNC_symbol filled in.

        TODO: change the order in which E! is queried when we would trust the E! gene ids:
        - E! + OMIM
        - E! + HGNC
        - OMIM
            - OMIM + HGNC on multiple
        - HGNC
        All queries include the chromosome

        """
        func_name = sys._getframe().f_code.co_name
        for line in data:
            omim_morbid = there(line, 'OMIM_morbid')
            hgnc_symbol = there(line, 'HGNC_symbol')
            chromosome = there(line, 'Chromosome')
            ensembl_gene_id = there(line, 'Ensembl_gene_id')

            ensembl_lines = []

            if ensembl_gene_id and omim_morbid:
                ensembl_lines = self.ensembldb.query(ensembl_gene_id=ensembl_gene_id, omim_morbid=omim_morbid, chromosome=chromosome)
                if ensembl_lines:
                    self.info('[{}] Found E! with {} {} {}'.format(func_name, ensembl_gene_id, omim_morbid, chromosome))

            if not ensembl_lines:
                if ensembl_gene_id and hgnc_symbol:
                    ensembl_lines = self.ensembldb.query(ensembl_gene_id=ensembl_gene_id, hgnc_symbol=hgnc_symbol, chromosome=chromosome)
                    if ensembl_lines:
                        self.info('[{}] Found E! with {} {} {}'.format(func_name, ensembl_gene_id, hgnc_symbol, chromosome))

            if not ensembl_lines and omim_morbid:
                ensembl_lines = self.ensembldb.query(omim_morbid=omim_morbid, chromosome=chromosome)
                if ensembl_lines:
                    self.info('[{}] Found E! with {}'.format(func_name, omim_morbid, chromosome))

                # multiple hits? WTF. Check with the hgnc symbol and omim morbid
                if len(ensembl_lines) > 1 and hgnc_symbol:
                    ensembl_lines = self.ensembldb.query(hgnc_symbol=hgnc_symbol, omim_morbid=omim_morbid, chromosome=chromosome)
                    self.info('[{}] Found E! with {} {} {}'.format(func_name, omim_morbid, hgnc_symbol, chromosome))

            if not ensembl_lines and hgnc_symbol:
                # then with the HGNC symbol only
                ensembl_lines = self.ensembldb.query(hgnc_symbol=hgnc_symbol, chromosome=chromosome)
                if ensembl_lines:
                    self.info('[{}] Found E! with {}'.format(func_name, hgnc_symbol))

            if ensembl_lines:
                if len(ensembl_lines) > 1:
                    e_ids = [entry['Ensembl_gene_id'] for entry in ensembl_lines]
                    #if ensembl_gene_id in e_ids:
                    self.info('[{}] Multiple E! entries: {}.'.format(func_name, e_ids))
                for ensembl_line in ensembl_lines:
                    yield self.merge_line(ensembl_line, line)
            else:
                self.error('[{}] {}: No E! entries!'.format(func_name, omim_morbid))
                yield line

    def query_transcripts(self, data):
        """Queries EnsEMBL for all transcripts.

        Args
            data (list of dicts): representing the lines and columns in a gene list.
                The keys of the dicts must match the column names of the EnsEMBLdb query.

        Yields (dict):
            a row with transcript data from ensEMBLdb filled in.
        """

        func_name = sys._getframe().f_code.co_name
        for line in data:
            omim_morbid = there(line, 'OMIM_morbid')
            ensembl_gene_id = there(line, 'Ensembl_gene_id')

            transcripts = []
            if omim_morbid and ensembl_gene_id:
                transcripts = self.ensembldb.query_transcripts_omim(omim_morbid=omim_morbid, ensembl_gene_id=ensembl_gene_id)
                if transcripts:
                    self.info('[{}] Found E! transcripts with {} {}'.format(func_name, omim_morbid, ensembl_gene_id))

            if not transcripts and ensembl_gene_id:
                transcripts = self.ensembldb.query_transcripts_omim(ensembl_gene_id=ensembl_gene_id)
                if transcripts:
                    self.info('[{}] Found E! transcripts with {}'.format(func_name, ensembl_gene_id))

            if transcripts:
                line = self.merge_line(transcripts, line)
            else:
                self.warn('[{}] No transcripts on E!'.format(func_name))

            yield line

    def add_uniprot(self, data):
        """ Add the UniProt ID and UniProt protein name based on the official HGNC symbol.

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: now with the UniProt information.
        """
        genenames = Genenames()
        uniprot = Uniprot()
        for line in data:
            uniprot_ids = genenames.uniprot(line['HGNC_symbol'])
            uniprot_ids = uniprot_ids if uniprot_ids != None else ''
            uniprot_ids_joined = self.delimiter.join(uniprot_ids)
            uniprot_description = ''
            if len(uniprot_ids) > 1:
                self.info('Multiple UniProt IDs: ' + uniprot_ids_joined)
            for uniprot_id in uniprot_ids:
                uniprot_description = uniprot.fetch_description(uniprot_id)

            yield self.merge_line(
                {
                    'Uniprot_protein_name': uniprot_description,
                    'UniProt_id': uniprot_ids_joined
                },
                line
            )

    def add_refseq(self, data):
        """ Add the RefSeq ID based on the official HGNC symbol.

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: now with the RefSeq information.
        """
        genenames = Genenames()
        for line in data:
            refseq = genenames.refseq(line['HGNC_symbol'])
            refseq = self.delimiter.join(refseq) if refseq != None else ''
            yield self.merge_line({'HGNC_RefSeq_NM': refseq}, line)

    def query_omim(self, data):
        """Queries OMIM to fill in the inheritance models

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: with the added HGNC symbol prepended to the HGNC_symbol column.
        """
        omim = OMIM(api_key=self.config['OMIM']['api_key'])
        for line in data:
            omim_morbid = there(line, 'OMIM_morbid')
            if omim_morbid and 'Chromosome' in line:
                entry = omim.gene(mim_number=omim_morbid)
            elif 'HGNC_symbol' in line and 'Chromosome' in line:
                entry = omim.gene(hgnc_symbol=line['HGNC_symbol'])
            else:
                func_name = sys._getframe().f_code.co_name
                self.warn('[{}] No entry in OMIM!')
                yield line
                continue

            phenotypic_disease_models = omim.\
                parse_phenotypic_disease_models(entry['phenotypes'], line['Chromosome'])

            # extract the inheritance model
            line_phenotypic_disease_models = []
            # if any inheritance models and omim numbers are present, use them!
            for omim_number, inheritance_models in phenotypic_disease_models.items():
                if omim_number is not None:
                    inheritance_models_str = ''
                    if inheritance_models is not None:
                        inheritance_models_str = '>' + '/'.join(inheritance_models)
                    line_phenotypic_disease_models.append('%s%s' % ( \
                        omim_number,
                        inheritance_models_str))

            new_line = {}

            if len(line_phenotypic_disease_models) > 0:
                new_line['Phenotypic_disease_model'] = '|'.join(line_phenotypic_disease_models)
            else:
                new_line['Phenotypic_disease_model'] = ''

            new_line['OMIM_morbid'] = entry['mim_number']
            new_line['Gene_locus'] = entry['gene_location']
            if entry['gene_location']:
                if any([x for x in entry['gene_location'] if x in ('q', 'p')]):
                    chromosome = re.compile('p|q').split(entry['gene_location'])[0]
                    new_line['Chromosome'] = str(chromosome)

            yield self.merge_line(new_line, line)

    def redpen2symbol(self, data):
        """If reduced penetrance is set, replace it with the HGNC symbol

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: with the replaced red pen to HGNC symbol
        """
        for line in data:
            if 'Reduced_penetrance' in line and line['Reduced_penetrance'].lower() == 'yes':
                line['Reduced_penetrance'] = line['HGNC_symbol']
            yield line

    def fill_alias(self, data):
        """If the starting HGNC ends up beign different than the endin HGNC symbol,
        add it to the alias column. Chck if the gene symbol isn't already present in the
        alias list.
        """
        for line in data:
            hgnc_symbol = there(line, 'HGNC_symbol')
            hgnc_symbol_start = there(line, 'HGNC_symbol_start')
            if hgnc_symbol != hgnc_symbol_start:
                alias = there(line, 'Alias')
                aliases = alias.split(',') if alias else []
                if hgnc_symbol_start not in aliases:
                    aliases.append(hgnc_symbol_start)
                    self.info('Adding {} to Alias'.format(hgnc_symbol_start))
                    line['Alias'] = ','.join(aliases)
            yield line

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

        if remove_non_genes:
            self.remove_non_genes = True

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

        # remove none genes
        reduced_data = self.remove_from_mim2gene(clean_data)

        # pick one HGNC symbol
        hgnc_data = self.pick_hgnc_symbol(reduced_data)

        # Get OMIM morbid number
        # Get E!
        # all from mim2gene, the magical file
        mim2gene_data = self.fill_from_mim2gene(hgnc_data)

        # fill in the inheritance models, chromosome
        omim_data = self.query_omim(mim2gene_data)

        # fill in info from ensembl
        ensembl_data = self.fill_from_ensembl(omim_data)

        # aggregate transcripts
        transcript_data = self.query_transcripts(ensembl_data)

        ## add uniprot
        uniprot_data = self.add_uniprot(transcript_data)

        ## add refseq
        refseq_data = self.add_refseq(uniprot_data)

        ## do some replacements
        redpen_data = self.redpen2symbol(refseq_data)

        # fill in missing values with ''
        completed_data = self.fill(redpen_data)

        # add the alias if any
        aliased_data = self.fill_alias(completed_data)

        # at last, clean up the output
        cleaner_data = self.cleanup(aliased_data)

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
