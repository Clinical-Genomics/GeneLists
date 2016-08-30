""" Provide all annotation functionality in one place."""
# encoding: utf-8

# TODO move all services init to the constructor

from __future__ import print_function
import pymysql
import re
import os
import yaml
from io import StringIO
from collections import Counter

from ..services.omim import OMIM
from ..services.ensembl import Ensembl
from ..services.genenames import Genenames
from ..services.uniprot import Uniprot
from ..services.mim2gene import Mim2gene

class Genelist(object):

    """Provide all annotation functionality in one class. """

    def __init__(self, config):
        self.gl_header = ['Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_symbol', 'Protein_name',
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
        self.prefix_header = ['Protein_name',
                              'Symptoms', 'Biochemistry', 'Imaging',
                              'Trivial_name_short',
                              'Phenotypic_disease_model', 'OMIM_morbid', 'UniProt_id',
                              'Ensemble_transcript_ID', 'Reduced_penetrance',
                              'Disease_associated_transcript',
                              'Ensembl_transcript_to_refseq_transcript', 'Gene_description',
                              'HGNC_RefSeq_NM', 'Uniprot_protein_name']

        self.config = yaml.load(config)

        # EnsEMBL connection
        self.conn = pymysql.connect(
            host=self.config['ensembl']['host'],
            port=self.config['ensembl']['port'],
            user=self.config['ensembl']['user'],
            db=self.config['ensembl']['db'])

        self.reset()

    def reset(self):
        """ Reset state for a next genelist to annotate """

        self.delimiter = '|' # join elements of a field
        self.current_hgnc_id = ''
        self.line_nr = 0
        self.verbose = False
        self.errors_only = False
        self.mim2gene = False
        self.contigs = set()
        self.error_buffer = StringIO()

    def print_header(self, header=None):
        """
        Prints the contigs meta data headers.
        Prints the gl_header.

        Args:
                header (list, optional): a list of strings
        Note:
                Prints to STDOUT
        """
        if not header:
            header = self.gl_header
        if not self.errors_only:
            print(self.get_header(header=header))

    def get_header(self, header=None):
        """
        Returns the gl_header.

        Args:
                header (list, optional): a list of strings
        Return (str): the header
        """
        if not header:
            header = self.gl_header

        return '#' + "\t".join(self.gl_header)

    def get_contigs(self):
        """ Returns the contigs meta data headers.

        Yields (list): one item per header
        """

        for contig in sorted(self.contigs, key=lambda item:
                             (int(item) if item.isdigit() else float('inf'), item)):
            yield '##contig=<ID={}'.format(contig)

    def format_line(self, line):
        """Formats a line based on the order of the gl_headers

        Args:
            line (dict): dict with values for one line of a gene list.
                         All values of gl_headers should be present as keys.

        Returns (str): a tab-delim line
        """
        ordered_line = list()
        for column_name in self.gl_header:
            ordered_line.append(str(line[column_name]))
        return "\t".join(ordered_line)

    def print_line(self, line):
        """Prints a line based on the order of the gl_headers

        Args:
                line (dict): dict with values for one line of a gene list.
                             All values of gl_headers should be present as keys

        Returns:
                pass

        Note:
                Will print the STDOUT
        """
        if not self.errors_only:
            print(self.format_line(line))

    def warn(self, line):
        """print only if the verbose switch has been set

        Args:
                line (str): line to print to STDOUT

        Returns:
                pass
        """
        if self.verbose:
            warning = '>>> #{} [{}] {}'.format(self.line_nr, self.current_hgnc_id, line)

            print('\033[93m' + warning + '\033[0m')
            self.error_buffer.write(warning + '\n')

    def get_context(self, data):
        """Increments the global line_nr for each passing line.
        Gets the HGNC identifier so we can make more sensical warning messages.

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys
        Yields: line

        """
        for line in data:
            self.line_nr += 1
            self.current_hgnc_id = line['HGNC_symbol']

            yield line

    def query(self, data, try_hgnc_again=False):
        """Queries EnsEMBL. Parameters are HGNC_symbol and/or Ensembl_gene_id, whatever is
        available. Data from EnsEMBLdb will overwrite the client data.
        A(n) identifier(s) should yield one result from EnsEMBLdb.
        It will be reported if a(n) identifier(s) don't yield any or multiple results.

        Args:
                data (list of dicts): representing the lines and columns in a gene list.
                    The keys of the dicts must match the column names of the EnsEMBLdb query.
                try_hgnc_again (bool): when providing multiple HGNC ids, try until you find a
                match on EnsEMBLdb. Only one HGNC ID will be used in final result.

        Yields:
                dict: a row with data from ensEMBLdb filled in.
        """

        cur = self.conn.cursor(pymysql.cursors.DictCursor)
        base_query = """
        SELECT g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop,
        x.display_label AS HGNC_symbol, g.stable_id AS Ensembl_gene_id,
        seq_region.name AS Chromosome
        FROM gene g JOIN xref x ON x.xref_id = g.display_xref_id
        join seq_region USING (seq_region_id)
        """
        keys_conds = {
            'HGNC_symbol': 'x.display_label',
            'Ensembl_gene_id': 'g.stable_id',
            'Chromosome': 'seq_region.name'
        }
        # these columns will be put into the condition statement if they have a value
        keys = ['HGNC_symbol', 'Ensembl_gene_id', 'Chromosome']
        for line in data:
            hgnc_symbols = [symbol for symbol in line['HGNC_symbol'].split(',') if len(symbol) > 0]
            hgnc_symbol_i = 1
            for hgnc_symbol in hgnc_symbols:
                line['HGNC_symbol'] = hgnc_symbol # actually replace the entry
                conds = ["%s = %%s" % keys_conds[key] for key in keys
                         if key in line and line[key] != None and len(line[key]) > 0]
                cond_values = [line[key] for key in keys
                               if key in line and line[key] != None and len(line[key]) > 0]
                # check on length of the region name to exclude scaffolds and patches
                query = "%s where length(seq_region.name) < 3 and %s" % \
                        (base_query, " and ".join(conds))
                cur.execute(query, cond_values)
                rs = cur.fetchall() # result set
                if len(rs) == 0:
                    if hgnc_symbol_i == len(hgnc_symbols):
                        not_found_id = hgnc_symbol if len(hgnc_symbols) == 1 else hgnc_symbols
                        self.warn("Not found: %s %s" % (not_found_id, cond_values))
                        yield line # evenif we don't find an entry for it on ensEMBL
                    if not try_hgnc_again:
                        break
                elif len(rs) > 1:
                    if hgnc_symbol_i > 1:
                        self.warn("Took %s/%s" % (hgnc_symbol, hgnc_symbols))
                    # we couldn't resolve this with genenames.org
                    if 'Chromosome' in line:
                        self.warn("Multiple entries: %s, chromosome: %s => " %
                               (hgnc_symbol, line['Chromosome']))
                    else:
                        self.warn("Multiple entries: %s => " % (hgnc_symbol))
                    self.warn("Adding: %s" % ', '.join((entry['Ensembl_gene_id'] for entry in rs)))
                    for entry in rs:
                        yield self.merge_line(entry, line)
                    break
                else:
                    if len(hgnc_symbols) > 1:
                        self.warn("Took %s/%s" % (hgnc_symbol, hgnc_symbols))
                    for entry in rs:
                        yield self.merge_line(entry, line)
                    break
                hgnc_symbol_i += 1

    def remove_hgnc_prefix(self, line):
        """ Removes the prefixed HGNC symbol from all fields

        Args
            line (dict): representing the line and columns in a gene list.

        Yields (dict):
            a line with the prefixed HGNC identifiers removed
        """

        for column in self.gl_header:
            if column in line:
                re.sub(r'^.*?:', '', str(line[column]))
        return line

    def prepend_hgnc(self, data):
        for line in data:
            for column in self.prefix_header:
                if column in line and line[column]:
                    line[column] = '%s:%s' % (line['HGNC_symbol'], line[column])

            yield line

    def query_transcripts(self, data):
        """Queries EnsEMBL for all transcripts.

        Args
            data (list of dicts): representing the lines and columns in a gene list.
                The keys of the dicts must match the column names of the EnsEMBLdb query.

        Yields (dict):
            a row with transcript data from ensEMBLdb filled in.
        """

        with Ensembl() as ensembldb:
            for line in data:
                if 'Ensembl_gene_id' in line and line['Ensembl_gene_id']:
                    transcripts = ensembldb.query_transcripts(line['Ensembl_gene_id'])
                    if transcripts is not None:
                        line = self.merge_line(transcripts, line)
                yield line

    def merge_line(self, line, client):
        """Will merge line with client.
           ens will take precedence over client. Changes will be reported.

        Args:
                ens (dict): dict with values for one line of a gene list
                client (dict): with values for one line of a gene list

        Yields:
                dict: merged ens and client dict
        """
        for key, value in client.items():
            # don't report start/stop mismatches
            if key in ('Gene_start', 'Gene_stop'):
                continue
            elif key in line:
                if str(line[key]) != str(value):
                    self.warn("%s: line '%s' differs from client '%s'" %
                           (key, line[key], value))

        merged = client.copy()
        merged.update(line)
        return merged

    def fill(self, data):
        """Removes #NA's

        Args:
                data (list of dicts): representing the lines and columns in a gene list

        Yields:
                dict: with all missing columns filled in with ''
        """
        defaults = dict((column_name, '') for column_name in self.gl_header)
        for line in data:
            d = defaults.copy()
            d.update(line)
            yield d

    def munge(self, data):
        """Make sure the data we got from EnsEMBL is good enough for the gene lists
           swap start and stop gene coordinates if start if bigger than stop

        Args:
                data (list of dicts): representing the lines and columns in a gene list

        Yields:
                dict: a row with some munged data
        """
        # swap coordinates if start > stop
        for line in data:
            if line['Gene_start'] > line['Gene_stop']:
                line['Gene_stop'], line['Gene_start'] = line['Gene_start'], line['Gene_stop']

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

    def zero2one(self, data):
        """Fix 0-based coordinates to 1-based coordinates

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: the fixed coordinates data dict
        """
        for line in data:
            for key in ('Gene_start', 'Gene_stop'):
                line[key] = int(line[key]) + 1
            yield line

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
                    self.warn("Add mim2gene HGNC symbol %s" % hgnc_symbol)
                    hgnc_symbols.insert(0, hgnc_symbol)
                if ensembl_gene_id != False and 'Ensembl_gene_id' in line.keys() \
                   and line['Ensembl_gene_id'] != ensembl_gene_id:
                    self.warn("morbidmap '{}' differs from local '{}'".\
                           format(line['Ensembl_gene_id'], ensembl_gene_id))
            line['HGNC_symbol'] = ','.join(hgnc_symbols)
            yield line

    def add_genome_build(self, data, genome_build):
        """Fills in the genome release version in the Clinical_db_genome_build column

        Args:
                data (list of dicts): Inner dict represents a row in a gene list
                genome_build (str): The genome build, e.g. GRCh37

        Yields:
                dict: with the filled in genome build in the Clinical_db_genome_build column
        """
        for line in data:
            line['Clinical_db_genome_build'] = genome_build
            yield line

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

    def remove_non_genes(self, data):
        """Based on mim2gene.txt, you can remove all non genes. The global type_of dict provides the type of the hgnc symbol.

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: with all-non genes removed
        """
        for line in data:
            if 'HGNC_symbol' not in line: # can't use mim2gene if there is no HGNC symbol
                yield line
            else:
                yielded = False
                hgnc_symbols = line['HGNC_symbol'].split(',')
                for hgnc_symbol in hgnc_symbols:
                    if self.mim2gene.is_gene(hgnc_symbol):
                        yield line
                        yielded = True
                        break
                if not yielded:
                    self.warn('Removed: {%s}' % ', '.join(["'%s':'%s'" % (k, line[k]) for k in sorted(line)]))

    def put_official_hgnc_symbol(self, data):
        """Resolve the official HGNC symbol from OMIM (mim2gene) and replace line['HGNC_symbol']

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: now with the official HGNC_symbol
        """
        for line in data:
            if 'OMIM_morbid' in line:
                hgnc_symbol = self.mim2gene.resolve_gene(line['OMIM_morbid'])
                if hgnc_symbol != False and line['HGNC_symbol'] != hgnc_symbol:
                    self.warn('Took official symbol {} over {}'.\
                           format(hgnc_symbol, line['HGNC_symbol']))
                    line['HGNC_symbol'] = hgnc_symbol
            yield line

    def add_official_hgnc_symbol(self, data):
        """Add the official HGNC symbol fetched from genenames.org to the field
        Official_HGNC_symbol. Also prepend the official HGNC symbol to the HGNC_symbol field.

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: now with the official HGNC_symbol prepended
        """
        genenames = Genenames()
        for line in data:
            OMIM_morbid=None
            if 'OMIM_morbid' in line and line['OMIM_morbid']:
                OMIM_morbid = line['OMIM_morbid']
            #HGNC_symbol = line['HGNC_symbol'].split(',')[-1] # take the last symbol
            official_symbols = []
            HGNC_symbols = line['HGNC_symbol'].split(',') # take the last symbol
            for HGNC_symbol in HGNC_symbols:
                official_symbol = genenames.official(HGNC_symbol, OMIM_morbid)
                if official_symbol:
                    if official_symbol not in HGNC_symbols:
                        self.warn('Add official HGNC symbol %s' % official_symbol)
                        line['HGNC_symbol'] = ','.join((official_symbol, line['HGNC_symbol']))
                    official_symbols.append(official_symbol)
                else:
                    official_symbols.append(HGNC_symbol)

            self.warn('Official symbols %s' % official_symbols)
            # ok, now we have several 'official' symbols. Take the most present one.
            # Create a Counter object, which will 'count' the items in the list and present a
            # list with (value, count) tuples.
            # Sort that list ascendingly.
            sorted_count = sorted(Counter(official_symbols).items(), key=lambda x: x[1])
            # take the last element (the most present one) and of that the value of the tuple
            line['Official_HGNC_symbol'] = sorted_count[-1][0]
            self.warn('Took %s as official symbol' % line['Official_HGNC_symbol'])

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
            uniprot_ids = genenames.uniprot(line['Official_HGNC_symbol'])
            uniprot_ids = uniprot_ids if uniprot_ids != None else ''
            uniprot_ids_joined = self.delimiter.join(uniprot_ids)
            if len(uniprot_ids) > 1:
                self.warn('Multiple UniProt IDs: ' + uniprot_ids_joined)
            if 'UniProt_id' in line and line['UniProt_id']:
                self.warn('Replaced Uniprot ID %s with %s' % (line['UniProt_id'], uniprot_ids_joined))

            for uniprot_id in uniprot_ids:
                uniprot_description = uniprot.fetch_description(uniprot_id)
                if 'Uniprot_protein_name' in line and line['Uniprot_protein_name']:
                    self.warn('Replaced Uniprot ID %s with %s' %
                           (line['Uniprot_protein_name'], uniprot_description))

            line['Uniprot_protein_name'] = uniprot_description
            line['UniProt_id'] = uniprot_ids_joined

            yield line

    def add_refseq(self, data):
        """ Add the RefSeq ID based on the official HGNC symbol.

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: now with the RefSeq information.
        """
        genenames = Genenames()
        for line in data:
            refseq = genenames.refseq(line['Official_HGNC_symbol'])
            refseq = self.delimiter.join(refseq) if refseq != None else ''
            if 'HGNC_RefSeq_NM' in line and line['HGNC_RefSeq_NM']:
                self.warn('Replaced HGNC_RefSeq_NM %s with %s' % (line['HGNC_RefSeq_NM'], refseq))
            line['HGNC_RefSeq_NM'] = refseq

            yield line

    def query_omim(self, data):
        """Queries OMIM to fill in the inheritance models

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: with the added HGNC symbol prepended to the HGNC_symbol column.
        """
        omim = OMIM(api_key=self.config['OMIM']['api_key'])
        for line in data:
            if 'OMIM_morbid' in line and line['OMIM_morbid'] and 'Chromosome' in line:
                entry = omim.gene(mim_number=line['OMIM_morbid'])
            elif 'HGNC_symbol' in line and 'Chromosome' in line:
                entry = omim.gene(hgnc_symbol=line['Official_HGNC_symbol'])
            else:
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

            if len(line_phenotypic_disease_models) > 0:
                line['Phenotypic_disease_model'] = '|'.join(line_phenotypic_disease_models)
            else:
                line['Phenotypic_disease_model'] = ''
            # add OMIM morbid
            if entry['mim_number'] is not None:
                if 'OMIM_morbid' in line \
                and len(line['OMIM_morbid']) > 0 \
                and str(line['OMIM_morbid']) != str(entry['mim_number']):
                    self.warn('%s %s > %s client OMIM number differs from OMIM query' % \
                           (line['HGNC_symbol'], line['OMIM_morbid'], entry['mim_number']))
                line['OMIM_morbid'] = entry['mim_number']

            # add Gene_locus
            if entry['gene_location'] is not None:
                if 'Gene_locus' in line and len(line['Gene_locus']) > 0 and \
                   line['Gene_locus'] != entry['gene_location']:
                    self.warn('%s > %s client Gene locus differs from OMIM query' % \
                           (line['Gene_locus'], entry['gene_location']))
                line['Gene_locus'] = entry['gene_location']
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

    def annotate(self, lines, verbose=False, errors=False, download_mim2gene=False, mim2gene=False, zero=False):
        """ Annotate a gene list """

        self.reset()

        # make sure we print if we are asked to
        if verbose:
            self.verbose = True

        # show only the EnsEMBLdb conflicts - so verbose, but supress printing of the gene list
        if errors:
            self.verbose = True
            self.errors_only = True

        raw_data = (line.strip() for line in lines) # sluuuurp
        parsable_data = (line.split("\t") for line in raw_data)

        # check mem2gene.txt for HGNC symbol resolution
        if mim2gene or download_mim2gene:
            mim2gene_filename = os.path.join(os.path.dirname(__file__), 'mim2gene.txt')
            if download_mim2gene:
                self.warn('Downloading {} ... '.format(mim2gene_filename))
                self.mim2gene = Mim2gene(filename=mim2gene_filename, download=True)
            else:
                self.mim2gene = Mim2gene()

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
        if zero:
            # fix 0-based coordinates to be 1-based
            fixed_data = self.zero2one(clean_data)
        else:
            fixed_data = clean_data

        # add the mim2gene alias to HGNC_symbol
        # remove non-genes based on mim2gene.txt
        if self.mim2gene:
            aliased_data = self.add_mim2gene_alias(fixed_data)
            reduced_data = self.remove_non_genes(aliased_data)
        else:
            reduced_data = fixed_data
        aliased_data = self.add_official_hgnc_symbol(reduced_data)

        # add uniprot
        uniprot_data = self.add_uniprot(aliased_data)

        # add refseq
        refseq_data = self.add_refseq(uniprot_data)

        # fill in missing blanks
        ensembld_data = self.query(refseq_data, try_hgnc_again=True)

        # put the official HGNC symbol
        if mim2gene:
            hgnc_official_data = self.put_official_hgnc_symbol(ensembld_data)
        else:
            hgnc_official_data = ensembld_data

        # aggregate transcripts
        transcript_data = self.query_transcripts(hgnc_official_data)

        # fill in the inheritance models
        omim_data = self.query_omim(transcript_data)

        # do some replacements
        redpen_data = self.redpen2symbol(omim_data)

        # fill in missing values with ''
        completed_data = self.fill(redpen_data)

        # clean up the data from EnsEMBL a bit
        munged_data = self.munge(completed_data)

        # at last, clean up the output
        cleaner_data = self.cleanup(munged_data)

        # prepend the HGNC symbol to some fields
        prefixed_data = self.prepend_hgnc(cleaner_data)

        # get all contigs
        final_data = self.gather_contig(prefixed_data)
        print_data = []
        for line in final_data:
            if self.verbose:
                print(self.format_line(line))
            print_data.append(line)

        # print the errors and warnings
        if verbose:
            # split the lines for easier unit testing
            for line in self.error_buffer.getvalue().split('\n'):
                yield line

        # print the gene list
        for comment in comments:
            yield '\t'.join(comment)
        for line in self.get_contigs():
            yield line
        yield self.get_header()
        for line in print_data:
            yield self.format_line(line)
