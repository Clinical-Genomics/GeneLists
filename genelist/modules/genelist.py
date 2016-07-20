""" Provide all annotation functionality in one place."""
# encoding: utf-8

from __future__ import print_function
import pymysql
import argparse
import re
import os
from urllib.request import urlretrieve, Request, urlopen
from collections import Counter

from ..services.omim import OMIM
from ..services.ensembl import Ensembl
from ..services.genenames import Genenames

class Genelist(object):

    """Provide all annotation functionality in one class. """

    def __init__(self):
        self.gl_header = ['Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_symbol', 'Protein_name',
                          'Symptoms', 'Biochemistry', 'Imaging', 'Disease_trivial_name',
                          'Trivial_name_short',
                          'Phenotypic_disease_model', 'OMIM_morbid', 'Gene_locus', 'UniProt_id',
                          'Ensembl_gene_id', 'Ensemble_transcript_ID', 'Reduced_penetrance',
                          'Clinical_db_gene_annotation', 'Disease_associated_transcript',
                          'Ensembl_transcript_to_refseq_transcript', 'Gene_description',
                          'Genetic_disease_model', 'HGNC_RefSeq_NM', 'Uniprot_protein_name',
                          'Database_entry_version', 'Curator', 'Alias', 'Group_or_Pathway']

        self.delimiter = '|' # join element of a field

        # EnsEMBL connection
        self.conn = None
        self.verbose = False # to print or not to print
        # print only errors. Does not print the gene list. Needs verbose=True to work.
        self.errors_only = False
        self.mim2gene = False # resolve HGNC symbol with mim2gene.txt
        self.contigs = set() # a list of al contigs in the list
        self.outfile = None # where to write to

        self.symbol_of = {} # omim_id: official hgnc_symbol
        self.type_of = {} # hgnc_symbol: type
        self.ensembl_gene_id_of = {} # hgnc_symbol: EnsEMBL_gene_id

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

    def p(self, line):
        """print only if the verbose switch has been set

        Args:
                line (str): line to print to STDOUT

        Returns:
                pass
        """
        if self.verbose:
            print('\033[93m', '>>> ', line, '\033[0m')
            self.outfile.write('\033[93m>>> {} \033[0m\n'.format(line))

    def resolve_ensembl_id(self, hgnc_id):
        """Query genenames.org for the EnsEMBL gene id based on the HGNC symbol.

        Args:
                hgnc_id (str): the HGNC symbol

        Returns (str): The ensEMBL gene id

        """
        import json
        response = urlopen(Request("http://rest.genenames.org/fetch/symbol/%s" % hgnc_id,
                                   None, {'Accept':'application/json'}))
        data = response.read().decode('UTF-8')
        data = json.loads(data)
        try:
            return data['response']['docs'][0]['ensembl_gene_id']
        except KeyError:
            return False

    # TODO: djees, put this in a separate package so we don't have to rely on a global var
    def cache_mim2gene(self,
                       mim2gene_file=os.path.dirname(os.path.abspath(__file__))+
                       os.path.sep+
                       'mim2gene.txt'):
        """Read in the mim2gene file and store it as a dict of OMIM id: HGNC_symbol.
        Only gene and gene/phenotype types will be saved.

        Kwargs:
                mim2gene_file (str): the aboslute path to the mim2gene.txt file

        Returns: None
        """

        mim2gene_fh = open(mim2gene_file, 'r')
        lines = (line for line in mim2gene_fh)
        for line in lines:
            if line.startswith('#'):
                continue
            (file_omim_id, omim_type, gene_id, hgnc_symbol, ensembl_gene_id) = line.split("\t")
            if omim_type in ('gene', 'gene/phenotype') and hgnc_symbol:
                self.symbol_of[file_omim_id] = hgnc_symbol
                self.ensembl_gene_id_of[file_omim_id] = ensembl_gene_id
            self.type_of[hgnc_symbol] = omim_type

    def resolve_gene(self, omim_id):
        """Looks up the omim_id in the mim2gene.txt file.
        If found and the omim type is 'gene', return the official HGNC symbol

        Args:
            omim_id (int): the omim id

        Returns:
            on omim morbid match, official HGNC symbol if type of gene or
            gene/phenotype otherwise False
        """

        if omim_id in self.symbol_of:
            return self.symbol_of[omim_id]
        return False

    def resolve_ensembl_gene_id(self, omim_id):
        """Looks up the EnsEMBL gene id in the mim2gene.txt file.
        If found and the omim type is 'gene', return it

        Args:
                omim_id (int): the omim id

        Returns: on omom id match, EnsEMBL gene id if type of gene or gene/phenotype otherwise False
        """

        if omim_id in self.ensembl_gene_id_of and \
           self.ensembl_gene_id_of[omim_id] != None and \
           self.ensembl_gene_id_of[omim_id] != '-':
            return self.ensembl_gene_id_of[omim_id]
        return False

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
            HGNC_symbols = [symbol for symbol in line['HGNC_symbol'].split(',') if len(symbol) > 0]
            HGNC_symbol_i = 1
            for HGNC_symbol in HGNC_symbols:
                line['HGNC_symbol'] = HGNC_symbol # actually replace the entry
                conds = ["%s = %%s" % keys_conds[key] for key in keys
                         if key in line and line[key] != None and len(line[key]) > 0]
                cond_values = [line[key] for key in keys
                               if key in line and line[key] != None and len(line[key]) > 0]
                # check on length of the region name to exclude scaffolds and patches
                query = "%s where length(seq_region.name) < 3 and %s" % (base_query, " and ".join(conds))
                cur.execute(query, cond_values)
                rs = cur.fetchall() # result set
                if len(rs) == 0:
                    if HGNC_symbol_i == len(HGNC_symbols):
                        not_found_id = HGNC_symbol if len(HGNC_symbols) == 1 else HGNC_symbols
                        self.p("Not found: %s %s" % (not_found_id, cond_values))
                        yield line # evenif we don't find an entry for it on ensEMBL
                    if not try_hgnc_again:
                        break
                elif len(rs) > 1:
                    if HGNC_symbol_i > 1:
                        self.p("Took %s/%s" % (HGNC_symbol, HGNC_symbols))
                    # we couldn't resolve this with genenames.org
                    if 'Chromosome' in line:
                        self.p("Multiple entries: %s, chromosome: %s => " %
                               (HGNC_symbol, line['Chromosome']))
                    else:
                        self.p("Multiple entries: %s => " % (HGNC_symbol))
                    self.p("Adding: %s" % ', '.join((entry['Ensembl_gene_id'] for entry in rs )))
                    for entry in rs:
                        yield self.merge_line(entry, line)
                    break
                else:
                    if len(HGNC_symbols) > 1:
                        self.p("Took %s/%s" % (HGNC_symbol, HGNC_symbols))
                    for entry in rs:
                        yield self.merge_line(entry, line)
                    break
                HGNC_symbol_i += 1

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
                if 'Ensembl_gene_id' in line.keys():
                    transcripts = ensembldb.query_transcripts(line['Ensembl_gene_id'])
                    if transcripts is not None:
                        line.update(transcripts)
                        if len(line['Gene_description']) > 0:
                            line['Gene_description'] = \
                                '%s:%s' % (line['HGNC_symbol'], line['Gene_description'])
                yield line

    def get_transcript(self, start, end, ensembl_gene_id=None, hgnc_id=None):
        """Queries EnsEMBL. Parameters are HGNC_symbol and/or Ensembl_gene_id, whatever is
           available. It will return one hit with the ensembl trasncript id.

        Args:
                ensembl_gene_id (str, optional): ensembl_gene_id and/or hgnc_id should be provided.
                hgnc_id (str, optional): ensembl_gene_id and/or hgnc_id should be provided.
                start (int): start coordinate of the possible transcript of this gene.
                end (int): stop coordinate of the possible transcript of this gene.

        Yields:
                dict: with following keys: ensembl_gene_id, hgnc_id, start, end, transcript_id
        """

        cur = self.conn.cursor(pymysql.cursors.DictCursor)
        base_query = "select g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop, x.display_label AS HGNC_symbol, g.stable_id AS Ensembl_gene_id, seq_region.name AS Chromosome, t.seq_region_start AS Transcript_start, t.seq_region_end AS Transcript_stop, t.stable_id AS Transcript_id from gene g join xref x on x.xref_id = g.display_xref_id join seq_region using (seq_region_id) join transcript t using (gene_id)"
        conds = {'t.seq_region_start': start, 't.seq_region_end': end}
        if ensembl_gene_id != None:
            conds.update({'g.stable_id': ensembl_gene_id})
        if hgnc_id != None:
            conds.update({'x.display_label': hgnc_id})
        query = "%s where %s" % (base_query, " and ".\
                join(['%s = %%s' % column for column in conds.keys()]))
        cur.execute(query, [str(value) for value in conds.values()])
        rs = cur.fetchall() # result set
        # O-oh .. for now this still means manual intervention!
        if len(rs) > 1:
            self.p("Getting '%s' ... " % conds.values())
            self.p('Multiple entries found!')
        elif len(rs) == 0:
            return None
        return rs[0]

    def merge_line(self, ens, client):
        """Will merge dict ens (EnsEMBL data) with client (data).
           ens will take precedence over client. Changes will be reported.

        Args:
                ens (dict): dict with values for one line of a gene list
                client (dict): with values for one line of a gene list

        Yields:
                dict: merged ens and client dict
        """
        for key, value in client.items():
            if key in ('Gene_start', 'Gene_stop'):
                pass
            elif key in ens:
                if str(ens[key]) != str(value):
                    self.p("%s > %s: ens '%s' diff from client '%s'" %
                           (ens['Ensembl_gene_id'], key, ens[key], value))
            #        else:
            #                p("%s: ens '%s' eq to client '%s'" % (key, ens[key], value))
            #else:
            #        p("%s not in ens!" % key)

        # Check the Gene_start and Gene_stop for being Transcript coordinates
        if 'Gene_start' in client and 'Gene_stop' in client:
            if str(ens['Gene_start']) != str(client['Gene_start']) \
               or str(ens['Gene_stop']) != str(client['Gene_stop']):
                # get the transcript, compare those coordinates and report
                transcript = self.get_transcript(client['Gene_start'], client['Gene_stop'],
                                                 ens['Ensembl_gene_id'], ens['HGNC_symbol'])
                for key in ('Gene_start', 'Gene_stop'):
                    if str(ens[key]) != str(client[key]):
                        if transcript and len(transcript) > 1:
                            self.p("%s > %s: ens '%s' diff from client '%s',"
                                   "but matches Transcript %s: %s" %
                                   (ens['Ensembl_gene_id'], key, ens[key], client[key],
                                    transcript['Transcript_id'],
                                    transcript[key.replace('Gene', 'Transcript')]))
                        else:
                            self.p("%s > %s: ens '%s' diff from client '%s'" %
                                   (ens['Ensembl_gene_id'], key, ens[key], client[key]))

        merged = client.copy()
        merged.update(ens)
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
                # swap start and stop gene coordinates if start if bigger than stop

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
            HGNC_symbols = line['HGNC_symbol'].split(',')
            if 'OMIM_morbid' in line:
                OMIM_id    = line['OMIM_morbid']
                HGNC_symbol = self.resolve_gene(OMIM_id)
                EnsEMBL_gene_id = self.resolve_ensembl_gene_id(OMIM_id)
                if HGNC_symbol != False and HGNC_symbol not in HGNC_symbols:
                    self.p("Add mim2gene HGNC symbol %s" % HGNC_symbol)
                    HGNC_symbols.insert(0, HGNC_symbol)
                if EnsEMBL_gene_id != False and 'Ensembl_gene_id' in line.keys() \
                   and line['Ensembl_gene_id'] != EnsEMBL_gene_id:
                    self.p("morbidmap '{}' differs from local '{}'".\
                           format(line['Ensembl_gene_id'], EnsEMBL_gene_id))
            line['HGNC_symbol'] = ','.join(HGNC_symbols)
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
                HGNC_symbols = line['HGNC_symbol'].split(',')
                for HGNC_symbol in HGNC_symbols:
                    if HGNC_symbol in self.type_of and \
                       self.type_of[ HGNC_symbol ] in ('gene', 'gene/phenotype'):
                        yield line
                        yielded = True
                        break
                if not yielded:
                    self.p('Removed: %s' % line)

    def put_official_hgnc_symbol(self, data):
        """Resolve the official HGNC symbol from OMIM and replace line['HGNC_symbol']

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: now with the official HGNC_symbol
        """
        for line in data:
            if 'OMIM_morbid' in line:
                HGNC_symbol = self.resolve_gene(line['OMIM_morbid'])
                if HGNC_symbol != False and line['HGNC_symbol'] != HGNC_symbol:
                    self.p('Took official symbol {} over {}'.\
                           format(HGNC_symbol, line['HGNC_symbol']))
                    line['HGNC_symbol'] = HGNC_symbol
            yield line

    def add_official_hgnc_symbol(self, data):
        """Add the official HGNC symbol fetched from genenames.org to the field Official_HGNC_symbol.
        Also prepend the official HGNC symbol to the HGNC_symbol field.

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
                        self.p('Add official HGNC symbol %s' % official_symbol)
                        line['HGNC_symbol'] = ','.join( (official_symbol, line['HGNC_symbol']) )
                    official_symbols.append(official_symbol)
                else:
                    official_symbols.append(HGNC_symbol)

            self.p('Official symbols %s' % official_symbols)
            # ok, now we have several 'official' symbols. Take the most present one.
            # Create a Counter object, which will 'count' the items in the list and present a
            # list with (value, count) tuples.
            # Sort that list ascendingly.
            sorted_count = sorted(Counter(official_symbols).items(), key=lambda x: x[1])
            # take the last element (the most present one) and of that the value of the tuple
            line['Official_HGNC_symbol'] = sorted_count[-1][0]
            self.p('Took %s as official symbol' % line['Official_HGNC_symbol'])

            yield line

    def add_uniprot(self, data):
        genenames = Genenames()
        for line in data:
            uniprot_ids = self.delimiter.join(genenames.uniprot(line['Official_HGNC_symbol']))
            if len(genenames.uniprot(line['Official_HGNC_symbol'])) > 1:
                self.p('Multiple UniProt IDs: ' + ','.join(genenames.uniprot(line['Official_HGNC_symbol'])))
            if 'UniProt_id' in line and line['UniProt_id']:
                self.p('Replaced Uniprot ID %s with %s' % (line['UniProt_id'], uniprot_ids))
            line['UniProt_id'] = uniprot_ids

            yield line

    def add_refseq(self, data):
        genenames = Genenames()
        for line in data:
            refseq = self.delimiter.join(genenames.refseq(line['Official_HGNC_symbol']))
            if len(genenames.refseq(line['Official_HGNC_symbol'])) > 1:
                print('REFSEQ ' + genenames.refseq(line['Official_HGNC_symbol']))
            if 'HGNC_RefSeq_NM' in line and line['HGNC_RefSeq_NM']:
                self.p('Replaced HGNC_RefSeq_NM %s with %s' % (line['HGNC_RefSeq_NM'], refseq))
            line['HGNC_RefSeq_NM'] = refseq

            yield line

    def download_mim2gene(self):
        """Download the mim2gene.txt file from omim ftp server. By default the file
        will be downloaded to the location of the script.

        Returns (str): full path and filename of the mim2gene.txt
        """
        filename = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + 'mim2gene.txt'
        (dl_filename, headers) = urlretrieve('http://omim.org/static/omim/data/mim2gene.txt',
                                             filename)

        return dl_filename

    def query_omim(self, data):
        """Queries OMIM to fill in the inheritance models

        Args:
                data (list of dicts): Inner dict represents a row in a gene list

        Yields:
                dict: with the added HGNC symbol prepended to the HGNC_symbol column.
        """
        omim = OMIM(api_key='<fill in key>')
        for line in data:
            if 'OMIM_morbid' in line and 'Chromosome' in line:
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
                line['Phenotypic_disease_model'] = '%s:%s' % \
                    (line['HGNC_symbol'], '|'.join(line_phenotypic_disease_models))
            else:
                line['Phenotypic_disease_model'] = ''
            # add OMIM morbid
            if entry['mim_number'] is not None:
                if 'OMIM_morbid' in line \
                and len(line['OMIM_morbid']) > 0 \
                and str(line['OMIM_morbid']) != line['HGNC_symbol']+':'+str(entry['mim_number']) \
                and str(line['OMIM_morbid']) != str(entry['mim_number']):
                    self.p('%s %s > %s client OMIM number differs from OMIM query' % \
                           (line['HGNC_symbol'], line['OMIM_morbid'], entry['mim_number']))
                line['OMIM_morbid'] = '%s:%s' % (line['HGNC_symbol'], entry['mim_number'])

            # add Gene_locus
            if entry['gene_location'] is not None:
                if 'Gene_locus' in line and len(line['Gene_locus']) > 0 and \
                   line['Gene_locus'] != entry['gene_location']:
                    self.p('%s > %s client Gene locus differs from OMIM query' % \
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

    def annotate(self, infile, outfile, verbose=False, errors=False, download_mim2gene=False, mim2gene=False, zero=False):
        """ Annotate a gene list """

        # make sure we print if we are asked to
        if verbose:
            self.verbose = True

        # show only the EnsEMBLdb conflicts - so verbose, but supress printing of the gene list
        if errors:
            self.verbose = True
            self.errors_only = True

        # read in the TSV file
        tsvfile = open(infile, 'r')
        self.outfile = open(outfile, 'w')
        raw_data = (line.strip() for line in tsvfile) # sluuuurp
        parsable_data = (line.split("\t") for line in raw_data)
        # download a new version of mim2gene.txt
        if download_mim2gene:
            self.p('Downloading mim2gene.txt ... ')
            dl_filename = self.download_mim2gene()
            self.p('Done: %s' % dl_filename)

        # check mem2gene.txt for HGNC symbol resolution
        if mim2gene:
            self.mim2gene = True
            self.cache_mim2gene()

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
        dict_data = self.list2dict(header, parsable_data)

        # clean up the input
        clean_data = self.cleanup(dict_data)
        fixed_data = None
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
        self.conn = pymysql.connect(host='localhost', port=3306,
                                    user='anonymous', db='homo_sapiens_core_75_37')
        ensembld_data = self.query(refseq_data, try_hgnc_again=True)

        # put the official HGNC symbol
        hgnc_official_data = self.put_official_hgnc_symbol(ensembld_data)

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

        # get all contigs
        final_data = self.gather_contig(cleaner_data)
        print_data = []
        for line in final_data:
            if self.verbose:
                print(self.format_line(line))
            print_data.append(line)

        # print the gene list
        for comment in comments:
            self.outfile.write('\t'.join(comment) + '\n')
        for line in self.get_contigs():
            self.outfile.write(line + '\n')
        self.outfile.write(self.get_header() + '\n')
        for line in print_data:
            self.outfile.write(self.format_line(line) + '\n')
