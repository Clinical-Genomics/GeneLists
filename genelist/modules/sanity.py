""" Methods to check validity of a gene list """
# encoding: utf-8

from __future__ import print_function
import sys
import re

class Sanity(object):
    """ Check the validty of a gene list. report back inconsistancies """

    def __init__(self):
        """TODO: Docstring for __init__.
        Returns: TODO

        """

        self.gl_header = ['Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_symbol', 'Protein_name',
                          'Symptoms', 'Biochemistry', 'Imaging', 'Disease_trivial_name',
                          'Trivial_name_short', 'Phenotypic_disease_model', 'OMIM_morbid',
                          'Gene_locus', 'UniProt_id', 'Ensembl_gene_id', 'Ensemble_transcript_ID',
                          'Reduced_penetrance', 'Clinical_db_gene_annotation',
                          'Disease_associated_transcript',
                          'Ensembl_transcript_to_refseq_transcript', 'Gene_description',
                          'Genetic_disease_model']
        self.mandatory_fields = {
            'Clinical_db_gene_annotation': re.compile(r'.+'),
            'Chromosome': re.compile(r'([\dXY]|MT)+'),
            'Gene_start': re.compile(r'\d+'),
            'Gene_stop': re.compile(r'\d+'),
            'HGNC_symbol': re.compile(r'.+'),
            # 'Genetic_disease_model': re.compile(r'.+'),
            # 'Gene_locus': re.compile(r'.+'),
            'Ensembl_gene_id': re.compile(r'ENSG\d{11}'),
        }

        # holds regex's with forbidden chars per column
        self.forbidden_chars = {
            # Can't have empty inheritance models:
            # Forbidden: SYMBOL:>, SYMBOL:OMIM||, SYMBOL:OMIM>, SYMBOL:OMIM>AR|>
            'Phenotypic_disease_model': re.compile(r'(:>|>$|:\d+>\||\|>)')
        }

        self.line_nr = 0 # hold on to the line nr
        self.warned = 0 # exit code

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

    def warn(self, msg):
        """Pretty print a warning message

        Args:
            msg (str): The message to print

        """
        print("#{}: {}".format(self.line_nr, msg))
        self.warned = 1

    def inc_line_nr(self, lines):
        """Increments the global line_nr for each passing line

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys
        Returns: pass

        """
        for line in lines:
            self.line_nr += 1
            yield line

    def check_forbidden_chars(self, lines):
        """Checks the precense of forbidden chars

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys

        yields: a line of the lines

        """
        for line in lines:
            for field, forbidden_re in self.forbidden_chars.items():
                if re.search(forbidden_re, line[field]):
                    self.warn("'{}' ('{}') has a forbidden char combination '{}'".\
                              format(field, line[field], forbidden_re))
            yield line

    def check_delimiter(self, lines, delimiters=re.compile(r', '), whitelist_keys=[]):
        """Check for the presence of a delimiter. Will use the 'warn' function to print presence of a delimiter

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys
            delimiters (re.compile): a compiled regex with what to match

        yields: a line of the lines

        """
        for line in lines:
            for key, value in line.items():
                if key in whitelist_keys: continue
                if re.search(delimiters, value):
                    self.warn("{} ('{}') holds a forbidden delimiter".format(key, value))
            yield line

    def check_trimming(self, lines):
        """Check if fields are trimmed

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys

        yields: a line of the lines
        """
        for line in lines:
            for key, value in line.items():
                if value != value.strip():
                    self.warn("{} ('{}') is not trimmed!".format(key, value))
            yield line

    def check_nr_fields(self, lines, header=None):
        """Checks if the nr of fields == nr of headers

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys

        yields: a line of the lines
        """
        if not header:
            header = self.gl_header
        len_header = len(header)
        for line in lines:
            len_line = len(line)
            if len_line != len_header:
                self.warn("Len fields ({}) != len headers ({})".format(len_line, len_header))
            yield line

    def check_mandatory_fields(self, lines):
        """Check if certain fields are filled in

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys

        yields: a line of the lines
        """
        for line in lines:
            for mandatory_field, mandatory_re in self.mandatory_fields.items():
                if not re.search(mandatory_re, line[mandatory_field]):
                    self.warn("'{}' ('{}') fails mandatory requirement '{}'".format(mandatory_field, line[mandatory_field], mandatory_re))
            yield line

    def check_coordinates(self, lines):
        """Checks if gene length is above zero

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys

        yields: a line of the lines
        """
        for line in lines:
            if line['Gene_stop'] and line['Gene_start']:
                if int(line['Gene_stop']) - int(line['Gene_start']) <= 0:
                    self.warn('Gene coordinates are not above zero.')
            yield line

    def check_duplicates(self, lines):
        """Check for duplicates based on HGNC_symbol, EnsEMBL_gene_id and coordinates

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys

        yields: a line of the lines
        """
        fields = {
            'HGNC_symbol': {}, # HGNC_symbol: line nr
            'Ensembl_gene_id': {},
        }
        Gene_start = {}
        Gene_stop = {}
        for line in lines:
            for field_key in fields.keys():
                if line[field_key] in fields[field_key]:
                    if (fields[field_key][line[field_key]] + 1) != self.line_nr:
                        self.warn("'{}' already listed at #{}".\
                                  format(line[field_key], fields[field_key][line[field_key]]))
                fields[field_key][line[field_key]] = self.line_nr

            if line['Gene_start'] in Gene_start and line['Gene_stop'] in Gene_stop:
                if (Gene_start[line['Gene_start']] + 1) != self.line_nr:
                    self.warn("'{}-{}' already listed at #{}".\
                              format(line['Gene_start'], line['Gene_stop'],
                                     Gene_start[line['Gene_start']]))
            Gene_start[line['Gene_start']] = self.line_nr
            Gene_stop[line['Gene_stop']]   = self.line_nr

            yield line

    def check_chromosome(self, lines):
        """Checks if the Gene_locus corresponds with the listed chromosome

        Args:
            lines (list of dicts): each dict contains a dict with gl_header as keys

        yields: a line of the lines
        """
        for line in lines:
            if 'Gene_locus' in line:
                # check if the locus has an q or p
                if any([x for x in line['Gene_locus'] if x in ('q', 'p')]):
                    gene_locus_chromosome = re.compile('p|q').split(line['Gene_locus'])[0]

                    if gene_locus_chromosome != line['Chromosome']:
                        self.warn("Chromosome '{}' differs from gene locus '{}'".\
                                  format(line['Chromosome'], line['Gene_locus']))

            yield line

    def check(self, infile):
        """Main program

        Args:
            infile (str): input gene list

        Returns: 0 on success, otherwise error code

        """
        infile = open(infile, 'r')
        lines = (line.strip('\r\n') for line in infile) # sluuuurp

        # get the line in parts
        parsable_data = (line.split('\t') for line in lines)

        # skip parsing of leading comments
        comments = []
        line = next(parsable_data)
        self.line_nr = 1
        while line[0].startswith('##'):
            comments.append(line)
            line = next(parsable_data)
            self.line_nr += 1

        # OK - start the sanity checks

        # header starts with #
        header = line # get the header
        if not header[0].startswith('#'):
            self.warn('Header does not start with #')
        else:
            header[0] = header[0][1:] # remove the leading #

        dict_data = self.list2dict(header, parsable_data)

        # are all mandatory headers there?
        header_diff = set(self.mandatory_fields.keys()).difference(header)
        if len(header_diff) > 0:
            self.warn("Missing columns {}".format(header_diff))

        # header should not contain white space
        for head in header:
            if re.search(' ', head) != None:
                self.warn("Header '{}' contains white space".format(head))

        # file contains ';' delimiter
        # elements in a field should be seperated by ',', not ', '
        # fields should be trimmed
        # #fields = #fields in header
        # following fields should be filled in:
        # - clinical_db_gene_annotation
        # - clinical_db_genome_build
        # - chromosome
        # - gene_start
        # - gene_stop
        # - HGNC_symbol
        # - Genetic_model
        # - Gene_locus
        # - Ensembl_gene_id
        # length of coordinates should be above zero
        # no duplicates based on coordinates tuple, EnsEMBL gene ID or HGNC symbol
        [line for line in \
              self.check_duplicates(
              self.check_coordinates(
              self.check_chromosome(
              self.check_mandatory_fields(
              self.check_forbidden_chars(
              self.check_trimming(
              self.check_delimiter(
              self.check_nr_fields(
              self.inc_line_nr(
                  dict_data
              ), header))))))))
        ]

        return self.warned
