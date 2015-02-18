#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function
import sys
import pymysql
import argparse
import re
import os
from time import sleep
from omim import OMIM
from urllib.request import urlretrieve, Request, urlopen

gl_header=['Chromosome', 'Gene_start', 'Gene_stop', 'HGNC_symbol', 'Protein_name', 'Symptoms', 'Biochemistry', 'Imaging', 'Disease_trivial_name', 'Trivial_name_short', 'Phenotypic_disease_model', 'OMIM_morbid', 'Gene_locus', 'UniProt_id', 'Ensembl_gene_id', 'Ensemble_transcript_ID', 'Reduced_penetrance', 'Clinical_db_gene_annotation', 'Disease_associated_transcript', 'Ensembl_transcript_to_refseq_transcript', 'Gene_description']

# EnsEMBL connection
# TODO make this prettier
conn = None

verbose = False # to print or not to print
errors_only = False # print only errors. Does not print the gene list. Needs verbose=True to work.
mim2gene = False # resolve HGNC symbol with mim2gene.txt

def print_header(header=gl_header):
    """Prints the gl_header

    Args:
        header (list, optional): a list of strings
    Returns:
        pass
    Note:
        Prints to STDOUT
    """
    if not errors_only:
        print('#' + "\t".join(gl_header))

def print_line(line):
    """Prints a line based on the order of the gl_headers

    Args:
        line (dict): dict with values for one line of a gene list. All values of gl_headers should be present as keys

    Returns:
        pass

    Note:
        Will print the STDOUT
    """
    if not errors_only:
        ordered_line = list()
        for column_name in gl_header:
            ordered_line.append(str(line[column_name]))
        print("\t".join(ordered_line))

def p(line, end=os.linesep):
    """print only if the verbose switch has been set

    Args:
        line (str): line to print to STDOUT

    Returns:
        pass
    """
    if verbose:
        print('\033[93m', '>>> ', line, '\033[0m', end=end)

def resolve_ensembl_id(hgnc_id):
    """Query genenames.org for the EnsEMBL gene id based on the HGNC symbol.

    Args:
        hgnc_id (str): the HGNC symbol

    Returns (str): The ensEMBL gene id

    """
    import json
    response = urlopen(Request("http://rest.genenames.org/fetch/symbol/%s" % hgnc_id, None, {'Accept':'application/json'}))
    data = response.read().decode('UTF-8')
    data = json.loads(data)
    try:
        return data['response']['docs'][0]['ensembl_gene_id']
    except KeyError as ke:
        return False

symbol_of = {} # omim_id: hgnc_symbol
type_of = {} # hgnc_symbol: type

# TODO: djees, put this in a separate package so we don't have to rely on a global var
def cache_mim2gene(mim2gene_file=os.path.dirname(os.path.abspath(__file__))+os.path.sep+'mim2gene.txt'):
    """Read in the mim2gene file and store it as a dict of OMIM id: HGNC_symbol. Only gene and gene/phenotype types will be saved.

    Kwargs:
        mim2gene_file (str): the aboslute path to the mim2gene.txt file

    Returns: None
    """
    global symbol_of
    mim2gene_fh = open(mim2gene_file, 'r')
    lines = ( line.strip() for line in mim2gene_fh )
    for line in lines:
        (file_omim_id, omim_type, gene_id, hgnc_symbol) = line.split("\t")
        if omim_type in ('gene', 'gene/phenotype') and hgnc_symbol != '-':
            symbol_of[file_omim_id] = hgnc_symbol
        type_of[hgnc_symbol] = omim_type

def resolve_gene(omim_id):
    """Looks up the omim_id in the mim2gene.txt file. If found and the omim type is 'gene', return the HGNC symbol

    Args:
        omim_id (int): the omim id

    Returns: on omom id match, HGNC symbol if type of gene or gene/phenotype otherwise False

    """
    global symbol_of
    if omim_id in symbol_of:
        return symbol_of[omim_id]
    return False

def query(data, try_hgnc_again=False):
    """Queries EnsEMBL. Parameters are HGNC_symbol and/or Ensembl_gene_id, whatever is available. Data from EnsEMBLdb will overwrite the client data.
    A(n) identifier(s) should yield one result from EnsEMBLdb. It will be reported if a(n) identifier(s) don't yield any or multiple results.

    Args:
        data (list of dicts): representing the lines and columns in a gene list. The keys of the dicts must match the column names of the EnsEMBLdb query.
        try_hgnc_again (bool): when providing multiple HGNC ids, try until you find a match on EnsEMBLdb. Only one HGNC ID will be used in final result.
    Yields:
        dict: a row with data from ensEMBLdb filled in.
    """
    conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')
    cur = conn.cursor(pymysql.cursors.DictCursor)

    base_query = """
    SELECT g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop,
    x.display_label AS HGNC_symbol, g.stable_id AS Ensembl_gene_id,
    seq_region.name AS Chromosome
    FROM gene g JOIN xref x ON x.xref_id = g.display_xref_id
    join seq_region USING (seq_region_id)
    """
    keys_conds = { 'HGNC_symbol': 'x.display_label', 'Ensembl_gene_id': 'g.stable_id', 'Chromosome': 'seq_region.name' }

    keys = ['HGNC_symbol', 'Ensembl_gene_id', 'Chromosome'] # these columns will be put into the condition statement if they have a value
    for line in data:
        HGNC_symbols=line['HGNC_symbol'].split(',')
        HGNC_symbol_i=1
        for HGNC_symbol in HGNC_symbols:
            line['HGNC_symbol'] = HGNC_symbol # actually replace the entry
            conds = [ "%s = %%s" % keys_conds[ key ] for key in keys if key in line and line[key] != None and len(line[key]) > 0 ]
            cond_values = [ line[ key ] for key in keys if key in line and line[key] != None and len(line[key]) > 0 ]

            # check on length of the region name to exclude scaffolds and patches
            query = "%s where length(seq_region.name) < 3 and %s" % ( base_query, " and ".join(conds) )
            cur.execute(query, cond_values)

            rs = cur.fetchall() # result set

            if len(rs) == 0:
                if HGNC_symbol_i == len(HGNC_symbols):
                    not_found_id = HGNC_symbol if len(HGNC_symbols) == 1 else HGNC_symbols
                    p("Not found: %s %s" % (not_found_id, cond_values))
                if not try_hgnc_again: break
            elif len(rs) > 1:
                if HGNC_symbol_i > 1:
                    p("Took %s/%s" % (HGNC_symbol, HGNC_symbols))

                # This happens when multiple genes are overlapping with a HGNC ID
                # Query genenames.org to pick one EnsEMBL gene id
#                gn_ensembl_id = resolve_ensembl_id(HGNC_symbol)
#                if gn_ensembl_id:
#                    matching_lines = [ line for line in rs if line['Ensembl_gene_id'] == gn_ensembl_id ] # yields one entry
#                    if len(matching_lines) > 0:
#                        p("Took %s" % gn_ensembl_id)
#                        for entry in matching_lines:
#                            yield merge_line(entry, line)
#                        break

                # we couldn't resolve this with genenames.org
                if 'Chromosome' in line:
                    p("Multiple entries: %s, chromosome: %s => " % (HGNC_symbol, line['Chromosome']), end='')
                else:
                    p("Multiple entries: %s => " % (HGNC_symbol), end='')
                p("Adding: %s" % ', '.join(( entry['Ensembl_gene_id'] for entry in rs )) )
                for entry in rs:
                    yield merge_line(entry, line)
                break
            else:
                if len(HGNC_symbols) > 1:
                    p("Took %s/%s" % (HGNC_symbol, HGNC_symbols))
                for entry in rs:
                    yield merge_line(entry, line)
                break
            HGNC_symbol_i += 1

def query_transcripts(data):
    """Queries EnsEMBL for all transcripts.

    Args
        data (list of dicts): representing the lines and columns in a gene list. The keys of the dicts must match the column names of the EnsEMBLdb query.
    Yields (dict):
        a row with transcript data from ensEMBLdb filled in.

    """
    global conn
    cur = conn.cursor(pymysql.cursors.DictCursor)

    def _cleanup_description(description):
        """Remove the comment in the description and clean up invalid characters: ,:;|>

        Args:
            description (str): text to clean up

        Returns: str or None

        """
        if description:
            description = description.strip()
            description = re.sub(r'\[.*\]', '', description)
            description = re.sub(r'[,:;>| ]', '_', description)
            if description.endswith('_'):
                description = description[:-1]
            return description
        return None

    def _process_transcripts(data):
        """Processes raw data:
        * aggregates transcripts, RefSeq IDs

        Args:
            data (dict): dictionary with following keys: EnsEMBL_ID, description, Transcript_ID, RefSeq_ID

        yields (str): A string with transcripts, RefSeq IDs aggregated

        """
        row = data.pop(0)

        # init
        Ensembl_ID = row['Ensembl_ID']
        line = {} # will only hold two keys: Ensembl_transcript_to_refseq_transcript and Gene_description
        prev_description = _cleanup_description(row['description'])
        transcripts = ['%s>%s' % (row['Transcript_ID'], row['RefSeq_ID'])]

        for row in data:
            if row['Ensembl_ID'] != Ensembl_ID:

                if len(transcripts) == 0:
                    p(Ensembl_ID + ' has no transcripts!')

                line['Ensembl_transcript_to_refseq_transcript'] = '%s:%s' % (Ensembl_ID, '|'.join(transcripts))
                line['Gene_description'] = prev_description

                yield line

                # reset
                transcripts = []
                Ensembl_ID = row['Ensembl_ID']

                line = {}
                prev_description = _cleanup_description(row['description'])

            if row['RefSeq_ID'] == None:
                p('%s:%s has no RefSeqID' % (Ensembl_ID, row['Transcript_ID']))

            transcripts.append('%s>%s' % (row['Transcript_ID'], row['RefSeq_ID']))

        # yield last one
        line['Ensembl_transcript_to_refseq_transcript'] = '%s:%s' % (Ensembl_ID, '|'.join(transcripts))
        line['Gene_description'] = prev_description
        yield line

    """
    external_db_id = 1801
    select * from xref where display_label like 'NM\_%' limit 10;
    """

    for line in data:
        base_query = """
        SELECT g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop,
        x.display_label AS HGNC_symbol, g.stable_id AS Ensembl_ID,
        seq_region.name AS Chromosome, t.stable_id AS Transcript_ID, g.description,
        tx.dbprimary_acc AS RefSeq_ID
        FROM gene g
        JOIN xref x ON x.xref_id = g.display_xref_id
        JOIN seq_region USING (seq_region_id)
        LEFT JOIN transcript t ON t.gene_id = g.gene_id
        LEFT JOIN object_xref ox ON ox.ensembl_id = t.transcript_id
        LEFT JOIN xref tx ON tx.xref_id = ox.xref_id
        WHERE length(seq_region.name) < 3
        AND tx.external_db_id in (1801, 1806, 1810)
        AND g.stable_id = %s
        """

        cur.execute(base_query, line['Ensembl_gene_id'])
        rs = cur.fetchall()
        transcripts = _process_transcripts(rs)

        for transcript in transcripts:
            line['Ensembl_transcript_to_refseq_transcript'] = transcript['Ensembl_transcript_to_refseq_transcript']
            line['Gene_description'] = transcript['Gene_description']

        yield line

def get_transcript(start, end, ensembl_gene_id=None, hgnc_id=None):
    """Queries EnsEMBL. Parameters are HGNC_symbol and/or Ensembl_gene_id, whatever is available. It will return one hit with the ensembl trasncript id.

    Args:
        ensembl_gene_id (str, optional): ensembl_gene_id and/or hgnc_id should be provided.
        hgnc_id (str, optional): ensembl_gene_id and/or hgnc_id should be provided.
        start (int): start coordinate of the possible transcript of this gene.
        end (int): stop coordinate of the possible transcript of this gene.

    Yields:
        dict: with following keys: ensembl_gene_id, hgnc_id, start, end, transcript_id
    """
    cur = conn.cursor(pymysql.cursors.DictCursor)

    base_query = "select g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop, x.display_label AS HGNC_symbol, g.stable_id AS Ensembl_gene_id, seq_region.name AS Chromosome, t.seq_region_start AS Transcript_start, t.seq_region_end AS Transcript_stop, t.stable_id AS Transcript_id from gene g join xref x on x.xref_id = g.display_xref_id join seq_region using (seq_region_id) join transcript t using (gene_id)"

    conds = {'t.seq_region_start': start, 't.seq_region_end': end}
    if ensembl_gene_id != None:
        conds.update({'g.stable_id': ensembl_gene_id})
    if hgnc_id != None:
        conds.update({'x.display_label': hgnc_id})

    query = "%s where %s" % ( base_query, " and ".join([ '%s = %%s' % column for column in conds.keys() ]) )
    cur.execute(query, [ str(value) for value in conds.values() ])

    rs = cur.fetchall() # result set

    # O-oh .. for now this still means manual intervention!
    if len(rs) > 1:
        p("Getting '%s' ... " % conds.values())
        p('Multiple entries found!')
    elif len(rs) == 0:
        return None
    return rs[0]

def merge_line(ens, client):
    """Will merge dict ens (EnsEMBL data) with client (data). ens will take precedence over client. Changes will be reported.

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
                p("%s > %s: ens '%s' diff from client '%s'" % (ens['Ensembl_gene_id'], key, ens[key], value))
        #    else:
        #        p("%s: ens '%s' eq to client '%s'" % (key, ens[key], value))
        #else:
        #    p("%s not in ens!" % key)

    # Check the Gene_start and Gene_stop for being Transcript coordinates
    if 'Gene_start' in client and 'Gene_stop' in client:
        if str(ens['Gene_start']) != str(client['Gene_start']) or str(ens['Gene_stop']) != str(client['Gene_stop']):
            # get the transcript, compare those coordinates and report
            transcript = get_transcript(client['Gene_start'], client['Gene_stop'], ens['Ensembl_gene_id'], ens['HGNC_symbol'])
            for key in ('Gene_start', 'Gene_stop'):
                if str(ens[key]) != str(client[key]):
                    if transcript and len(transcript) > 1:
                        p("%s > %s: ens '%s' diff from client '%s', but matches Transcript %s: %s" % (ens['Ensembl_gene_id'], key, ens[key], client[key], transcript['Transcript_id'], transcript[key.replace('Gene', 'Transcript')]))
                    else:
                        p("%s > %s: ens '%s' diff from client '%s'" % (ens['Ensembl_gene_id'], key, ens[key], client[key]))

    merged = client.copy()
    merged.update(ens)
    return merged

def fill(data):
    """Removes #NA's

    Args:
        data (list of dicts): representing the lines and columns in a gene list
    Yields:
        dict: with all missing columns filled in with ''

    """
    defaults=dict((column_name, '') for column_name in gl_header)
    for line in data:
        d = defaults.copy()
        d.update(line)
        yield d

def munge(data):
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

def cleanup(data):
    """Will clean the data according to rules set by MIP
        # replace white space seperated comma's with just a comma
        # replace ; with comma
        # remove leading and trailing white space

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
                    line[key] = re.sub(r'\s*[,;]\s*', ',', value.strip())
        yield line

def list2dict(header, data):
    """Will convert each row in the data from a list to dict using the header list as keys.

    Args:
        header (list): A list containing the keys for the dict generation
        data (list of lists): Inner list represents a row in a gene list
    Yields:
        dict: the next dictified line

    """
    for line in data:
        yield dict(zip(header, line))

def zero2one(data):
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

def add_mim2gene_alias(data):
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
            OMIM_id  = line['OMIM_morbid']
            HGNC_symbol = resolve_gene(OMIM_id)
            if HGNC_symbol != False and HGNC_symbol not in HGNC_symbols:
                HGNC_symbols.insert(0, HGNC_symbol)
        line['HGNC_symbol'] = ','.join(HGNC_symbols)
        yield line

def add_genome_build(data, genome_build):
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

def redpen2symbol(data):
    """If reduced penetrance is set, replace it with the HGNC symbol

    Args:
        data (list of dicts): Inner dict represents a row in a gene list

    Yields:
        dict: with the replaced red pen to HGNC symbol

    """
    for line in data:
        if 'Reduced_penetrance' in line and line['Reduced_penetrance'] == 'Yes':
            line['Reduced_penetrance'] = line['HGNC_symbol']
        yield line

def remove_non_genes(data):
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
                if HGNC_symbol in type_of and type_of[ HGNC_symbol ] in ('gene', 'gene/phenotype'):
                    yield line
                    yielded = True
                    break
            if not yielded:
                p('Removed: %s' % line)

def download_mim2gene():
    """Download the mim2gene.txt file from omim ftp server. By default the file
    will be downloaded to the location of the script.


    Returns (str): full path and filename of the mim2gene.txt
    """
    filename = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + 'mim2gene.txt'
    (dl_filename, headers) = urlretrieve('ftp://anonymous:kennybilliau%40scilifelab.se@ftp.omim.org/OMIM/mim2gene.txt', filename)
    return dl_filename

def query_omim(data):
    """Queries OMIM to fill in the inheritance models

    Args:
        data (list of dicts): Inner dict represents a row in a gene list

    Yields:
        dict: with the added HGNC symbol prepended to the HGNC_symbol column.
    """
    TERMS_MAPPER = {
      'Autosomal recessive': 'AR',
      'Autosomal dominant': 'AD',
      'X-linked dominant': 'XD',
      'X-linked recessive': 'XR',
    }

    TERMS_BLACKLIST = [
      'Isolated cases',
    ]

    omim = OMIM(api_key='<fill in key>')
    for line in data:
        if 'HGNC_symbol' in line:
            entry = omim.gene(line['HGNC_symbol'])

            phenotypic_disease_model = []
            models = set()
            for phenotype in entry['phenotypes']:
                if phenotype['inheritance'] is not None:
                    models.update([model.strip('? ') for model in phenotype['inheritance'].split(';')])
                    models = models.difference(TERMS_BLACKLIST) # remove blacklisted terms
                    models = set([TERMS_MAPPER.get(model_human, model_human) for model_human in models]) # rename them if possible

                    phenotypic_disease_model.append('%s>%s' % (phenotype['phenotype_mim_number'], '/'.join(models)))
                else:
                    if (phenotype['phenotype_mim_number'] is not None):
                        phenotypic_disease_model.append(str(phenotype['phenotype_mim_number']))

            if len(phenotypic_disease_model) > 0:
                line['Phenotypic_disease_model'] = '%s:%s' % (line['HGNC_symbol'], '|'.join(phenotypic_disease_model))

            # add OMIM morbid
            if entry['mim_number'] is not None:
                if 'OMIM_morbid' in line and len(line['OMIM_morbid']) > 0 and line['OMIM_morbid'] != entry['mim_number']:
                    p('%s > %s client OMIM number differs from OMIM query' % (line['OMIM_morbid'], entry['mim_number']))
                line['OMIM_morbid'] = '%s:%s' % (line['HGNC_symbol'], entry['mim_number'])

            # add Gene_locus
            if entry['gene_location'] is not None:
                if 'Gene_locus' in line and len(line['Gene_locus']) > 0 and line['Gene_locus'] != entry['gene_location']:
                    p('%s > %s client Gene locus differs from OMIM query' % (line['Gene_locus'], entry['gene_location']))
                line['Gene_locus'] = entry['gene_location']

            sleep(0.25) # wait for 250ms as according to OMIM specs
        yield line

def main(argv):
    # set up the argparser
    parser = argparse.ArgumentParser(description='Queries EnsEMBL and fills in the blanks of a gene list. Only columns headers found in gl_headers will be used')
    parser.add_argument('infile', type=argparse.FileType('r'), help='the tsv file with correct headers')
    parser.add_argument('--zero', default=False, action='store_true', dest='zero_based', help='if set, will convert 0-based coordinates to 1-based')
    parser.add_argument('--verbose', default=False, action='store_true', dest='verbose', help='if set, will show conflict messages from EnsEMBLdb inbetween the gene list lines')
    parser.add_argument('--errors-only', default=False, action='store_true', dest='errors_only', help='if set, will not output the gene list, but only the conflict messages from EnsEMBLdb.')
    parser.add_argument('--download-mim2gene', default=False, action='store_true', dest='download_mim2gene', help='if set, will download a new version of the mim2gene.txt file, used to check the OMIM type')
    parser.add_argument('--mim2gene', default=False, action='store_true', dest='mim2gene', help='if set, will try to resolve an HGNC symbol with the help of mim2gene.txt')
    parser.add_argument('--genome-build', default=None, dest='genome_build', help='Sets the genome release version, e.g. GRCh37')
    args = parser.parse_args(argv)

    global verbose
    # make sure we print if we are asked to
    if args.verbose:
        verbose = True

    # show only the EnsEMBLdb conflicts - so verbose, but supress printing of the gene list
    if args.errors_only:
        global errors_only
        verbose = True
        errors_only = True

    # download a new version of mim2gene.txt
    if args.download_mim2gene:
        p('Downloading mim2gene.txt ... ', end='')
        dl_filename = download_mim2gene()
        p('Done: %s' % dl_filename)

    # check mem2gene.txt for HGNC symbol resolution
    if args.mim2gene:
        global mim2gene
        mim2gene = True
        cache_mim2gene()

    # read in the TSV file
    tsvfile = args.infile
    raw_data = ( line.strip() for line in tsvfile ) # sluuuurp
    parsable_data = ( line.split("\t") for line in raw_data )

    # skip parsing of leading comments
    comments = []
    line = next(parsable_data)
    while line[0].startswith('##'):
        comments.append(line)
        line = next(parsable_data)

    # list to dict
    header = line # get the header
    dict_data = list2dict(header, parsable_data)

    # clean up the input
    clean_data = cleanup(dict_data)

    # add genome build, if any
    if args.genome_build is not None:
        genome_data = add_genome_build(clean_data, args.genome_build)
    else:
        genome_data = clean_data

    fixed_data = None
    if args.zero_based:
        # fix 0-based coordinates to be 1-based
        fixed_data = zero2one(genome_data)
    else:
        fixed_data = genome_data

    # add the mim2gene alias to HGNC_symbol
    # remove non-genes based on mim2gene.txt
    if mim2gene:
        aliased_data = add_mim2gene_alias(fixed_data)
        reduced_data = remove_non_genes(aliased_data)
    else:
        reduced_data = fixed_data

    # fill in missing blanks
    global conn
    conn = pymysql.connect(host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37')
    ensembld_data = query(reduced_data, try_hgnc_again=True)

    # aggregate transcripts
    transcript_data = query_transcripts(ensembld_data)

    # fill in the inheritance models
    omim_data = query_omim(transcript_data)

    # do some replacements
    redpen_data = redpen2symbol(omim_data)

    # fill in missing values with ''
    completed_data = fill(redpen_data)

    # clean up the data from EnsEMBL a bit
    munged_data = munge(completed_data)

    # at last, clean up the output
    cleaner_data = cleanup(munged_data)

    # print the gene list
    for comment in comments:
        print('\t'.join(comment))
    print_header()
    for line in cleaner_data:
        print_line(line)

if __name__ == '__main__':
    main(sys.argv[1:])
