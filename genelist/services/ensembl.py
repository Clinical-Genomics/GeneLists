#!/usr/bin/env python
# encoding: utf-8

import pymysql

from ..utils import cleanup_description

class Ensembl:

    def __init__(self, host='localhost', port=3306, user='anonymous', db='homo_sapiens_core_85_37'):
        self.conn = pymysql.connect(host=host, port=port, user=user, db=db)

    def __enter__(self, host='localhost', port=3306, user='anonymous', db='homo_sapiens_core_85_37'):
        self.conn = pymysql.connect(host=host, port=port, user=user, db=db) # TODO find out how to combine init with enter
        return self

    def __exit__(self, type, value, traceback):
        self.conn.close()

    def query(self, omim_morbid=None, ensembl_gene_id=None, hgnc_symbol=None, chromosome=None):
        """Queries EnsEMBL based on the Ensembl_gene_id. Data from EnsEMBLdb will overwrite
        the client data.
        An identifiers should yield one result from EnsEMBLdb.

        Args:
            ensembl_id (str): the EnsEMBL gene id.

        Yields (dict):
            { gene start,
            gene stop,
            chromosome,
            hgnc symbol }
            
        """

        # add 'x.display_label AS HGNC_symbol,' if oyu want to have the HGNC_symbol
        base_query = """
        SELECT DISTINCT g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop,
        g.stable_id AS Ensembl_gene_id,
        seq_region.name AS Chromosome
        FROM gene g JOIN xref x ON x.xref_id = g.display_xref_id
        join seq_region USING (seq_region_id)
        LEFT join object_xref ox on ox.ensembl_id = g.gene_id and ensembl_object_type = 'Gene'
        LEFT join xref xx on xx.xref_id = ox.xref_id and xx.external_db_id IN (1500, 1510, 1520)
        where length(seq_region.name) < 3
        """

        cond_values = []
        if omim_morbid:
            base_query += " AND xx.dbprimary_acc = %s"
            cond_values.append(str(omim_morbid))
        if ensembl_gene_id:
            base_query += " AND g.stable_id = %s"
            cond_values.append(ensembl_gene_id)
        if hgnc_symbol:
            base_query += " AND x.display_label = %s"
            cond_values.append(hgnc_symbol)
        if chromosome:
            base_query += " AND seq_region.name = %s"
            cond_values.append(chromosome)

        # execute the query
        cur = self.conn.cursor(pymysql.cursors.DictCursor)
        cur.execute(base_query, cond_values)
        rs = cur.fetchall() # result set

        if len(rs) == 0:
            return []
        else:
            return rs

    def query_transcripts_omim(self, omim_morbid=None, ensembl_gene_id=None):
        """Queries EnsEMBL for all transcripts.

        Args
            gene_id: an ensembl gene id e.g. ENS00000124433
        Returns:
            dict: with keys Ensembl_transcript_to_refseq_transcript and Gene_description
                  Ensembl_transcript_to_refseq_transcript is formatted like this:
                  HGNC_symbol:ensembl_transcript_id>ref_seq_id/ref_seq_id|

        """

        def _join_refseqs(transcripts):
            transcripts_refseqs = []
            for transcript in sorted(transcripts.keys()):
                refseqs = '/'.join(sorted([refseq for refseq in transcripts[transcript]
                                           if refseq != None]))

                if len(refseqs) == 0:
                    transcripts_refseqs.append(transcript)
                else:
                    transcripts_refseqs.append('%s>%s' % (transcript, refseqs))

            return transcripts_refseqs

        def _process_transcripts(data):
            """Processes raw data:
            * aggregates transcripts, RefSeq IDs

            Args:
                data (dict): dictionary with following keys: EnsEMBL_ID,
                             description, Transcript_ID, RefSeq_ID

            yields (str): A string with transcripts, RefSeq IDs aggregated

            """
            row = data.pop(0)

            # init
            ensembl_gene_id = row['Ensembl_gene_id']
            line = { # keys: Ensembl_transcript_to_refseq_transcript, Gene_description,
                     # Gene_start, Gene_stop, Chromosome, HGNC_symbol, Ensembl_gene_id
                'Gene_description': cleanup_description(row['description']),
                'Gene_start': row['Gene_start'],
                'Gene_stop': row['Gene_stop'],
                'Chromosome': row['Chromosome'],
                'Ensembl_gene_id': ensembl_gene_id
            }
            transcripts = {row['Transcript_ID']: [row['RefSeq_ID']]}

            for row in data:
                if row['Ensembl_gene_id'] != ensembl_gene_id:

                    line['Ensembl_transcript_to_refseq_transcript'] = \
                            '|'.join(_join_refseqs(transcripts))
                    yield line

                    # reset
                    transcripts = {}
                    ensembl_gene_id = row['Ensembl_gene_id']
                    line = {
                        'Gene_description': cleanup_description(row['description']),
                        'Gene_start': row['Gene_start'],
                        'Gene_stop': row['Gene_stop'],
                        'Chromosome': row['Chromosome'],
                        'Ensembl_gene_id': ensembl_gene_id
                    }

                if row['Transcript_ID'] not in transcripts:
                    transcripts[row['Transcript_ID']] = []
                transcripts[row['Transcript_ID']].append(row['RefSeq_ID'])

            # yield last one
            line['Ensembl_transcript_to_refseq_transcript'] = '|'.join(_join_refseqs(transcripts))
            yield line

        """
        external_db_id = 1801
        select * from xref where display_label like 'NM\_%' limit 10;
        """

        base_query = """
        SELECT DISTINCT g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop,
        g.stable_id AS Ensembl_gene_id, sr.name AS Chromosome,
        t.stable_id AS Transcript_ID, g.description, tx.dbprimary_acc AS RefSeq_ID
        FROM gene g JOIN xref x ON x.xref_id = g.display_xref_id
        JOIN seq_region sr ON sr.seq_region_id = g.seq_region_id
        LEFT JOIN object_xref ox on ox.ensembl_id = g.gene_id AND ensembl_object_type = 'Gene'
        LEFT JOIN xref xx on xx.xref_id = ox.xref_id AND xx.external_db_id IN (1500, 1510, 1520)
        LEFT JOIN transcript t ON t.gene_id = g.gene_id
        LEFT JOIN object_xref tox ON tox.ensembl_id = t.transcript_id AND tox.ensembl_object_type = 'Transcript'
        LEFT JOIN xref tx ON tx.xref_id = tox.xref_id AND tx.external_db_id in (1801, 1806, 1810)
        WHERE length(sr.name) < 3
        """

        cond_values = []
        if omim_morbid:
            base_query += " AND xx.dbprimary_acc = %s"
            cond_values.append(str(omim_morbid))
        if ensembl_gene_id:
            base_query += " AND g.stable_id = %s"
            cond_values.append(ensembl_gene_id)
            
        base_query += " ORDER BY g.gene_id, t.transcript_id"

        cur = self.conn.cursor(pymysql.cursors.DictCursor)
        cur.execute(base_query, cond_values)
        rs = cur.fetchall()
        if len(rs) > 0:
            transcripts = _process_transcripts(rs)
            return next(transcripts)
        return None
