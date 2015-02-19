#!/usr/bin/env python
# encoding: utf-8

import pymysql
import re

class Ensembl:

    def __init__(self, host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37'):
        self.conn = pymysql.connect(host=host, port=port, user=user, db=db)

    def __enter__(self, host='ensembldb.ensembl.org', port=5306, user='anonymous', db='homo_sapiens_core_75_37'):
        self.conn = pymysql.connect(host=host, port=port, user=user, db=db) # TODO find out how to combine init with enter
        return self

    def __exit__(self,type, value, traceback):
        self.conn.close()

    def query_transcripts(self, gene_id):
        """Queries EnsEMBL for all transcripts.

        Args
            gene_id: an ensembl gene id e.g. ENS00000124433
        Returns:
            dict: with keys Ensembl_transcript_to_refseq_transcript and Gene_description
                  Ensembl_transcript_to_refseq_transcript is formatted like this: ensembl_id:ensembl_transcript_id>ref_seq_id/ref_seq_id|

        """
        cur = self.conn.cursor(pymysql.cursors.DictCursor)

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

        def _join_refseqs(transcripts):
            transcripts_refseqs = []
            for transcript in sorted(transcripts.keys()):
                refseqs = '/'.join(sorted([ refseq for refseq in transcripts[ transcript ] if refseq != None ]))

                if len(refseqs) == 0:
                    transcripts_refseqs.append(transcript)
                else:
                    transcripts_refseqs.append('%s>%s' % (transcript, refseqs))

            return transcripts_refseqs

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
            transcripts = { row['Transcript_ID']: [ row['RefSeq_ID'] ] }

            for row in data:
                if row['Ensembl_ID'] != Ensembl_ID:

                    if len(transcripts) == 0:
                        p(Ensembl_ID + ' has no transcripts!')

                    line['Ensembl_transcript_to_refseq_transcript'] = '%s:%s' % (Ensembl_ID, '|'.join(_join_refseqs(transcripts)))
                    line['Gene_description'] = prev_description

                    yield line

                    # reset
                    transcripts = {}
                    Ensembl_ID = row['Ensembl_ID']

                    line = {}
                    prev_description = _cleanup_description(row['description'])

                if row['Transcript_ID'] not in transcripts:
                    transcripts[ row['Transcript_ID'] ] = []
                transcripts[ row['Transcript_ID'] ].append(row['RefSeq_ID'])

            # yield last one
            line['Ensembl_transcript_to_refseq_transcript'] = '%s:%s' % (Ensembl_ID, '|'.join(_join_refseqs(transcripts)))
            line['Gene_description'] = prev_description
            yield line

        """
        external_db_id = 1801
        select * from xref where display_label like 'NM\_%' limit 10;
        """

        base_query = """
        SELECT DISTINCT g.seq_region_start AS Gene_start, g.seq_region_end AS Gene_stop,
        x.display_label AS HGNC_symbol, g.stable_id AS Ensembl_ID,
        seq_region.name AS Chromosome, t.stable_id AS Transcript_ID, g.description,
        tx.dbprimary_acc AS RefSeq_ID
        FROM gene g
        JOIN xref x ON x.xref_id = g.display_xref_id
        JOIN seq_region USING (seq_region_id)
        LEFT JOIN transcript t ON t.gene_id = g.gene_id
        LEFT JOIN object_xref ox ON ox.ensembl_id = t.transcript_id AND ox.ensembl_object_type = 'Transcript'
        LEFT JOIN xref tx ON tx.xref_id = ox.xref_id AND tx.external_db_id in (1801, 1806, 1810)
        WHERE length(seq_region.name) < 3
        AND g.stable_id = %s
        ORDER BY g.gene_id, t.transcript_id
        """

        cur.execute(base_query, gene_id)
        rs = cur.fetchall()
        if len(rs) > 0:
            transcripts = _process_transcripts(rs)

            for transcript in transcripts:
                return {
                    'Ensembl_transcript_to_refseq_transcript': transcript['Ensembl_transcript_to_refseq_transcript'],
                    'Gene_description': transcript['Gene_description'],
                }
