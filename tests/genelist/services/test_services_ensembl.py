from genelist.services.ensembl import Ensembl
import yaml

def init(config_stream):
    config = yaml.load(config_stream)
    return Ensembl(
        host=config['ensembl']['host'],
        port=config['ensembl']['port'],
        user=config['ensembl']['user'],
        db=config['ensembl']['db']
    )

def test_query_ensembl(config_stream):
    ensembl = init(config_stream)

    # QUERY
    # TRMT10A: found with HGNC not OMIM
    assert ensembl.query(omim_morbid='616013') == []
    assert ensembl.query(hgnc_symbol='TRMT10A') == [{'Chromosome': '4', 'Ensembl_gene_id': 'ENSG00000145331', 'Gene_start': 100467866, 'Gene_stop': 100485189}]

    # test with OMIM
    assert ensembl.query(omim_morbid='609300') == [{'Chromosome': '10', 'Ensembl_gene_id': 'ENSG00000148795', 'Gene_start': 104590288, 'Gene_stop': 104597290}]

    # test with EnsEBML gene id
    assert ensembl.query(ensembl_gene_id='ENSG00000136872') == [{'Chromosome': '9', 'Ensembl_gene_id': 'ENSG00000136872', 'Gene_start': 104182860, 'Gene_stop': 104198105}]

    # test multiple hits: PIK3R2] Multiple E! entries: ['ENSG00000268173', 'ENSG00000105647']
    assert ensembl.query(hgnc_symbol='PIK3R2') == [
        {'Gene_start': 18263968, 'Chromosome': '19', 'Gene_stop': 18288927, 'Ensembl_gene_id': 'ENSG00000268173'},
        {'Chromosome': '19', 'Ensembl_gene_id': 'ENSG00000105647', 'Gene_start': 18263928, 'Gene_stop': 18281350} 
    ]

def test_query_transcripts(config_stream):
    ensembl = init(config_stream)

    # TRMT10A: should return 7 transcripts!
    assert ensembl.query_transcripts(ensembl_gene_id='ENSG00000145331') == {
        'Chromosome': '4',
        'Ensembl_gene_id': 'ENSG00000145331',
        'Ensembl_transcript_to_refseq_transcript': 'ENST00000273962>NM_152292/XM_005263352|ENST00000394876|ENST00000394877>NM_001134665/NM_001134666|ENST00000455368|ENST00000507394|ENST00000514547|ENST00000515831',
        'Gene_description': 'tRNA_methyltransferase_10_homolog_A_(S._cerevisiae)',
        'Gene_start': 100467866,
        'Gene_stop': 100485189
    }
    
    # ZSWIM6: found with E! gene id not OMIM
    assert ensembl.query_transcripts(omim_morbid='615951') == {}
    assert ensembl.query_transcripts(ensembl_gene_id='ENSG00000130449') == {
        'Chromosome': '5',
        'Ensembl_gene_id': 'ENSG00000130449',
        'Ensembl_transcript_to_refseq_transcript': 'ENST00000252744>NM_020928',
        'Gene_description': 'zinc_finger__SWIM-type_containing_6',
        'Gene_start': 60628100,
        'Gene_stop': 60841997
    }

    # BMPR1B: found with OMIM and E! gene id
    assert ensembl.query_transcripts(omim_morbid='603248', ensembl_gene_id='ENSG00000138696') == {
        'Chromosome': '4',
        'Ensembl_gene_id': 'ENSG00000138696',
        'Ensembl_transcript_to_refseq_transcript': 'ENST00000264568>NM_001256794|ENST00000394931|ENST00000440890|ENST00000502683|ENST00000506363|ENST00000509540>NM_001256793|ENST00000512312>NM_001256792|ENST00000515059>NM_001203/XM_005263180/XM_005263181',
        'Gene_description': 'bone_morphogenetic_protein_receptor__type_IB',
        'Gene_start': 95679119,
        'Gene_stop': 96079599
    }
