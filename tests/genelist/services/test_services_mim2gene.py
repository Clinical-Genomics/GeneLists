from genelist.services.mim2gene import Mim2gene

mim2gene = Mim2gene(download=True)

def test_resolve_gene():
    assert mim2gene.resolve_gene('611592') == 'FARS2' # FARS2
    assert mim2gene.resolve_gene('102300') is False # random phenotype

def test_resolve_gene_id():
    assert mim2gene.resolve_ensembl_gene_id('611592') == 'ENSG00000145982' # FARS
