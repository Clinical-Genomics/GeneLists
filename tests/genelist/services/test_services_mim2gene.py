from genelist.services.mim2gene import Mim2gene

mim2gene = Mim2gene(download=True)

def test_get_hgnc():
    assert mim2gene.get_hgnc('611592') == 'FARS2' # FARS2
    assert mim2gene.get_hgnc('102300') is False # random phenotype

def test_get_ensembl():
    assert mim2gene.get_ensembl('611592') == 'ENSG00000145982' # FARS
    assert mim2gene.get_ensembl('102300') is False # phenotype

def test_get_omim():
    assert mim2gene.get_omim('FARS2') == '611592'

def test_is_gene():
    assert mim2gene.is_gene('102300') is False # phenotype
    assert mim2gene.is_gene('102570') is False # removed
    assert mim2gene.is_gene('611592') is True  # gene
    assert mim2gene.is_gene('100650') is True  # gene/phenotype
    assert mim2gene.is_gene('ACTN1') is True   # gene
    assert mim2gene.is_gene('ACTN3') is True   # gene/phenotype

