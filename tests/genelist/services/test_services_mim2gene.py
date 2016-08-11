from genelist.services.mim2gene import Mim2gene

def test_resolve_gene():
    mim2gene = Mim2gene()
    mim2gene.resolve_gene('FARS2')

def test_resolve_gene_id():
    pass
