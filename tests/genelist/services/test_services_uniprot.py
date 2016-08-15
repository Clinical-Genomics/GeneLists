from genelist.services.uniprot import Uniprot
from genelist.utils import cleanup_description

uniprot = Uniprot()

def test_resolve_gene():
    assert uniprot.fetch_description('O95363') == cleanup_description('Phenylalanine--tRNA ligase, mitochondrial')
