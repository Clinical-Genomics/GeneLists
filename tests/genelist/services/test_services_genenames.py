from genelist.services.genenames import Genenames

genenames = Genenames()

def test_uniprot():
    assert genenames.uniprot('FARS2') == ['O95363']
    assert genenames.uniprot('GNAS') == ['O95467', 'P63092', 'P84996', 'Q5JWF2']

def test_refseq():
    assert genenames.refseq('FARS2') == ['NM_006567']

def test_aliases():
    assert genenames.aliases('FARS2') == ['dJ236A3.1']
    assert genenames.aliases('ALG9') == None

def test_official():
    assert genenames.official('FARS2') == 'FARS2'
    assert genenames.official('FARS2', '611592') == 'FARS2' # test with OMIM morbid 
    assert genenames.official('ALG9') == 'ALG9'
    assert genenames.official('DIBD1') == 'ALG9' # test previous symbol

def test_omim():
    assert genenames.omim('SLC2A1') == [138140, 143090]
    assert genenames.omim('APOA1BP') == None # NAXE alias
    assert genenames.omim('NAXE') == [608862]
