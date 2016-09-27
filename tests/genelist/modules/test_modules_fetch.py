from genelist.modules.fetch import Fetch

def test_query_ensembl():
    pass

def test_annotate(config_stream):
    genelist = Fetch(config_stream, download_mim2gene=False)

    # read in the an annotated, albeit bare list create with the 'genelist merge' command
    cmms_file = open('tests/fixtures/cmms.txt', 'r')
    cmms_lines = [line for line in cmms_file]

    # read in the an annotated
    cmms_complete_file = open('tests/fixtures/cmms-complete.txt', 'r')
    cmms_complete_lines = [line for line in cmms_complete_file]

    lines_out = [line + '\n' for line in genelist.annotate(lines=cmms_lines)]
    #from pprint import pprint
    #pprint(lines_out)

    assert cmms_complete_lines == lines_out
