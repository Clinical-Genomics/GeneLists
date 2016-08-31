from genelist.modules.genelist import Genelist

def test_annotate(config_stream):
    genelist = Genelist(config_stream)

    # read in the an annotated, albeit bare list create with the 'genelist merge' command
    bare_file = open('tests/fixtures/cust000-Clinical_master_list-bare.txt', 'r')
    bare_lines = [line for line in bare_file]

    # read in all result files
    annotated_file = open('tests/fixtures/cust000-Clinical_master_list.txt', 'r')
    annotated_lines = [line for line in annotated_file]

    annotated_verbose_file = open('tests/fixtures/cust000-Clinical_master_list-verbose.txt', 'r')
    annotated_verbose_lines = [line for line in annotated_verbose_file]

    mim2gene_file = open('tests/fixtures/cust000-Clinical_master_list-mim2gene.txt', 'r')
    mim2gene_lines = [line for line in mim2gene_file]

    mim2gene_verbose_file = open('tests/fixtures/cust000-Clinical_master_list-mim2gene-verbose.txt', 'r')
    mim2gene_verbose_lines = [line for line in mim2gene_verbose_file]

    # test if we get the same gene list with and without verbosity
    lines_out = [line + '\n' for line in genelist.annotate(bare_lines, info=False, warn=False,
                                                           error=False, download_mim2gene=False,
                                                           mim2gene=False, zero=False)]
    assert annotated_lines == lines_out

    lines_out = [line + '\n' for line in genelist.annotate(bare_lines, info=True, warn=True,
                                                           error=True, download_mim2gene=False,
                                                           mim2gene=False, zero=False)]
    assert annotated_verbose_lines == lines_out

    # test if we get the same gene list with and without mim2gene enabled
    lines_out = [line + '\n' for line in genelist.annotate(bare_lines, info=False, warn=False,
                                                           error=False, download_mim2gene=False,
                                                           mim2gene=True, zero=False)]
    assert mim2gene_lines == lines_out

    lines_out = [line + '\n' for line in genelist.annotate(bare_lines, info=True, warn=True,
                                                           error=True, download_mim2gene=False,
                                                           mim2gene=True, zero=False)]

    assert mim2gene_verbose_lines == lines_out
