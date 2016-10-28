from genelist.modules.merge import merge_panels

def test_annotate():
    merge_file_1 = open('tests/fixtures/merge-1.txt', 'r')
    merge_file_2 = open('tests/fixtures/merge-2.txt', 'r')
    merged_file  = open('tests/fixtures/merged.txt', 'r')

    merged_file_lines = [ line.rstrip('\n') for line in merged_file ]

    # skip the database line as it will never match
    merged_file_lines = merged_file_lines[1:]

    infiles = [ merge_file_1, merge_file_2 ]
    databases = ('ID', 'OMIM')
    merge_panels_out = []
    for line in merge_panels(infiles, databases):
        merge_panels_out.append(line)

    # skip the database line as it will never match
    merge_panels_out = merge_panels_out[1:]

    assert merged_file_lines == merge_panels_out
