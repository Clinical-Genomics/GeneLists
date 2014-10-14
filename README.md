# GeneLists

Scripts to automate the handling, populating and validating of gene lists


## Scripts

### fetch_info.py

```bash
usage: fetch_info.py [-h] infile

Queries EnsEMBL and fills in the blanks of a gene list. Only columns headers
found in gl_headers will be used.

positional arguments:
  infile      the tsv file with correct headers

optional arguments:
  -h, --help  show this help message and exit
```

This works by querying EnsEMBLdb and filling in missing information. The HGNC_ID and/or the Ensembl_gene_id are used as parameters in the query, whatever is available. Information coming from EnsEMBL will overwrite already present information in the infile, but it will be reported if doing so.

#### The Input File (infile)

This is a TSV file.
The infile should contain all the already available information for the gene list.
Each of the column headers should be one of the following:
Chromosome, Gene_start, Gene_stop, HGNC_ID, Disease_group_pathway, Protein_name, Symptoms, Biochemistry, Imaging, Disease_trivial_name, Trivial_name_short, Genetic_model, OMIM_gene, OMIM_morbid, Gene_locus, Genome_build, UniPort_ID, Ensembl_gene_id, Ensemble_transcript_ID, Red_pen, Database

Columns with other headers will be ignored and left out of the output.

### create_list.py

To autogenerate a genelist with only gene identifiers filled in. This was only used for testing purposes.

```bash
$ python create_list.py -h
usage: create_list.py [-h] [--hgnc] infile

Will output a genelist with the right right amount of column filled in with
#NA values. The right column with identifiers from the csv-infile will be
populated.

positional arguments:
  infile      a one column csv file with either HGNC or EnsEMBL IDs

optional arguments:
  -h, --help  show this help message and exit
  --hgnc      Indicate that the csv-infile is filled with HGNC identifiers
```

Example:

```
$ cat infile.csv
PMG3
GPD1L

```


```bash
$ python create_list.py --hgnc infile.csv

#Chromosome	Gene_start	Gene_stop	HGNC_ID	Disease_group_pathway	Protein_name	Symptoms	Biochemistry	Imaging	Disease_trivial_name	Trivial_name_short	Genetic_model	OMIM_gene	OMIM_morbid	Gene_locus	Genome_build	UniPort_ID	Ensembl_gene_id	Ensemble_transcript_ID	Red_pen	Database
#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	 PMG3 	#NA	#NA	#NA
#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	#NA	 GPD1L 	#NA	#NA	#NA
```


### find_hg_version.py

Script to figure out which versions of EnsEMBLdb are usable for a certain gene and gene coordinates. Would be used in a situation where a lot of coordinate mismatches are found between the client provided gene list and what is returned from EnsEMBLdb.

#### Usage
```bash
$ python find_hg_version.py 144989321 145050902 ENSG00000178209
Skipping homo_sapiens_core_48_36j with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_49_36k with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_50_36l with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_51_36m with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_52_36n with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_53_36o with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_54_36p with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_55_37 with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_56_37a with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_57_37b with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_58_37c with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_59_37d with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_60_37e with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_61_37f with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_62_37g with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_63_37 with Unknown column 'g.stable_id' in 'field list'
Skipping homo_sapiens_core_64_37 with Unknown column 'g.stable_id' in 'field list'
homo_sapiens_core_65_37 is the one
homo_sapiens_core_66_37 is the one
homo_sapiens_core_67_37 is the one
homo_sapiens_core_68_37 is the one
homo_sapiens_core_69_37 is the one
homo_sapiens_core_70_37 is the one
homo_sapiens_core_71_37 is the one
homo_sapiens_core_72_37 is the one
homo_sapiens_core_73_37 is the one
homo_sapiens_core_74_37 is the one
homo_sapiens_core_75_37 is the one
homo_sapiens_core_76_38 is NOT the one
homo_sapiens_core_77_38 is NOT the one
```
