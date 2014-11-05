#!/bin/bash

# exit on errr
set -e

source log.bash
log $(getversion)

MIPREFDIR=${1-'/mnt/hds/proj/bioinfo/mip/mip_references/'}
BITBUCKETDIR=${2-"${MIPREFDIR}/GeneLists/"}

CWD=`pwd`
for CUSTDIR in `ls -1 -d ${BITBUCKETDIR}/cust???`; do
    log "cd ${CUSTDIR}"
    cd "${CUSTDIR}"

    # pull the repo's
    log "git pull"
    git pull
    [[ -n $! ]] && exit 1 # if there is any conflict, exit

    for LIST in `ls -1 ${CUSTDIR}/*.txt`; do
        LIST=`basename "${LIST}"`
        # link the lists
        if [[ ! -e "${MIPREFDIR}/${LIST}" ]]; then
            log "ln -s \"${CUSTDIR}/${LIST}\" \"${MIPREFDIR}/${LIST}\""
            ln -s "${CUSTDIR}/${LIST}" "${MIPREFDIR}/${LIST}"
        fi

        LISTNAME=${LIST%.txt}
        # create the select*txt files
        SELECTFILE="${MIPREFDIR}/select_${LISTNAME}_vairants_db_master.txt"
        if [[ ! -e "${SELECTFILE}" ]]; then
            log "Creating ${SELECTFILE} file"
            echo "outinfo:Chromosome=>0_0,Variant_start=>0_1,Variant_stop=>0_2,Reference_allele=>0_3,Alternative_allele=>0_4,Ensemble_gene_id=>0_5,HGNC_symbol=>0_6,HGNC_approved_name=>0_7,HGNC_synonyms=>0_8,HGMD_accession=>0_9,HGMD_variant_type=>0_10,HGMD_variant_pmid=>0_11,Gene_annotation=>0_12,Functional_annotation=>0_13,HGNC_transcript_info=>0_14,Phast_cons_elements=>0_15,GERP_elements=>0_16,Genomic_super_dups=>0_17,Pseudogene=>0_18,1000G=>0_19,Dbsnp129=>0_20,Dbsnp_freq=>0_21,Dbsnp_nonflagged=>0_22,Dbsnp_rs_nr=>0_23,Esp6500=>0_24,HBVDB=>0_25,SIFT=>0_26,Poly_phen_hdiv=>0_27,Poly_phen_hvar=>0_28,Mutation_taster=>0_29,GERP=>0_30,LRT=>0_31,Phylo_p=>0_32,Transfac=>0_33,SnoRNA_miRNA_annotations=>0_34,Unscaled_C_score_1000G=>0_35,Scaled_C_score_1000G=>0_36,Unscaled_C_score_SNV=>0_37,Scaled_C_score_SNV=>0_38,Disease_group=>1_4,Clinical_db_genome_build=>1_15,Disease_gene_model=>1_11,Reduced_penetrance=>1_19,Clinical_db_gene_annotation=>1_20,GT_call_filter=>0_39,IDN!
ODF!/FDN!/ALIGNER!/GATK/FDN!_FILEENDING!_CALLTYPE!.txt	\t	5	Na	exact	0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,IDN!	small
RD!/${LIST}	\t	17	Na	exact	4,11,15,19,20	small" > "${SELECTFILE}"
        fi
    done
done
