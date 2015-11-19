#!/bin/bash

# exit on errr
set -e

source log.bash
log $(getversion)

MIPREFDIR=${1-'/mnt/hds/proj/bioinfo/MIP_ANALYSIS/references'}
BITBUCKETDIR=${2-"/mnt/hds/proj/bioinfo/MIP_ANALYSIS/GeneLists/"}

SCRIPTDIR=`dirname $(readlink -f $0)`
CWD=`pwd`

for CUSTDIR in `ls -1 -d ${BITBUCKETDIR}/cust???`; do
    log "cd ${CUSTDIR}"
    cd "${CUSTDIR}"

    # pull the repo's
    log "git pull"
    git pull
    #[[ -n $? ]] && exit 1 # if there is any conflict, exit

    for LIST in `ls -1 ${CUSTDIR}/*.txt`; do

        # skip cust000
        LISTFILENAME=`basename "${LIST}"`
        CUSTNAME=${LISTFILENAME%-*}
        if [[ ${CUSTNAME} == 'cust000' ]]; then
            log "Skipping ${LIST}"
            continue
        fi

        # validate
        log "python ${SCRIPTDIR}/../scripts/sanity_check.py ${LIST}"
        python ${SCRIPTDIR}/../scripts/sanity_check.py ${LIST}
        
        # link the lists
        if [[ ! -e "${MIPREFDIR}/${LISTFILENAME}" ]]; then
            log "ln -s \"${CUSTDIR}/${LISTFILENAME}\" \"${MIPREFDIR}/${LISTFILENAME}\""
            ln -s "${CUSTDIR}/${LISTFILENAME}" "${MIPREFDIR}/${LISTFILENAME}"
        fi
    done
done
