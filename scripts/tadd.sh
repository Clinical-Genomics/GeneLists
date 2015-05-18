#!/bin/bash

# will add a gene list to a repo
# tag it with a version
# and retrieve the version + branch of the software it was generated with

# exit on errr
set -e

if [[ ${#@} < 2 ]]; then
    echo "USAGE: $0 genelist repo --minor|major"
    echo "	$0 ~/git/cust000 --minor"
    exit 1
fi

GENELIST=$1
TAG_BUMP=$2
OLD_WD=$(pwd)

SCRIPT_PATH=$(dirname $(readlink -nm $0))
GL_PATH=$(dirname $GENELIST)
GENELIST_NAME=$(basename $GENELIST)

# create the tag
cd $GL_PATH
TAG=$(git describe --abbrev=0 | tail -1 2> /dev/null)
case "$TAG_BUMP" in
    --minor) 
        MAJORPART=${TAG//.*}
        MINORPART=${TAG//*.}
        TAG=$(python -c "print('%d.%d' % (${MAJORPART}, ${MINORPART}+1))")
        ;;
    --major) TAG=$(python -c "import math; print('%.1f' % math.floor(${TAG}+1.0))")
        ;;
    *) >&2 echo "Invalid option"
       exit 1
       ;;
esac

# get current software version and branch of generated software repo
cd $SCRIPT_PATH
VERSION=$(git describe --tags | tail -1 2> /dev/null)
BRANCH=$(git symbolic-ref --short HEAD 2>/dev/null || git rev-parse --short HEAD 2>/dev/null)
cd $OLD_WD

# say something
read -p "Commit message: " MSG

cd "$(dirname $GENELIST)"

# get all panels
PANELS=( $(python $SCRIPT_PATH/get_panels.py $GENELIST) )
echo ${PANELS[@]}
TMP_GL=$(mktemp)

# remove the previous meta data headers
if [[ ${GENELIST_NAME} != 'cust000-Clinical_master_list.txt' ]]; then
    grep -v '^##' ${GENELIST} > ${TMP_GL} && mv ${TMP_GL} ${GENELIST}

    # add a header to the commited gene list
    for PANEL in "${PANELS[@]}"; do

        # we cat to avoid a positive exit code
        COMPLETE_NAME=$(grep -h "^${PANEL}:" $GL_PATH/LISTS | cat)
        COMPLETE_NAME=${COMPLETE_NAME#*: }

        if [[ ! -z "$COMPLETE_NAME" ]]; then
            echo "##Database=<ID=${GENELIST_NAME},Version=${TAG},Date=$(date +'%Y%m%d'),Acronym=${PANEL},Complete_name=${COMPLETE_NAME},Clinical_db_genome_build=GRCh37.p13" | cat - ${GENELIST} > ${TMP_GL} && mv ${TMP_GL} ${GENELIST}
        else
            echo "##Database=<ID=${GENELIST_NAME},Version=${TAG},Date=$(date +'%Y%m%d'),Acronym=${PANEL},Clinical_db_genome_build=GRCh37.p13" | cat - ${GENELIST} > ${TMP_GL} && mv ${TMP_GL} ${GENELIST}
        fi
    done
fi

# add the version to a changelog
DATE=$(date +"%y/%m/%d %H:%M")
echo "$DATE :: Generated with version $BRANCH:$VERSION" > CHANGELOG

# commit + tag
git add "$GENELIST_NAME"
git add CHANGELOG
git commit -m "$MSG"
git tag -a "$TAG" -m "$MSG"
git push
git push --tags origin

# update clinicalgenomics.se
cd $SCRIPT_PATH/..
python -m scripts.update_cg $(dirname $(dirname $(readlink -nm $GENELIST)))/cust00[1234]/cust*.txt > ~/git/clinical-genomics.github.io/_pages/namnpagenlistor.md
cd ~/git/clinical-genomics.github.io
git pull
git add _pages/namnpagenlistor.md
git commit -m "Update to $(basename $GENELIST)"
git push
cd $OLD_WD
