#!/bin/bash

# will add a gene list to a repo
# tag it with a version
# and retrieve the version + branch of the software it was generated with

# exit on errr
set -e

if [[ ${#@} < 2 ]]; then
    echo "USAGE: $0 genelist --minor|major"
    echo "	$0 ~/git/cust000 --minor"
    exit 1
fi

GENELIST=$1
TAG_BUMP=$2
OLD_WD=$(pwd)

SCRIPT_PATH=$(dirname $(readlink -nm $0))
GL_PATH=$(dirname $GENELIST)
GENELIST_NAME=$(basename $GENELIST)
GENELIST_SHORTNAME=${GENELIST_NAME%%.*} # remove .txt
GENELIST_SHORTNAME=${GENELIST_SHORTNAME#*-} # remove cust???-

# create the tag
cd $GL_PATH
# first find a tag of this list, if any
set +e
FULL_TAG=$(git tag --sort version:refname -l "${GENELIST_SHORTNAME}*" 2> /dev/null)
if [[ $? -ne 0 ]]; then # on fail, get the last tag
    TAG=$(git tag --sort version:refname | grep '^[0-9]\+\.[0-9]\+$' | tail -1 2> /dev/null)
    if [[ $? -ne 0 ]]; then
        TAG=0.0
    fi
else
    IFS=- read -a TAG_PARTS <<< "$FULL_TAG"
    TAG=${TAG_PARTS[1]}
    unset IFS
fi
set -e

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
#MSG="Update $GENELIST_NAME with OMIM"

cd "$(dirname $GENELIST)"

# get all panels
PANELS=( $(cd $SCRIPT_PATH/..; python -m genelist.cli panels $GENELIST) )
echo ${PANELS[@]}
TMP_GL=$(mktemp)

# remove the previous meta data headers
if [[ ${GENELIST_NAME} != 'cust000-Clinical_master_list.txt' ]]; then
    grep -v '^##Database' ${GENELIST} > ${TMP_GL} && mv ${TMP_GL} ${GENELIST}

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
else
    grep -v '^##Database=<ID=cust000-Clinical_master_list.txt' ${GENELIST} > ${TMP_GL} && mv ${TMP_GL} ${GENELIST}
    PANEL='FullList'
    COMPLETE_NAME=$(grep -h "^${PANEL}:" $GL_PATH/LISTS | cat)
    COMPLETE_NAME=${COMPLETE_NAME#*: }
    echo "##Database=<ID=${GENELIST_NAME},Version=${TAG},Date=$(date +'%Y%m%d'),Acronym=${PANEL},Complete_name=${COMPLETE_NAME},Clinical_db_genome_build=GRCh37.p13" | cat - ${GENELIST} > ${TMP_GL} && mv ${TMP_GL} ${GENELIST}
fi

# add the version to a changelog
DATE=$(date +"%y/%m/%d %H:%M")
echo "$DATE :: Generated with version $BRANCH:$VERSION" > VERSION

# commit + tag
git add "$GENELIST_NAME"
git add CHANGELOG
git add VERSION
git commit -m "$MSG"
git tag -a "${GENELIST_SHORTNAME}-${TAG}" -m "$MSG"
git push
git push --tags origin

# update clinicalgenomics.se
cd ~/git/clinical-genomics.github.io
git checkout source
cd -
cd $SCRIPT_PATH/..
python -m scripts.update_cg $(dirname $(dirname $(readlink -nm $GENELIST)))/cust00[01234]/cust*.txt > ~/git/clinical-genomics.github.io/_topics/namnpagenlistor.md
cd ~/git/clinical-genomics.github.io
git pull
git add _topics/namnpagenlistor.md
git commit -m "Update to $(basename $GENELIST)"
git push
cd $OLD_WD
