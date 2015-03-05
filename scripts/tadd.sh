#!/bin/bash

# will add a gene list to a repo
# tag it with a version
# and retrieve the version + branch of the software it was generated with

# exit on errr
set -e

if [[ ${#@} < 2 ]]; then
    echo "USAGE: $0 genelist repo tag"
    echo "	$0 ~/bitbucket/GeneList/cust000 5.1"
    exit 1
fi

GENELIST=$1
TAG=$2
OLD_WD=$(pwd)

# get current software version and branch of generated software repo
SCRIPT_PATH=$(dirname $(readlink -nm $0))
cd $SCRIPT_PATH
VERSION=$(git describe | tail -1 2> /dev/null)
BRANCH=$(git symbolic-ref --short HEAD 2>/dev/null || git rev-parse --short HEAD 2>/dev/null)
cd $OLD_WD

# say something
read -p "Commit message: " MSG

cd "$(dirname $GENELIST)"

# add the version to a changelog
DATE=$(date +"%y/%m/%d %H:%M")
echo "$DATE :: Generated with version $BRANCH:$VERSION" > CHANGELOG

# commit + tag
git add "$(basename $GENELIST)"
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
#git push
cd $OLD_WD
