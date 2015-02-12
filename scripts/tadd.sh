#!/bin/bash

# will add a gene list to a repo
# tag it with a version
# and retrieve the version + branch of the software it was generated with


if [[ ${#@} < 3 ]]; then
    echo "USAGE:"
    exit 1
fi

GENELIST=$1
REPO=$2
TAG=$3
OLD_WD=$(pwd)

# get current software version and branch
VERSION=$(git describe | tail -1 2> /dev/null)
BRANCH=$(git symbolic-ref --short HEAD 2>/dev/null || git rev-parse --short HEAD 2>/dev/null)

read -p "Commit message: " MSG

cp "$GENELIST" "$REPO"
cd "$REPO"
git add "$GENELIST"
git commit -m $MSG
git tag -a "$TAG" -m $MSG

DATE=$(date +"%y/%m/%d %H:%M")
echo "$DATE :: SOFTWARE VERSION $BRANCH:$VERSION" >> CHANGELOG
