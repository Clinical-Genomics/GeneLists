#!/bin/bash

# Creates one master list based from multiple genelists.
#
# usage:
#    merge_lists.sh path-to-repo [outputfile [database name]]
#

# exit on errr
set -e

cleanup() {
    if [[ -e $PROCESSDIR ]]; then
        rm -Rf $PROCESSDIR
    fi
}
trap cleanup EXIT

REPODIR=$1
OUTPUTFILE=${2-'FullList.txt'}
FULLLISTNAME=${3-'FullList'} # combine multiple gene lists into one big list with this name

# copy them to easily work with them - CD to $REPODIR/original/ClinGeneLists/src:
echo -n "Copying lists ..."
PROCESSDIR=`mktemp -d`
cd $REPODIR
for GENELIST in $REPODIR/cust???/cust*txt
do
    cp $GENELIST $PROCESSDIR
done
echo "Lists are at '$PROCESSDIR'"

# cd to $PROCESSDIR - make sure the rest of the script runs on these files
cd $PROCESSDIR

# copy them all into one gene list file
echo -n "Copying to the one list ..."
FULLLISTALLCOLUMNS=`mktemp -p $PROCESSDIR`
cat * | sort > $FULLLISTALLCOLUMNS
echo "Done."

# pick the HGNC_ID and database column
echo -n "Picking HGNC_ID and Database column ..."
FULLLIST=`mktemp -p $PROCESSDIR`
cut -f4,21 $FULLLISTALLCOLUMNS > $FULLLIST
echo "Done."

# merge duplicated entries
echo -n "Merging duplicated entries ..."
PREVSYMBOL=''
DATABASES=()
TMPFILE=`mktemp -p $PROCESSDIR`
while read LINE; do
    IFS=$'\t' read -a LINE <<< "$LINE"
    if [[ "${LINE[0]}" == "$PREVSYMBOL" ]]; then
        IFS=$',' read -a CUR_DATABASES <<< ${LINE[1]}
        for DATABASE in ${CUR_DATABASES[*]}; do
            if [[ "$DATABASE" != "FullList" ]]; then
                DATABASES+=($DATABASE)
            fi
        done
    else
        if [[ -n $PREVSYMBOL ]]; then
            CUR_DATABASES=${DATABASES[*]} | sort -u
            DATABASES+=($FULLLISTNAME)
            DATABASE=$(IFS=,; echo "${DATABASES[*]}")
            echo "$PREVSYMBOL	$DATABASE" >> $TMPFILE
        fi
        DATABASES=(${LINE[1]})
    fi
    PREVSYMBOL=${LINE[0]}
done < $FULLLIST
DATABASES+=($FULLLISTNAME)
DATABASE=$(IFS=,; echo "${DATABASES[*]}")
echo "$PREVSYMBOL	$DATABASE" >> $TMPFILE

mv $TMPFILE $REPODIR/$OUTPUTFILE
echo "Done."

# add the right headers
echo -n "Adding the headers ..."
TMPFILE=`mktemp -p $PROCESSDIR`
echo "HGNC_ID	Database" | cat - $REPODIR/$OUTPUTFILE > $TMPFILE && mv $TMPFILE $REPODIR/$OUTPUTFILE
echo "Done."

# cleanup
trap - EXIT
cleanup
