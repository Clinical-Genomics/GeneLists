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

uniqify() {
    echo "$1" | tr ' ' '\n' | sort -u | tr '\n' ' '
}


REPODIR=$1
OUTPUTFILE=${2-'FullList.txt'}
FULLLISTNAME=${3-'FullList'} # combine multiple gene lists into one big list with this name

# cd to $PROCESSDIR - make sure the rest of the script runs on these files
PROCESSDIR=`mktemp -d`
cd $PROCESSDIR

# copy them all into one gene list file
echo -n "Copying to the one list ..."
FULLLISTALLCOLUMNS=`mktemp -p $PROCESSDIR`
cat $REPODIR/cust???/cust*txt > $FULLLISTALLCOLUMNS
echo "Done."

# pick the HGNC_ID and database column
echo -n "Picking HGNC_ID, EnsEMBL Gene ID and Database column ..."
FULLLIST=`mktemp -p $PROCESSDIR`
cut -f4,18,21 $FULLLISTALLCOLUMNS > $FULLLIST
echo "Done."

# sort the list so we can merge duplicates in next step
# Also: remove all the headers from the files with sed :)
echo -n "Sorting ..."
FULLLISTSORTED=`mktemp -p $PROCESSDIR`
sort $FULLLIST | sed '/HGNC_ID/d' > $FULLLISTSORTED
echo "Done."

# merge duplicated entries
echo -n "Merging duplicated entries ..."
PREVSYMBOL=''
DATABASES=''
ENSIDS=() # hold on to all EnsEMBL IDs
TMPFILE=`mktemp -p $PROCESSDIR`
while read LINE; do
    IFS=$'\t' read -a LINE <<< "$LINE"
    if [[ "${LINE[0]}" == "$PREVSYMBOL" ]]; then
        DATABASES+=','${LINE[2]}
    else
        if [[ -n $PREVSYMBOL ]]; then
            ENSIDS=( $( uniqify "${ENSIDS[*]}" ) )
            for ENSID in ${ENSIDS[@]}; do
                DATABASES+=','$FULLLISTNAME
                IFS=$',' read -a DATABASES_SPLIT <<< $DATABASES
                DATABASES_SPLIT=($( uniqify "${DATABASES_SPLIT[*]}" ))
                DATABASE=$(IFS=,; echo "${DATABASES_SPLIT[*]}")
                echo "$PREVSYMBOL	$ENSID	$DATABASE" >> $TMPFILE
            done
        fi
        DATABASES=''
        ENSIDS=()
    fi
    PREVSYMBOL=$(echo ${LINE[0]} | sed -e 's/^ *//' -e 's/ *$//') # arg .. trim!
    ENSIDS+=(${LINE[1]})
    DATABASES+=','${LINE[2]} # adds DBs again if the HGNC_ID == $PREVSYMBOL; gets uniq'ed out later on
done < $FULLLISTSORTED

ENSIDS=( $( uniqify "${ENSIDS[*]}" ) )
for ENSID in ${ENSIDS[@]}; do
    DATABASES+=','$FULLLISTNAME
    IFS=$',' read -a DATABASES_SPLIT <<< $DATABASES
    DATABASES_SPLIT=($( uniqify "${DATABASES_SPLIT[*]}" ))
    DATABASE=$(IFS=,; echo "${DATABASES_SPLIT[*]}")
    echo "$PREVSYMBOL	$ENSID	$DATABASE" >> $TMPFILE
done

mv $TMPFILE $REPODIR/$OUTPUTFILE
echo "Done."

# add the right headers
echo -n "Adding the headers ..."
TMPFILE=`mktemp -p $PROCESSDIR`
echo "HGNC_ID	Ensembl_gene_id	Database" | cat - $REPODIR/$OUTPUTFILE > $TMPFILE && mv $TMPFILE $REPODIR/$OUTPUTFILE
echo "Done."

# cleanup
trap - EXIT
cleanup
