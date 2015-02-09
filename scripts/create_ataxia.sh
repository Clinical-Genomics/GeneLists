#!/bin/bash

# Creates one master list out of several Ataxia gene lists.
#
# usage:
#    create_ataxia.sh path-to-cust-repo [outputfile [database name]]
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
OUTPUTFILE=${2-cust002-ATX.txt}
FULLLISTNAME=${3-'FullList'} # combine multiple gene lists into one big list with this name
GENOMEBUILD=${4-'GRCh37.p13'} # default column

# gene lists we need to process
declare -A FILES=( [neutropenia_DN.txt]=SCN [SpasticParaplegia_genelist_MP_DN.txt]=AD-HSP,SPG,Ataxi [Ataxia_AutRec_MP_DN.txt]=Ataxi [ataxia_list_MP.txt]=Ataxi [other_dominant_ataxia_genelist_MP_DN.txt]=Ataxi [SCA_genelist_MP_DN.txt]=Ataxi [spastic_paraplegia_related_140911_MP_DN.txt]=SPG,Ataxi)

echo -n "Checking if files exist ..."
for f in ${!FILES[@]}; do
    [[ ! -e $REPODIR/original/ClinGeneLists/src/$f ]] && >&2 echo "$f not found!" && exit
done
echo "Done."

# download the genelists - CD to $REPODIR
echo -n "Updating the git repo ..."
cd $REPODIR
git submodule foreach git pull
echo "Done."

# copy them to easily work with them - CD to $REPODIR/original/ClinGeneLists/src:
echo -n "Copying lists ..."
PROCESSDIR=`mktemp -d`
cd $REPODIR/original/ClinGeneLists/src
cp ${!FILES[@]} $PROCESSDIR
echo "Lists are at '$PROCESSDIR'"

# cd to $PROCESSDIR - make sure the rest of the script runs on these files
cd $PROCESSDIR

# add the right Database name
echo "Adding the right database ..."
for f in `ls $PROCESSDIR`; do
    TMPFILE=`mktemp`
    FILENAME=`basename $f`
    DATABASE=${FILES[$FILENAME]}
    echo "	$FILENAME: $DATABASE"
    while read line; do
        echo "$line	$DATABASE" >> $TMPFILE
    done < $f
    mv $TMPFILE $f
done
echo "Done."

# copy them all into one gene list file
echo -n "Copying to the one list ..."
cat ${!FILES[@]} | sort > $REPODIR/$OUTPUTFILE
echo "Done."

# merge duplicated entries
echo -n "Merging duplicated entries ..."
PREVSYMBOL=''
DATABASES=()
TMPFILE=`mktemp`
while read LINE; do
    IFS=$'\t' read -a LINE <<< "$LINE"
    if [[ "${LINE[0]}" == "$PREVSYMBOL" ]]; then
        IFS=',' read -a DBS <<< ${LINE[1]}
        DATABASES=( "${DATABASES[*]}" "${DBS[*]}" )
    else
        if [[ -n $PREVSYMBOL ]]; then
            DATABASES+=($FULLLISTNAME)
            DATABASES=($( uniqify "${DATABASES[*]}" ))
            DATABASE=$(IFS=,; echo "${DATABASES[*]}")
            echo "$PREVSYMBOL	$GENOMEBUILD	$DATABASE" >> $TMPFILE
        fi
        IFS=',' read -a DATABASES <<< ${LINE[1]}
    fi
    PREVSYMBOL=${LINE[0]}
done < $REPODIR/$OUTPUTFILE
DATABASES+=($FULLLISTNAME)
DATABASES=($( uniqify "${DATABASES[*]}" ))
DATABASE=$(IFS=,; echo "${DATABASES[*]}")
echo "$PREVSYMBOL	$GENOMEBUILD	$DATABASE" >> $TMPFILE

mv $TMPFILE $REPODIR/$OUTPUTFILE
echo "Done."

# add the right headers
echo -n "Adding the headers ..."
TMPFILE=`mktemp`
echo "HGNC_ID	Clinical_db_genome_build	Clinical_db_gene_annotation" | cat - $REPODIR/$OUTPUTFILE > $TMPFILE && mv $TMPFILE $REPODIR/$OUTPUTFILE
echo "Done."

# cleanup
trap - EXIT
cleanup
