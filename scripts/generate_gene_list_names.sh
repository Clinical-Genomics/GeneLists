#!/bin/bash

if [[ -z $1 ]]; then
    echo "Usage: $0 bitbucket-basedir"
    echo "	bitbucket-basedir: the directory with the bitbicket private gene list repositories"
    exit 1
fi

for dir in `ls $1`; do
    echo ${dir}
    if [[ -e ${1}/${dir}/LISTS ]]; then
        while read LINE; do
            echo "	$LINE"
        done < ${1}/${dir}/LISTS 
    else
        lists=(`ls ${1}/${dir}/*txt`)
        for list in ${dir}/${lists[@]}; do
            list=`basename "$list"`
            list=`echo ${list} | sed 's/.*-\(.*\).txt/\1/'`
            echo "	$list"
        done
    fi
done
