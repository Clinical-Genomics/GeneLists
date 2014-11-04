#!/bin/bash

if [[ -z $1 ]]; then
    echo "Usage: $0 bitbucket-basedir"
    echo "	bitbucket-basedir: the directory with the bitbicket private gene list repositories"
    exit 1
fi

for dir in `ls $1`; do
   lists=(`ls ${1}/${dir}/*txt`)
   echo ${dir}
   for list in ${dir}/${lists[@]}; do
       list=`basename "$list"`
       list=`echo ${list} | sed 's/.*-\(.*\).txt/\1/'`
       echo "     $list"
   done
done
