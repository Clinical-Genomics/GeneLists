#!/bin/bash

[[ -z $1 ]] && echo "Usage: $0 genelist" && exit
[[ ! -e $1 ]] && >&2 echo "$1 not found" && exit 1

CWD=`pwd`
cd $(dirname `readlink -n $1`)
git describe
cd $CWD
