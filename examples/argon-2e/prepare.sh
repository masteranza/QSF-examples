#!/bin/bash
set -e
# Uses template to create Prometeusz batch files
# example:
# ./prepare.sh ptemplates/prun-re-big.sh -n 2048 --dt 0.12 -f 0.22
org=$1
template=`basename ${org}`
template=`echo ${template} | cut -f 1 -d '.'`
shift;
compressed=`echo $@ | tr -cd "[:alnum:]"`
newname="${template}-${compressed}.sh"
change=$@ #`echo "$@"`
echo "Switching [PARAMS]->[$change] and [COMPRESSED]->[$compressed], saved in ./${newname}"
sed "s/PARAMS/$change/g" ${org} | sed "s/COMPRESSED/$compressed/g" > ${newname}
#  ${newname}  > ${newname}
chmod u+x ${newname}
echo "Done"