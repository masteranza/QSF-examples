#!/bin/bash
set -e
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