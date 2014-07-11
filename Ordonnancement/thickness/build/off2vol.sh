#!/bin/bash

OFFNAME=$1
RESOLUTION=$2

binvox -e -d $RESOLUTION -t raw $OFFNAME
NAME=${OFFNAME%%.off}  # Delete .off

mv $NAME.raw ${NAME}-${RESOLUTION}.raw
~/Documents/logiciels/DGtalTools/build/converters/raw2vol -x  $RESOLUTION  -y  $RESOLUTION  -z  $RESOLUTION  -i $NAME-$RESOLUTION.raw -o tmp.vol
~/Documents/volFillExterior/build/volFillExterior -i tmp.vol -o $NAME-$RESOLUTION.vol
rm tmp.vol
rm ${NAME}-${RESOLUTION}.raw
