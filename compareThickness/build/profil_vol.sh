#!/bin/bash

OFFNAME=$1
~/Documents/off2vol.sh $OFFNAME $2
VOLNAME=${OFFNAME%%.off}-$2.vol
FILENAME=${OFFNAME%%.off}-$2-data.txt

voladdborder $VOLNAME tmp.vol
rm $VOLNAME
mv tmp.vol $VOLNAME

thicknessComputationHildebrand -i $VOLNAME -o tmp.longvol > $FILENAME
rm tmp.longvol 
rm $VOLNAME

SDFFILE=${OFFNAME%%.off}-sdf.txt
FACES=$(wc -l $SDFFILE | cut -d" " -f1)
FACES=$(($FACES-2))  # because of the header of the .off file
MAX=$(tail -1 $FILENAME)
VOXELS=$(wc -l $FILENAME | cut -d" " -f1)

./compare_resol ${OFFNAME%%.off}-$2 $FACES $MAX $VOXELS 
rm $FILENAME
