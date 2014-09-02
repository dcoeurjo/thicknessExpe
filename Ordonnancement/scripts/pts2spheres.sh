#!/bin/zsh

## Read a PTS file and extract an approximation of the MA from power crust

PTSNAME=$1
ARITHM=$2
CRUSTPATH=../crust/


cp $PTSNAME tmp.pts

$CRUSTPATH/powercrust -m $ARITHM -i tmp.pts
echo "$CRUSTPATH/powercrust -m $ARITHM -i $PTSNAME"

mv inpball $1:r-crust-spheres.txt
mv axis.off $1:r-crust-spheres.off

rm pc.off axis.off  inpole outpole poleinfo

echo "Number of balls= " `wc  -l $1:r-crust-spheres.txt`

