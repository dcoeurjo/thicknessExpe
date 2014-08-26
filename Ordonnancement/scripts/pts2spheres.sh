#!/bin/zsh

## Read a PTS file and extract an approximation of the MA from power crust

PTSNAME=$1
CRUSTPATH=../crust/

$CRUSTPATH/powercrust -m 100000 -i $1

mv inpball $1:r-crust-spheres.txt

rm   pole inpball pc.off axis.off axisface.off inpole outpole poleinfo

echo "Number of balls= " `wc  -l $1:r-crust-spheres.txt`


