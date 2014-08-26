#!/bin/zsh

## Read a OFF, digitize it and extract the RDMA


OFFNAME=$1
off2vol.sh $OFFNAME $2
VOLNAME=$OFFNAME:r-$2.vol
THICKNESSPATH=../build/thickness

echo $OFFNAME
echo $VOLNAME


sedt -i $VOLNAME -o tmp.longvol
rdma -i tmp.longvol -o $VOLNAME:r-rdma
rm tmp.longvol

$THICKNESSPATH/longvol2sphere -i $VOLNAME:r-rdma.longvol >! $VOLNAME:r-spheres.txt
rm $VOLNAME:r-rdma.longvol

# Output : Spheres file : $VOLNAME:r-spheres.txt
