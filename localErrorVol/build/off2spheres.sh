#!/bin/bash

OFFNAME=$1
~/Bureau/off2vol.sh $OFFNAME $2
VOLNAME=${OFFNAME%%.off}-$2.vol

sedt -i $VOLNAME -o tmp.longvol
rdma -i tmp.longvol -o ${VOLNAME%%.vol}
rm tmp.longvol

~/Bureau/longvol2sphere ${VOLNAME%%.vol}.longvol
rm ${VOLNAME%%.vol}.longvol

# Output : Spheres file : ${VOLNAME%%.vol}-spheres.txt
