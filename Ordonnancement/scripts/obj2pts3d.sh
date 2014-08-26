#!/bin/zsh

## Read a OBJ file and extract the vertices to a PTS file (x y z)

OBJNAME=$1
NBVERT=`head -n 100 $1| grep "Vertices:" | awk '{print $3}'`
echo $NBVERT

#echo $NBVERT >! $1:r.pts

cat $1 | grep "v " | awk '{ print $2 " " $3 " " $4  }' >! $1:r.pts


