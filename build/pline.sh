#!/bin/bash
#
# pline: generate a uniformly-spaced line of points from one position
# to another.
#
# usage: pline n x1 y1 z1 x2 y2 z2
#
# $Id$

case $# in
	7) ;;
	*) echo 'Usage: pline n x1 y1 z1 x2 y2 z2'; exit 1;;
esac

if test $1 -ge 2
then
    echo $1 $2 $3 $4 $5 $6 $7 | awk '
    {
	N  = $1 ;
	X0 = $2 ;
	Y0 = $3 ;
	Z0 = $4 ;
	DX = $5 - $2 ;
	DY = $6 - $3 ;
	DZ = $7 - $4 
    }
    END {
	for (i = 0; i < N; i++)
	    printf ("%8g %8g %8g\n", X0+i*DX/(N-1),Y0+i*DY/(N-1),Z0+i*DZ/(N-1))

    } '
else
    echo 'number of points must be at least 2'
    exit 1
fi

exit 0
