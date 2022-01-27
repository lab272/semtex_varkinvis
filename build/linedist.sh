#!/bin/bash
#
# linedist: generate a set of 1D points from an initial to a final location
# according to a geometric growth factor (that may be 1).
#
# usage: linedist x1 x2 npoints growth
#
# $Id$

case $# in
4) ;;
*) echo 'Usage: linedist x1 x2 npoints growth'; exit 1;;
esac

echo $1 $2 $3 $4 | awk '
{
    X1      = $1      ;
    X2      = $2      ;
    NP      = $3 - 1  ;
    GR      = $4      ;
}
END {
len = X2 - X1;
printf ("%g\n", X1);
if ( ((GR-1)^2) < 1e-12 ) {
  for (i = 0; i < NP; i++) {
    x = (i + 1)/NP*len + X1;
    printf ("%8g\n", x) ;
  }
} else {
  lin = len * (1-GR)/(1-GR^NP); 
  del = lin;
  for (i = 0; i < NP; i++) {
    frac = del / len; 
    x    = frac * len + X1 ;
    lin *= GR ;
    del += lin ;
    printf ("%8g\n", x) ;
  }
}
} '

exit 0
