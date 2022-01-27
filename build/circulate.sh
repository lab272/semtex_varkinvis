#!/bin/sh
#
# Circulate a base flow file's sequence by $1, write to stdout.
#
# usage: circulate n file
#
# $Id$
###############################################################################

if [ $# -ne 2 ]
then
  echo usage: circulate n file
  exit 1
fi

C=$1
BASE=$2

N=`convert $BASE | grep Session | wc -l`
if [ $C -gt $N ]
then
  echo increment $C exceeds number of slices $N
  exit 1
fi

for i in `seq $N`
do
  j=`echo \($i + $C\) % $N | bc`
  if [ $j -eq 0 ]
  then
    j=$N
  fi
#  echo $i $j

  convert -a -n $i $BASE | chop -n 10 > $BASE.head
  convert -a -n $j $BASE | chop -s 11 > $BASE.data
  cat $BASE.head $BASE.data | convert  
done

rm $BASE.head $BASE.data
