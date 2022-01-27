#!/bin/bash
#
# Print a sequence of numbers embedded in a given number of zeros.
# This can be useful for generating a sequence of numbers that the
# shell (and other things) are guaranteed to consider in ascending
# order.
#
# Usage: index00 <start> <end> <inc> <width>
#
# E.g.
# $ index00 1 3 2 4
# 0001
# 0003
#
# $Id$

echo $1 $2 $3 $4 | awk '
  {
	start = $1 ;
	end   = $2 ;
	inc   = $3 ;
        width = $4
  }
  END {
	 for (i = start; i <= end; i += inc)
	 printf "%0'$4'd\n", i
  } '
