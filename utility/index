#!/bin/bash
#
# Print a sequence of numbers from "start" to "end" in steps "inc",
# using awk.  This functionality is usually present with the shell
# utility seq, but that does not always exist.
#
# Usage: index <start> <end> <inc>
#
# $Id$

echo $1 $2 $3 | awk '
  {
	start = $1 ;
	end   = $2 ;
	inc   = $3 
  }
  END {
	 for (i = start; i <= end; i += inc)
	 print i
  } '
