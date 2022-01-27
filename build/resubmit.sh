#!/bin/bash
#
# Usage: resubmit session time
#
# Returns 1 if the final time in session.rst exceeds time (i.e. no
# resubmit), else 0.  This can be used in many queuing systems as part
# of a termination test.
#
# If no awk, or gawk, try nawk.
#
# $Id$

head $1.rst | grep Time | grep -v step | \
	awk '{t = $1} \
	END {if (t >= '$2') exit 1; else exit 0 }'

