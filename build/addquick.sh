#!/bin/bash
#
# $Id$
#
# This is an example of how to use awk to do a simple addition of a
# new field variable, in this case perturbation KE from a 2D avg file.

convert $1 | head \
  | sed 's/ABCpuv /ABCpuvq/' > .tmp

#          123456789012345678901234567

convert $1 | chop -s 11 \
  | awk '{ print $1, $2, $3, $4, $5, $6, 0.5*($1+$3) }' >> .tmp

convert .tmp

rm .tmp

