#!/bin/sh
##############################################################################
# This is an example of the kind of shell script that would be
# combined with run_job to do a parametric sweep of a set of problems,
# where the parameter to be changed in each session file ($SSN) is
# given by the numeric name of the directory in which job is to be
# run.
#
# $Id$
##############################################################################

SSN=Msten50_50D
ROOTD=..

WD=`pwd` # -- Assume last extension of directory name is new parameter.
NP=`echo $WD | sed "s+.*./++"`

rm -rf * .sequence
sed "s/URED.*=.*/URED=$NP/" < $ROOTD/$SSN  > $SSN

ln -s $ROOTD/$SSN.num .
ln -s $ROOTD/$SSN.rst .
ln -s $ROOTD/$SSN.geom .
ln -s $ROOTD/$SSN.map .

nice /home/hmb/bin/dns_tbcs $SSN > /dev/null
/home/hmb/bin/save $SSN
