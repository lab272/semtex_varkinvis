#!/usr/bin/env python
#
# Interactive field file composition utility to combine two field
# files, with field selection.  Alternatively, given a single field
# file, choose (and perhaps rename) selected fields.  Write binary
# file.
#
# 1. compose.py file.fld
# 2. compose.py file1.fld file2.fld
#
# Input files are assumed to be in a machine-compatible binary format.
# In the second case, the field files must conform (nr, ns, nel, nz).
# ----------------------------------------------------------------------------

import numpy as np
import sys
import argparse
import fieldfile
import re

if (len(sys.argv) == 2):
    mode = 1
    file1 = fieldfile.Fieldfile(sys.argv[1], "r")
elif (len(sys.argv) == 3):
    mode = 2
    file1 = fieldfile.Fieldfile(sys.argv[1], "r")
    file2 = fieldfile.Fieldfile(sys.argv[2], "r")
    if ((file1.hdr.geometry.nr  != file2.hdr.geometry.nr)   or
        (file1.hdr.geometry.ns  != file2.hdr.geometry.ns)   or
        (file1.hdr.geometry.nz  != file2.hdr.geometry.nz)   or
        (file1.hdr.geometry.nel != file2.hdr.geometry.nel)) :
        print "The two input files do not conform"
        sys.exit(1)
else:
    print "Usage: compose.py file.fld [anotherfile.fld]"
    sys.exit(1)

if (mode == 1):
    print 'Input file contains these fields: indicate the ones you want:' 
    print file1.hdr.fields
    required_fields = raw_input()
    print required_fields
    wanted = re.findall(r'\S+', required_fields)
    print wanted
    
outfile = raw_input("type in an output file name: ")

ofile = fieldfile.Fieldfile(outfile, "w", file1.hdr)
ofile.hdr.fields = wanted[0]
print ofile.hdr
print ofile.ntot

data_dict = {} 
for i, field in enumerate(ff.fields):
    data_dict[field] = data[i]

    # -- Expect uvcp.
    
    if ff.hdr.fields != 'uvcp':
        sys.exit ("c2w.py: expected fields uvcp, found something else")
    
    # -- Zero everything but c.
    
    data_dict['u'].fill(0.0)
    data_dict['v'].fill(0.0)
    data_dict['p'].fill(0.0)
    
    # -- Write output.

    ff.hdr.fields = 'uvwp'
    ff.hdr.time   = 0.0
    ff.hdr.step   = 0
    ff_out = fieldfile.Fieldfile('/dev/stdout', 'w', ff.hdr)
    ff_out.write(data)
