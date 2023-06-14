#!/usr/bin/env python
#
# wavy.py -- generate a semtex session file with elements which are
# sinusoidally wavy in x but where wave amplitude varies linearly with
# y, much like the wavypipe meshes.  We use these semtex utilities:
# 
#   rectmesh
#   mapmesh
#
# The initial domain size is in [0,1] x [0,1], and the (x, y)
# locations are hardcoded here rather than bothering to read them in
# from a rectmesh input file.
# 
# Usage: wavy.py <beta_x> <eps_y> session
#   <beta_x> ... real number that maps the x locations
#   <eps_y>  ... real number that supplies that maximum wave amplitude
#   session  ... name of session file that will be created
# 
# NB: the SURFACES, BCS, GROUPS sections will likely need hand editing.
# ----------------------------------------------------------------------------

import sys
import os
import numpy as np
import subprocess as sp

if (len(sys.argv) == 4):
    beta_x  = float(sys.argv[1])
    eps_y   = float(sys.argv[2])
    session = sys.argv[3]
else:
    print ("Usage: wavy.py <beta_x> <eps_y> session")
    sys.exit(1)

#print beta_x
#print eps_y
#print (session)

# -- Here we create information equivalent to that in a rectmesh input file.

#x = np.linspace(0.0, 1.0, num=4)
#x = x * 2.0 * np.pi / beta_x
#y = np.array([0, 0.5, 0.7, 1.0])

x = np.linspace(0.0, 1.0, num=11)
x = x * 2.0 * np.pi / beta_x

y = np.array([0, 0.075, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.84, 0.9, 0.945, 0.975, 0.993, 1.0])

# -- And generate the file.


rfile = open (session+'.rect', 'w')

i = 0
while i < x.size:
    rfile.write ('%24.16f' % x[i] + '\n')
    i = i+1

rfile.write ('\n')

i = 0
while i < y.size:
    rfile.write ('%24.16f' % y[i] + '\n')
    i = i+1

rfile.close()

# -- Run rectmesh to generate session file.

sfile = open (session, 'w')

sp.run (['rectmesh', session+'.rect'], stdout=sfile)

sfile.close()

# -- Run mapmesh to make NODES in session file wavy in x.

sfile = open (session+'_wavy', 'w')

mapcommand = '-y ' + 'y*(1+' + str(eps_y) + '*sin(' + str(beta_x) + '*x))'

#print(mapcommand)

sp.run (['mapmesh', mapcommand, session], stdout=sfile)

sfile.close()

# -- Make the input files for the CURVES section of final session file.
#    Nx points in x (previously we have used 960).

Nx = 960
for i in range(1, len(y)):
    sfile = open ('wave' + str(i) + '.geo', 'w')
    for j in range(Nx):
        xl = (x[len(x)-1]-x[0]) * j/(Nx-1)
        yl = y[i]*(1.0+eps_y*np.sin(xl*beta_x))
        sfile.write (str(xl) + ' ' + str(yl) + '\n')
    sfile.close ()

# -- Make the CURVES section in a separate file.

sfile = open('curves.txt', 'w')
sfile.write ('\n<CURVES NUMBER=' + str((len(x)-1)*(2*(len(y)-2)+1)) + '>\n')

k = 1

for i in range(1,len(y)):
    for j in range(1,len(x)):
#        print ('i='+str(i)+', j='+str(j))
        sfile.write ('%5d' % k + '%5d' % ((i-1)*(len(x)-1)+j) +
                     '  3 <SPLINE> wave' + str(i)+'.geo </SPLINE>' + '\n')
        k = k + 1
    sfile.write ('\n')

for i in range(1,len(y)-1):
    for j in range(1,len(x)):
#        print ('i='+str(i)+', j='+str(j))
        sfile.write ('%5d' % k + '%5d' % (i*(len(x)-1)+j) +
                     '  1 <SPLINE> wave' + str(i)+'.geo </SPLINE>' + '\n')
        k = k + 1
    sfile.write ('\n')

sfile.write ('</CURVES>\n')
sfile.close ()

# -- Cat the wavy session file and the CURVES section onto the original session.

cmd = 'cat ' + session + '_wavy' + ' curves.txt > ' + session
failure = os.system (cmd)
cmd = 'rm curves.txt'
failure = os.system (cmd)
