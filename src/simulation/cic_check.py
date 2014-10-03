#!/opt/Ptools/Python/2.5/bin/python

import pylab as M
import numpy as N
import sys

tol = 0.0001

argv = sys.argv

if len(argv) < 3:
    print 'USAGE: cic_check.py <refName> <inName>'
    sys.exit(-1)

refName = argv[1]
inName = argv[2]

ref = M.load(refName)
refIndx = ref[:,0].astype(int)
refVal = ref[:,4]

a = M.load(inName)
aIndx = a[:,0].astype(int)
aVal = a[:,4]

for i in range(0,aIndx.size):
    indx = aIndx[i]
    val = aVal[i]
    rv = refVal[indx]
    if abs(val-rv) > tol:
        print '%d %f %f' % (indx,rv,val)

