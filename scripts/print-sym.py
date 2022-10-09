#!/usr/bin/python
#
import sys
import os
import openbabel as ob
from pybel import *

if len(sys.argv)==2:
    input_filename = sys.argv[1]
else:
    print("missing command line arg: filename")
    sys.exit(1)

for molecule in readfile("pdb",input_filename):

    canonlabels = ob.vectorUnsignedInt()
    symclasses = ob.vectorUnsignedInt()
    accountedfor = ob.vectorUnsignedInt(len(molecule.atoms)+1)

    gs=ob.OBGraphSym(molecule.OBMol)
    gs.GetSymmetry(symclasses)
    ob.CanonicalLabels(molecule.OBMol, symclasses, canonlabels)

    i=0
    for a in molecule.atoms:
        print symclasses[i]
        i+=1
