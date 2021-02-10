#!/usr/bin/env python3

import sys,os
from pyscf import scf,gto,dft
import numpy as np
from urllib.request import urlretrieve

#~~~~~~~ Input ~~~~~~~~~~~~~~~~~~~~~~~

#Set the current working directory
cwd = os.getcwd()
pplib = "http://pseudopotentiallibrary.org/recipes"
pptype = "ccECP"

#Obtain basis and ECP files from pseudopotentiallibrary.org
atom1 = "C"
bastype1 = "aug-cc-pVQZ"
basfile1 = "{0}.{1}.nwchem".format(atom1,bastype1)
ecpfile1 = "{0}.{1}.nwchem".format(atom1,pptype)
xmlfile1 = "{0}.{1}.xml".format(atom1,pptype) #grab qmcpack xml file for later
urlretrieve("{0}/{1}/{2}/{3}".format(pplib,atom1,pptype,basfile1),
        filename=basfile1)
urlretrieve("{0}/{1}/{2}/{3}".format(pplib,atom1,pptype,ecpfile1),
        filename=ecpfile1)
urlretrieve("{0}/{1}/{2}/{3}".format(pplib,atom1,pptype,xmlfile1),
        filename=xmlfile1)

atom2 = "H"
bastype2 = "aug-cc-pVQZ"
basfile2 = "{0}.{1}.nwchem".format(atom2,bastype2)
ecpfile2 = "{0}.{1}.nwchem".format(atom2,pptype)
xmlfile2 = "{0}.{1}.xml".format(atom2,pptype) #grab qmcpack xml file for later
urlretrieve("{0}/{1}/{2}/{3}".format(pplib,atom2,pptype,basfile2),
        filename=basfile2)
urlretrieve("{0}/{1}/{2}/{3}".format(pplib,atom2,pptype,ecpfile2),
        filename=ecpfile2)
urlretrieve("{0}/{1}/{2}/{3}".format(pplib,atom2,pptype,xmlfile2),
        filename=xmlfile2)

#~~~~ Build the molecule ~~~~

# Nexus expands this with Mole info
$system

with open(os.path.join(cwd,basfile1)) as f:
    bas1 = f.read()
with open(os.path.join(cwd,basfile2)) as f:
    bas2 = f.read()
with open(os.path.join(cwd,ecpfile1)) as f:
    ecp1 = f.read()
with open(os.path.join(cwd,ecpfile2)) as f:
    ecp2 = f.read()
mol.basis = {atom1: gto.basis.parse(bas1), atom2: gto.basis.parse(bas2)} 
mol.ecp = {atom1: gto.basis.parse_ecp(ecp1), atom2: gto.basis.parse_ecp(ecp2)}
mol.build()

#~~~Run HF on molecule~~~~
mf = scf.RKS(mol)
mf.xc = 'hf'
#mf.irrep_nelec = {
#'Ag' : (4,2),   # s    
#'B3u': (1,1),   # x    1
#'B1u': (1,1),   # z    0
#'B2u': (1,1),   # y   -1
#'B2g': (1,1),   # xz   1
#'B3g': (1,0),   # yz  -1
#'B1g': (1,0),   # xy  -2
#'Au' : (0,0)    # xyz  
#}
mf.verbose=4
mf.max_cycle=100
#dm=mf.init_guess_by_chkfile('../../scf_guess/guess/scf.chkfile')
mf.chkfile='scf.chkfile'
#en=mf.kernel(dm)
en=mf.kernel()

