***,Calculation of O
memory,1,g

gthresh,oneint=1.0E-13
gthresh,twoint=1.0E-13

!set,dkroll=1,dkho=99,dkhp=2

geometry={
   4
   C2H2
   C  0.0000  0.0000  0.6013
   C  0.0000  0.0000 -0.6013
   H  0.0000  0.0000  1.6644
   H  0.0000  0.0000 -1.6644
}

basis={
include,/global/cscratch1/sd/aannabe/repos/pseudopotentiallibrary/recipes/C/ccECP/C.ccECP.molpro
include,/global/cscratch1/sd/aannabe/repos/pseudopotentiallibrary/recipes/H/ccECP/H.ccECP.molpro

include,/global/cscratch1/sd/aannabe/repos/pseudopotentiallibrary/recipes/C/ccECP/C.aug-cc-pV5Z.molpro
include,/global/cscratch1/sd/aannabe/repos/pseudopotentiallibrary/recipes/H/ccECP/H.aug-cc-pV5Z.molpro
}

{rhf,nitord=30;
 maxit,100;
 wf,10,1,0
 print,orbitals=1
}

scf(i)=energy
_CC_NORM_MAX=2.0
{uccsd(t),shifts=0.1,shiftp=0.1,thrdis=1.0;diis,1,1,15,1;maxit,100;core}
ccsd(i)=energy

table,scf,ccsd
save, 5.csv, new

