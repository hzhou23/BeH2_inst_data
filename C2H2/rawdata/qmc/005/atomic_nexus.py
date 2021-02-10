#! /usr/bin/env python3

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf
from nexus import generate_convert4qmc
from nexus import generate_qmcpack,loop,linear,vmc,dmc

settings(
    pseudo_dir    = '/global/homes/h/haihan/Research/C2H2/qmc/001/pseudo',
    results       = '.',
    runs          = 'runs',
    status_only   = 0,
    generate_only = 0,
    sleep         = 10,      # In seconds
    machine       = 'cori', # Use lowercase letters
    account       = 'm2113',   # QMCPACK
    user          = 'haihan'
    )

scf_presub = '''
export HDF5_USE_FILE_LOCKING=FALSE
module unload cray-libsci
module load cray-hdf5
#module load python/3.7-anaconda-2019.10
source activate myenv
module list
'''

system = generate_physical_system(
    structure = 'geom.xyz',
    C          = 4,                    # Zeff=7 for ccECP
    H          = 1,
    net_spin   = 0,   # 2S
    net_charge = 0,
    )

sims = []
# perform Hartree-Fock
scf = generate_pyscf(
    identifier = 'scf',               # log output goes to scf.out
    path       = 'scf',            # directory to run in
    job          = job(serial=True, nodes=4, queue='debug', hours=0.5, constraint='knl', presub=scf_presub),
    #job        = job(serial=True, nodes=1, threads=4, presub=scf_presub),
    template   = 'scf_guess/scf_template.py', # pyscf template file
    system     = system,
    mole       = obj(                 # used to make Mole() inputs
        verbose  = 4,
        #symmetry = 'D2h',
        ),
    save_qmc   = True,                # save wfn data for qmcpack
    )
sims.append(scf)

#### convert orbitals to QMCPACK format
c4q = generate_convert4qmc(
    identifier   = 'c4q',
    path         = 'scf',
    job          = job(nodes=1, queue='debug', hours=0.50, constraint='knl', presub=scf_presub),
    #job        = job(serial=True, nodes=1, threads=4),
    no_jastrow   = True,
    hdf5         = True,              # use hdf5 format
    dependencies = (scf,'orbitals'),
    )
sims.append(c4q)

# collect dependencies relating to orbitals
orbdeps = [(c4q,'particles'), # pyscf changes particle positions
           (c4q,'orbitals' ),]

#######OPTIMIZATION methods
linopt1 = linear(
    energy               = 0.0,
    unreweightedvariance = 1.0,
    reweightedvariance   = 0.0,
    samples              = 2000,
    substeps             = 2,
    steps                = 20,
    blocks               = 100,
    nonlocalpp           = True,
    usedrift             = True,
    minmethod            = 'OneShiftOnly',
    minwalkers           = 1e-4,
    timestep             = 1.5,
#   twistnum             = 0,
    #usebuffer            = True,
    #beta                 = 0.0,
    #exp0                 = -15,
    #bigchange            = 20.0,
    #alloweddifference    = 1e-5,
    #stepsize             = 0.2,
    #stabilizerscale      = 2.0,
    #nstabilizers         = 3
    )
linopt2 = linopt1.copy()
linopt2.minwalkers = 0.10
linopt2.samples = linopt1.samples*2
linopt3 = linopt2.copy()
linopt3.unreweightedvariance = 0.0
linopt3.reweightedvariance = 0.10
linopt3.energy = 0.90
linopt3.minwalkers = 0.30
linopt3.samples = linopt2.samples*2
linopt4 = linopt3.copy()
linopt4.minwalkers = 0.85   # 0.5
linopt4.samples = linopt3.samples*4

#### optimize 12-body Jastrow
optJ12 = generate_qmcpack(
    identifier     = 'optJ12',
    path           = 'optJ12',
    job            = job(nodes=1, queue='debug', hours=0.5, constraint='knl', presub=scf_presub),
    #job            = job(nodes=1, threads=4),
    system         = system,
    J2             = True,         # 2-body B-spline Jastrow
    J1_rcut        = 5.0,          # 4 Bohr cutoff for J1
    J2_rcut        = 5.0,          # 7 Bohr cutoff for J2
    pseudos        = ['C.ccECP.xml', 'H.ccECP.xml'],
    calculations   = [
                loop(max=3, qmc=linopt1),
                loop(max=3, qmc=linopt2),
                loop(max=3, qmc=linopt3),
        ],
    dependencies = orbdeps,
    )
sims.append(optJ12)

##### optimize 3-body Jastrow
optJ123 = generate_qmcpack(
    identifier     = 'optJ123',
    path           = 'optJ123',
    job            = job(nodes=4, queue='regular', hours=1, constraint='knl', presub=scf_presub),
    #job            = job(nodes=1, threads=4),
    system         = system,
    J3             = True,         # 2-body B-spline Jastrow
    J3_rcut        = 5.0,          # 7 Bohr cutoff for J2
    pseudos        = ['C.ccECP.xml', 'H.ccECP.xml'],
    calculations   = [
                loop(max=3, qmc=linopt4),
        ],
    dependencies      = orbdeps+[(optJ12,'jastrow')],
    )
sims.append(optJ123)

NODES=1
### run DMC with QMCPACK
gs_vmc = generate_qmcpack(
    identifier   = 'gs_vmc',
    path         = 'gs_vmc',
    job          = job(nodes=NODES, queue='regular', hours=4.0, constraint='knl', presub=scf_presub),
#    job          = job(nodes=1, threads=8),
    system       = system,
    pseudos        = ['C.ccECP.xml', 'H.ccECP.xml'],
    jastrows     = [],
    calculations   = [ 
        vmc(
            walkers     = int(2048.0/(NODES*8)),   # Per MPI
            #walkers     =   1,
            warmupsteps =  20,
            blocks      = 100,
            steps       =  20,
            substeps    =   2,
            timestep    = 1.0,
            ),
        ],
    dependencies = orbdeps+[(optJ123,'jastrow')],
    )
sims.append(gs_vmc)

NODES=4
### run DMC with QMCPACK
gs_qmc = generate_qmcpack(
    identifier   = 'gs_qmc',
    path         = 'gs_qmc',
    job          = job(nodes=NODES, queue='regular', hours=4.0, constraint='knl', presub=scf_presub),
    #job          = job(nodes=1, threads=4),
    system       = system,
    pseudos        = ['C.ccECP.xml', 'H.ccECP.xml'],
    jastrows     = [],
    calculations   = [ 
        vmc(
            walkers     = int(4096.0/(NODES*8)),   # Per MPI
            #walkers     =   1,
            warmupsteps =  20,
            blocks      = 100,
            steps       =  20,
            substeps    =   2,
            timestep    = 1.0,
            ),
        dmc(targetwalkers = 4096,  # Total walkers
            timestep = 0.02,
            warmupsteps =  int(1.0/0.02),
            blocks = 20,
            steps = int(1.0/0.02),
            nonlocalmoves = 'no',
            checkpoint = 5
            ),
        dmc(targetwalkers = 4096,  # Total walkers
            timestep = 0.01,
            warmupsteps =  int(0.5/0.01),
            blocks = 20,
            steps = int(1.0/0.01),
            nonlocalmoves = 'no',
            checkpoint = 5
            ),
        dmc(targetwalkers = 4096,  # Total walkers
            timestep = 0.005,
            warmupsteps =  int(0.5/0.005),
            blocks = 20,
            steps = int(1.0/0.005),
            nonlocalmoves = 'no',
            checkpoint = 5
            ),

        ],
    dependencies = orbdeps+[(optJ123,'jastrow')],
    )
sims.append(gs_qmc)

run_project(sims)

