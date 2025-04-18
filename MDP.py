import os

def write_mdp(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(content.strip() + "\n")

def create_mdp_files(work_dir):
    # Energy minimization - first pass
    write_mdp(os.path.join(work_dir, "ions.mdp"), """\
integrator      = steep
emtol           = 1000.0     ; Initial loose minimization
emstep          = 0.01
nsteps          = 5000
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
""")

    # Energy minimization - refined
    write_mdp(os.path.join(work_dir, "minim.mdp"), """\
integrator      = steep
emtol           = 100.0      ; Stricter minimization
emstep          = 0.01
nsteps          = 50000
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
""")

    # Optional: NVT warm-up phase
    write_mdp(os.path.join(work_dir, "nvt_warmup.mdp"), """\
title           = NVT Warm-Up
define          = -DPOSRES
integrator      = md
nsteps          = 10000       ; 20 ps warm-up
dt              = 0.002
tcoupl          = V-rescale
tc_grps  = Protein Non-Protein
tau_t    = 0.1 0.1
ref_t    = 300 300
constraints             = h-bonds
constraint_algorithm    = lincs
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
verlet-buffer-tolerance = 0.0008
pbc = xyz
comm-mode = Linear
comm-grps = System
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000
nstlog      = 1000
gen_vel     = yes
gen_temp    = 100
gen_seed    = -1
""")

    # NVT Equilibration (Stable)
    nvt_mdp = """\
integrator  = md
nsteps      = 50000
dt          = 0.002
tcoupl      = V-rescale
tc_grps  = Protein Non-Protein
tau_t    = 0.1 0.1
ref_t    = 300 300
constraints         = h-bonds
constraint_algorithm= lincs
cutoff-scheme       = Verlet
coulombtype         = PME
rcoulomb            = 1.0
rvdw                = 1.0
verlet-buffer-tolerance = 0.0008
pbc = xyz
comm-mode = Linear
comm-grps = System
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000
nstlog      = 1000
define = -DPOSRES
gen_vel = no
continuation = yes
"""
    write_mdp(os.path.join(work_dir, "nvt.mdp"), nvt_mdp)
    write_mdp(os.path.join(work_dir, "nvt_300.mdp"), nvt_mdp)

    # NPT Equilibration
    npt_mdp = """\
title           = NPT Equilibration
define          = -DPOSRES
integrator      = md
nsteps          = 50000
dt              = 0.002
tcoupl          = V-rescale
tc_grps  = Protein Non-Protein
tau_t    = 0.1 0.1
ref_t    = 300 300
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
ref_p           = 1.0
compressibility = 4.5e-5
refcoord_scaling = com              
constraints             = h-bonds
constraint_algorithm    = lincs
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
verlet-buffer-tolerance = 0.0008
pbc = xyz
comm-mode = Linear
comm-grps = System
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000
nstlog      = 1000
continuation = yes
gen_vel = no
"""
    write_mdp(os.path.join(work_dir, "npt.mdp"), npt_mdp)
    write_mdp(os.path.join(work_dir, "npt_300.mdp"), npt_mdp)

    # Production MD
    md_mdp = """\
integrator              = md
nsteps                  = 250000         ; 500 ps
dt                      = 0.002
tcoupl                  = V-rescale
tc_grps  = Protein Non-Protein
tau_t    = 0.1 0.1
ref_t    = 300 300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
ref_p                   = 1.0
compressibility         = 4.5e-5
constraints             = h-bonds
constraint_algorithm    = lincs
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
verlet-buffer-tolerance = 0.0008
pbc                     = xyz
comm-mode               = Linear
comm-grps               = System
nstxout                 = 1000
nstvout                 = 1000
nstenergy               = 1000
nstlog                  = 1000
continuation            = yes
gen_vel                 = no
"""
    write_mdp(os.path.join(work_dir, "md.mdp"), md_mdp)
    write_mdp(os.path.join(work_dir, "md_300.mdp"), md_mdp)


