import os
import subprocess
import requests
import argparse

# Constants
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PDB_DIR = os.path.join(SCRIPT_DIR, "pdbs")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def run_command(cmd, input=None, cwd=RESULTS_DIR):
    print(f"⚙️ Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, input=input, text=True, capture_output=True, cwd=cwd)
    if result.returncode != 0:
        print(f"❌ Error:\n{result.stderr}")
        raise RuntimeError("Command failed")
    return result.stdout

def download_pdb(pdb_id, outdir=PDB_DIR):
    ensure_dir(outdir)
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    out_path = os.path.join(outdir, f"{pdb_id}.pdb")
    if not os.path.exists(out_path):
        print(f"📥 Downloading {pdb_id}...")
        r = requests.get(url)
        if r.ok:
            with open(out_path, "w") as f:
                f.write(r.text)
        else:
            raise ValueError(f"Could not download {pdb_id}")
    return out_path

def write_mdp(filename, content):
    full_path = os.path.join(RESULTS_DIR, filename)
    os.makedirs(os.path.dirname(full_path), exist_ok=True)
    if os.path.exists(full_path):
        os.remove(full_path)
    with open(full_path, "w") as f:
        f.write(content)

def prepare_gromacs_system(pdb_file):
    # Example logic — assumes pdb_file is used to generate outputs in RESULTS_DIR
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    output_gro = os.path.join(RESULTS_DIR, f"{base_name}.gro")
    output_top = os.path.join(RESULTS_DIR, f"{base_name}.top")

    # Example command (replace with real GROMACS ones)
    run_command(["echo", f"Simulating GROMACS prep for {pdb_file}"])
    # Simulate output creation
    with open(output_gro, "w") as f:
        f.write(f"; GRO file for {base_name}\n")
    with open(output_top, "w") as f:
        f.write(f"; TOP file for {base_name}\n")


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

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Water_and_ions
tau_t           = 0.5 0.5
ref_t           = 100 100

; Constraints
constraints             = h-bonds
constraint_algorithm    = lincs

; Cutoffs
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0

; Verlet settings
verlet-buffer-tolerance = 0.0008

; Periodic Boundary Conditions
pbc = xyz

; Center of mass motion removal
comm-mode = Linear
comm-grps = System

; Output
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000
nstlog      = 1000

; Velocity generation
gen_vel     = yes
gen_temp    = 100
gen_seed    = -1
""")

    # Stable NVT equilibration (updated as per your exact spec)
    write_mdp(os.path.join(work_dir, "nvt.mdp"), """\
integrator  = md
nsteps      = 50000
dt          = 0.002

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = Protein Water_and_ions
tau_t       = 0.5 0.5
ref_t       = 300 300

; Constraints
constraints         = h-bonds
constraint_algorithm= lincs

; Cutoffs
cutoff-scheme       = Verlet
coulombtype         = PME
rcoulomb            = 1.0
rvdw                = 1.0

; Verlet settings
verlet-buffer-tolerance = 0.0008

; Periodic Boundary Conditions
pbc = xyz

; Center of mass motion removal
comm-mode = Linear
comm-grps = System

; Output (minimal for NVT)
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000
nstlog      = 1000

; Position restraints
define = -DPOSRES

; Velocity generation (off by default, only use in warm-up)
gen_vel = no
continuation = yes
""")

    # NPT Equilibration (slightly tweaked for consistency)
    write_mdp(os.path.join(work_dir, "npt.mdp"), """\
title           = NPT Equilibration
define          = -DPOSRES
integrator      = md
nsteps          = 50000
dt              = 0.002

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Water_and_ions
tau_t           = 0.5 0.5
ref_t           = 300 300

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
ref_p           = 1.0
compressibility = 4.5e-5

refcoord_scaling = com              

; Constraints
constraints             = h-bonds
constraint_algorithm    = lincs

; Cutoffs
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0

; Verlet settings
verlet-buffer-tolerance = 0.0008

; Periodic Boundary Conditions
pbc = xyz

; Center of mass motion removal
comm-mode = Linear
comm-grps = System

; Output
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000
nstlog      = 1000

continuation = yes
gen_vel = no
""")

    # Production MD
    write_mdp(os.path.join(work_dir, "md.mdp"), """\
; Production MD parameters
integrator              = md
nsteps                  = 250000         ; 500 ps
dt                      = 0.002

; Temperature coupling
tcoupl                  = V-rescale
tc-grps                 = Protein Water_and_ions
tau_t                   = 0.5 0.5
ref_t                   = 300 300

; Pressure coupling
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
ref_p                   = 1.0
compressibility         = 4.5e-5

; Constraints
constraints             = h-bonds
constraint_algorithm    = lincs

; Cutoffs
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0

; Verlet buffer
verlet-buffer-tolerance = 0.0008

; Periodic boundary conditions
pbc                     = xyz

; Remove center of mass motion
comm-mode               = Linear
comm-grps               = System

; Output control
nstxout                 = 1000
nstvout                 = 1000
nstenergy               = 1000
nstlog                  = 1000

; Simulation continuation
continuation            = yes
gen_vel                 = no
""")

def prepare_gromacs_system(pdb_path, ff="amber99sb-ildn", water="tip3p"):
    base = os.path.splitext(os.path.basename(pdb_path))[0]
    work_dir = os.path.dirname(os.path.abspath(pdb_path))
    create_mdp_files(work_dir)

    processed = os.path.join(work_dir, f"{base}_processed.gro")
    boxed = os.path.join(work_dir, "boxed.gro")
    solvated = os.path.join(work_dir, "solvated.gro")
    solv_ions = os.path.join(work_dir, "solv_ions.gro")
    topol = os.path.join(work_dir, "topol.top")

    # 1. pdb2gmx
    run_command(["gmx", "pdb2gmx", "-f", pdb_path, "-o", processed, "-p", topol,
                 "-i", os.path.join(work_dir, "posre.itp"), "-ff", ff, "-water", water, "-ignh"])

    # 2. editconf
    run_command(["gmx", "editconf", "-f", processed, "-o", boxed, "-c", "-d", "1.0", "-bt", "cubic"])

    # 3. solvate
    run_command(["gmx", "solvate", "-cp", boxed, "-cs", "spc216.gro", "-p", topol, "-o", solvated])

    # 4. Add ions
    ions_tpr = os.path.join(work_dir, "ions.tpr")
    run_command(["gmx", "grompp", "-f", os.path.join(work_dir, "ions.mdp"),
                 "-c", solvated, "-p", topol, "-o", ions_tpr, "-maxwarn", "2"])
    run_command(["gmx", "genion", "-s", ions_tpr, "-o", solv_ions, "-p", topol,
                 "-pname", "NA", "-nname", "CL", "-neutral"], input="13\n")

    # 5. Energy minimization
    em_tpr = os.path.join(work_dir, "em.tpr")
    run_command(["gmx", "grompp", "-f", os.path.join(work_dir, "minim.mdp"),
                 "-c", solv_ions, "-p", topol, "-o", em_tpr])
    run_command(["gmx", "mdrun", "-deffnm", os.path.join(work_dir, "em")])

    # 6. NVT warm-up (gentle heating to 100 K)
    nvt_warmup_tpr = os.path.join(work_dir, "nvt_warmup.tpr")
    run_command(["gmx", "grompp", "-f", os.path.join(work_dir, "nvt_warmup.mdp"),
                 "-c", os.path.join(work_dir, "em.gro"), "-r", os.path.join(work_dir, "em.gro"),
                 "-p", topol, "-o", nvt_warmup_tpr, "-maxwarn", "2"])
    run_command(["gmx", "mdrun", "-deffnm", os.path.join(work_dir, "nvt_warmup")])

    # 7. NVT equilibration (at 300 K)
    nvt_tpr = os.path.join(work_dir, "nvt.tpr")
    run_command(["gmx", "grompp", "-f", os.path.join(work_dir, "nvt.mdp"),
                 "-c", os.path.join(work_dir, "nvt_warmup.gro"), "-r", os.path.join(work_dir, "nvt_warmup.gro"),
                 "-t", os.path.join(work_dir, "nvt_warmup.cpt"), "-p", topol, "-o", nvt_tpr])
    run_command(["gmx", "mdrun", "-deffnm", os.path.join(work_dir, "nvt")])

    # 8. NPT equilibration
    npt_tpr = os.path.join(work_dir, "npt.tpr")
    run_command(["gmx", "grompp", "-f", os.path.join(work_dir, "npt.mdp"),
                 "-c", os.path.join(work_dir, "nvt.gro"), "-r", os.path.join(work_dir, "nvt.gro"),
                 "-t", os.path.join(work_dir, "nvt.cpt"), "-p", topol, "-o", npt_tpr])
    run_command(["gmx", "mdrun", "-deffnm", os.path.join(work_dir, "npt")])

    # 9. Production run
    md_tpr = os.path.join(work_dir, "md.tpr")
    run_command(["gmx", "grompp", "-f", os.path.join(work_dir, "md.mdp"),
                 "-c", os.path.join(work_dir, "npt.gro"), "-t", os.path.join(work_dir, "npt.cpt"),
                 "-p", topol, "-o", md_tpr])
    run_command(["gmx", "mdrun", "-deffnm", os.path.join(work_dir, "md"), "-gpu_id", "0"])

    print(f"🎉 MD simulation complete for {base}")

def main(pdb_ids):
    ensure_dir(PDB_DIR)
    ensure_dir(RESULTS_DIR)
    for pdb_id in pdb_ids:
        try:
            pdb_file = download_pdb(pdb_id)
            prepare_gromacs_system(pdb_file)
        except Exception as e:
            print(f"❌ Failed for {pdb_id}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GROMACS Full Pipeline Runner")
    parser.add_argument("--pdbs", nargs="+", required=True, help="List of PDB IDs (e.g. 1HZH 4PE5)")
    args = parser.parse_args()
    main(args.pdbs)