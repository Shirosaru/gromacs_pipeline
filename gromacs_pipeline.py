import os
import subprocess
import requests
import argparse
import Immunebuilder
import MDP
import shutil

# 📁 Directory Constants (relative to this script's location)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PDB_DIR = os.path.join(SCRIPT_DIR, "pdbs")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
IMMUNE_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "immune_output")


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def run_command(cmd, input=None, cwd=RESULTS_DIR):
    print(f"⚙️ Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, input=input, text=True, capture_output=True, cwd=cwd)
    if result.returncode != 0:
        print(f"❌ Error:\n{result.stderr}")
        raise RuntimeError("Command failed")
    return result.stdout

def clean_pdb(input_pdb, output_pdb, keep_water=False):
    """
    Cleans a PDB file by removing HETATM entries and keeping only protein ATOM lines.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pdb (str): Path to write the cleaned PDB file.
        keep_water (bool): If True, keeps water molecules (residue name HOH).
    """
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            record_type = line[0:6].strip().upper()
            res_name = line[17:20].strip().upper()

            if record_type == "ATOM":
                f_out.write(line)
            elif keep_water and record_type == "HETATM" and res_name == "HOH":
                f_out.write(line)

    print(f"✅ Cleaned PDB written to {output_pdb}")


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


def edit_mdp(original_mdp, new_mdp, ref_t=None, gen_temp=None, nsteps=None):
    with open(original_mdp, 'r') as f:
        lines = f.readlines()

    with open(new_mdp, 'w') as f:
        for line in lines:
            if ref_t and "ref_t" in line:
                f.write(f"ref_t = {ref_t[0]} {ref_t[1]}\n")
            elif gen_temp and "gen_temp" in line:
                f.write(f"gen_temp = {gen_temp}\n")
            elif nsteps and "nsteps" in line:
                f.write(f"nsteps = {nsteps[0]}\n")
            else:
                f.write(line)

def prepare_gromacs_system(pdb_path, temps="300", time=100, ff="amber99sb-ildn", water="tip3p"):
    avail_temps = ['300', '310', '350', '373', '400']
    temps = [str(t.strip()) for t in temps.split(',')]

    base = os.path.splitext(os.path.basename(pdb_path))[0]
    work_dir = os.path.join(RESULTS_DIR, base)
    ensure_dir(work_dir)

    # Copy PDB to the working directory to keep it self-contained
    local_pdb_path = os.path.join(work_dir, os.path.basename(pdb_path))
    shutil.copy(pdb_path, local_pdb_path)

    MDP.create_mdp_files(work_dir)

    processed = os.path.join(work_dir, f"{base}_processed.gro")
    boxed = os.path.join(work_dir, "boxed.gro")
    solvated = os.path.join(work_dir, "solvated.gro")
    solv_ions = os.path.join(work_dir, "solv_ions.gro")
    topol = os.path.join(work_dir, "topol.top")

    # Preprocessing: pdb2gmx -> solvation + ions
    run_command(["gmx", "pdb2gmx", "-f", pdb_path, "-o", processed, "-p", topol,
                 "-i", os.path.join(work_dir, "posre.itp"), "-ff", ff, "-water", water, "-ignh"])
    run_command(["gmx", "editconf", "-f", processed, "-o", boxed, "-c", "-d", "1.0", "-bt", "cubic"])
    run_command(["gmx", "solvate", "-cp", boxed, "-cs", "spc216.gro", "-p", topol, "-o", solvated])
    ions_tpr = os.path.join(work_dir, "ions.tpr")
    run_command(["gmx", "grompp", "-f", os.path.join(work_dir, "ions.mdp"),
                 "-c", solvated, "-p", topol, "-o", ions_tpr, "-maxwarn", "2"])
    run_command(["gmx", "genion", "-s", ions_tpr, "-o", solv_ions, "-p", topol,
                 "-pname", "NA", "-nname", "CL", "-neutral"], input="13\n")

    # Energy minimization
    em_tpr = os.path.join(work_dir, "em.tpr")
    run_command(["gmx", "grompp", "-f", os.path.join(work_dir, "minim.mdp"),
                 "-c", solv_ions, "-p", topol, "-o", em_tpr])
    run_command(["gmx", "mdrun", "-deffnm", os.path.join(work_dir, "em")])

    # Multi-temperature loop
    for temp in temps:
        print(f"\n🔥 Starting simulation at {temp} K")

        # Set filenames
        nvt_mdp = os.path.join(work_dir, f"nvt_{temp}.mdp")
        npt_mdp = os.path.join(work_dir, f"npt_{temp}.mdp")
        md_mdp = os.path.join(work_dir, f"md_{temp}.mdp")

        # Use default MDPs or modify from 300K
        if temp not in avail_temps:
            edit_mdp(os.path.join(work_dir, "nvt_300.mdp"), new_mdp=nvt_mdp, ref_t=[temp, temp], gen_temp=temp)
            edit_mdp(os.path.join(work_dir, "npt_300.mdp"), new_mdp=npt_mdp, ref_t=[temp, temp])
            edit_mdp(os.path.join(work_dir, "md_300.mdp"), new_mdp=md_mdp, ref_t=[temp, temp])
        else:
            # If using standard ones, assume already in place
            pass

        # Modify MD time if needed
        if time != 100:
            steps = int(time * 1000 * 1000 / 2)  # assuming 0.002 ps timestep
            custom_md_mdp = os.path.join(work_dir, f"md_{temp}_{time}.mdp")
            edit_mdp(md_mdp, new_mdp=custom_md_mdp, nsteps=[steps])
            md_mdp = custom_md_mdp

        # NVT
        nvt_tpr = os.path.join(work_dir, f"nvt_{temp}.tpr")
        run_command(["gmx", "grompp", "-f", nvt_mdp, "-c", os.path.join(work_dir, "em.gro"),
                     "-r", os.path.join(work_dir, "em.gro"), "-p", topol, "-o", nvt_tpr])
        run_command(["gmx", "mdrun", "-deffnm", os.path.join(work_dir, f"nvt_{temp}")])

        # NPT
        npt_tpr = os.path.join(work_dir, f"npt_{temp}.tpr")
        run_command(["gmx", "grompp", "-f", npt_mdp, "-c", os.path.join(work_dir, f"nvt_{temp}.gro"),
                     "-t", os.path.join(work_dir, f"nvt_{temp}.cpt"), "-r", os.path.join(work_dir, f"nvt_{temp}.gro"),
                     "-p", topol, "-o", npt_tpr, "-maxwarn", "1"])
        run_command(["gmx", "mdrun", "-deffnm", os.path.join(work_dir, f"npt_{temp}")])

        # MD Production
        md_tpr = os.path.join(work_dir, f"md_{temp}.tpr")
        run_command(["gmx", "grompp", "-f", md_mdp,
                     "-c", os.path.join(work_dir, f"npt_{temp}.gro"),
                     "-t", os.path.join(work_dir, f"npt_{temp}.cpt"),
                     "-p", topol, "-o", md_tpr])
        run_command(["gmx", "mdrun", "-deffnm", os.path.join(work_dir, f"md_{temp}"), "-gpu_id", "0"])

        print(f"✅ Done with {temp} K")

    print(f"\n🎉 All simulations complete for {base}")



def process_pdbs(pdb_ids):
    ensure_dir(PDB_DIR)
    ensure_dir(RESULTS_DIR)
    for pdb_id in pdb_ids:
        try:
            print(f"\n📦 Processing {pdb_id} from RCSB...")
            pdb_file = download_pdb(pdb_id)  # Save as pdbs/XXXX.pdb
            clean_file = os.path.join(PDB_DIR, f"{pdb_id}_clean.pdb")
            clean_pdb(pdb_file, clean_file)
            prepare_gromacs_system(clean_file)
        except Exception as e:
            print(f"❌ Failed for {pdb_id}: {e}")

def process_local_pdbs(pdb_dir):
    ensure_dir(pdb_dir)
    for file_name in os.listdir(pdb_dir):
        if file_name.endswith(".pdb"):
            try:
                pdb_path = os.path.join(pdb_dir, file_name)
                print(f"\n📁 Processing local PDB: {file_name}")
                clean_file = os.path.join(pdb_dir, f"{os.path.splitext(file_name)[0]}_clean.pdb")
                print(f"\n📁 Creating local PDB: {clean_file}")
                clean_pdb(pdb_path, clean_file)
                print("cleaned pdb")
                prepare_gromacs_system(clean_file)
            except Exception as e:
                print(f"❌ Failed for {file_name}: {e}")

def process_fasta(csv_file):
    ensure_dir(IMMUNE_OUTPUT_DIR)
    try:
        print(f"\n🧬 Running ImmuneBuilder on sequences from '{csv_file}'...")
        Immunebuilder.run_pipeline(csv_file, IMMUNE_OUTPUT_DIR)
        print(f"✅ ImmuneBuilder finished. Output saved in '{IMMUNE_OUTPUT_DIR}'")
    except Exception as e:
        print(f"❌ ImmuneBuilder failed: {e}")

if __name__ == "__main__":
    print("running ImmuneBuilder Gromacs")
    parser = argparse.ArgumentParser(description="GROMACS / ImmuneBuilder Pipeline Runner")
    parser.add_argument("--pdbs", nargs="+", help="List of PDB IDs (e.g. 1HZH 4PE5)")
    parser.add_argument("--pdbdir", type=str, help="Directory containing local PDB files")
    parser.add_argument("--fasta", type=str, help="CSV file with paired sequences for ImmuneBuilder")

    args = parser.parse_args()

    if args.pdbs:
        print(f"processing pdbs {args.pdbs}")        
        process_pdbs(args.pdbs)

    if args.pdbdir:
        print(f"processing local directory {args.pdbdir}")
        process_local_pdbs(args.pdbdir)

    if args.fasta:
        print(f"processing fastas {args.fasta}")        
        process_fasta(args.fasta)

    if not args.pdbs and not args.fasta and not args.pdbdir:
        print("⚠️ Please provide at least one of: --pdbs, --pdbdir, or --fasta to run the pipeline.")