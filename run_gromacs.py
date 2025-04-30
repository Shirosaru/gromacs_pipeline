import MDP
import os
import os
import subprocess
import requests

def run_command(cmd, input=None, cwd=RESULTS_DIR):
    print(f"⚙️ Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, input=input, text=True, capture_output=True, cwd=cwd)
    if result.returncode != 0:
        print(f"❌ Error:\n{result.stderr}")
        raise RuntimeError("Command failed")
    return result.stdout

def prepare_gromacs_system(pdb_path, temps="300", time=100, ff="amber99sb-ildn", water="tip3p"):
    avail_temps = ['300', '310', '350', '373', '400']
    temps = [str(t.strip()) for t in temps.split(',')]
    base = os.path.splitext(os.path.basename(pdb_path))[0]
    work_dir = os.path.dirname(os.path.abspath(pdb_path))
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