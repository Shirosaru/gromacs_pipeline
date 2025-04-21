import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import gromacs

# ------------ Helper Functions ------------ #

def core_surface(temp, directory):
    print(f"Calculating core/surface SASA at {temp}K...")
    core_group = '19'
    surface_group = '20'

    gromacs.sasa(f=os.path.join(directory, f'md_final_{temp}.xtc'),
                 s=os.path.join(directory, f'md_final_{temp}.gro'),
                 o=f'sasa_core_{temp}.xvg', n='index.ndx', input=[core_group])
    
    gromacs.sasa(f=os.path.join(directory, f'md_final_{temp}.xtc'),
                 s=os.path.join(directory, f'md_final_{temp}.gro'),
                 o=f'sasa_surface_{temp}.xvg', n='index.ndx', input=[surface_group])

def order_lambda(master_dict, mab, temps, block_length, start):
    print("Calculating Lambda values across temperatures...")

    s2_means = []
    temps_float = [float(t) for t in temps]

    for temp in temps_float:
        s2_values = master_dict[int(temp)]
        mean_s2 = np.mean(list(s2_values.values()))
        s2_means.append(mean_s2)

    slope, intercept = np.polyfit(temps_float, s2_means, 1)
    lambda_value = -slope

    print(f"Lambda (disorder coefficient): {lambda_value:.5f}")

    plt.figure(figsize=(6,4))
    plt.plot(temps_float, s2_means, marker='o')
    plt.title(f'Lambda Estimation for {mab}')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Avg S²')
    plt.grid(True)
    plt.savefig(f'lambda_plot_{mab}.png')
    plt.close()

    with open(f'lambda_{mab}.txt', 'w') as f:
        f.write(f'Lambda: {lambda_value:.5f}\n')
        f.write(f'Temperatures: {temps_float}\n')
        f.write(f'S² values: {s2_means}\n')

# ------------ Main Analysis Function ------------ #

def find_trajectory_file(directory, temp):
    """Returns the path to the trajectory file (.xtc or .trr), or None if neither exists."""
    xtc = os.path.join(directory, f'md_{temp}.xtc')
    trr = os.path.join(directory, f'md_{temp}.trr')
    return xtc if os.path.exists(xtc) else (trr if os.path.exists(trr) else None)


def analyze_trajectories(directory, temps, eq_time, mabs, block_length=10):
    print(f"Analyzing trajectories in '{directory}' for {len(temps)} temperatures...")

    master_s2_dict = {int(temp): {} for temp in temps}

    for temp in temps:
        temp = str(temp)
        traj_file = find_trajectory_file(directory, temp)
        tpr = os.path.join(directory, f'md_{temp}.tpr')
        gro = os.path.join(directory, f'md_{temp}.gro')
        index = os.path.join(directory, 'index.ndx')

        if not traj_file:
            raise FileNotFoundError(f"No trajectory (.xtc or .trr) found for {temp}K in {directory}")

        # Global analyses
        gromacs.sasa(f=traj_file, s=gro, o=f'sasa_{temp}.xvg', input=['1'])
        gromacs.hbond(f=traj_file, s=tpr, num=f'bonds_{temp}.xvg', input=['1', '1'])
        gromacs.rms(f=traj_file, s=gro, o=f'rmsd_{temp}.xvg', input=['3', '3'])
        gromacs.gyrate(f=traj_file, s=gro, o=f'gyr_{temp}.xvg', n=index, input=['1'])

        for i in range(12, 19):
            gromacs.sasa(f=traj_file, s=gro, o=f'sasa_group{i}_{temp}.xvg', n=index, input=[str(i)])
            gromacs.gyrate(f=traj_file, s=gro, o=f'gyr_group{i}_{temp}.xvg', n=index, input=[str(i)])
            gromacs.rmsf(f=traj_file, s=gro, o=f'rmsf_group{i}_{temp}.xvg',
                         n=index, b=str(eq_time * 1000), input=[str(i), str(i)])

        gromacs.hbond(f=traj_file, s=tpr, num=f'bonds_lh_{temp}.xvg', n=index, input=['10', '11'])

        trjconv_out = os.path.join(directory, f'md_final_covar_{temp}.xtc')
        gromacs.trjconv(f=traj_file, s=tpr, dt='0', fit='rot+trans',
                        n=index, o=trjconv_out, input=['1', '1'])

        gromacs.covar(f=trjconv_out, s=tpr, n=index,
                      o=f'covar_{temp}.xvg',
                      av=f'avg_covar{temp}.pdb',
                      ascii=f'covar_matrix_{temp}.dat',
                      v=f'covar_{temp}.trr', input=['4', '4'])

        gromacs.anaeig(f=trjconv_out, v=f'covar_{temp}.trr', entropy=True,
                       temp=temp, s=tpr, nevskip='6', n=index,
                       b=str(eq_time * 1000), input=[f'> sconf_{temp}.log'])

        for i in range(12, 19):
            gromacs.potential(f=traj_file, s=tpr, spherical=True, sl='10',
                              o=f'potential_group{i}_{temp}.xvg',
                              oc=f'charge_group{i}_{temp}.xvg',
                              of=f'field_group{i}_{temp}.xvg',
                              n=index, input=[str(i)])

        gromacs.dipoles(f=traj_file, s=tpr, o=f'dipole_{temp}.xvg', n=index, input=['1'])

        # Order parameters
        s2_blocks_dict = order_s2(mab=mabs, temp=temp, block_length=block_length, start=eq_time)
        master_s2_dict[int(temp)] = avg_s2_blocks(s2_blocks_dict)

        # Core vs surface
        core_surface(temp, directory)

    # Lambda calc
    order_lambda(master_dict=master_s2_dict, mab=mabs, temps=temps,
                 block_length=block_length, start=eq_time)

# ------------ CLI Interface ------------ #

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analyze GROMACS MD trajectories.")
    parser.add_argument('--dir', required=True, help='Directory containing MD files')
    parser.add_argument('--temp', required=True, help='Comma-separated temperatures (e.g., 300,310,350)')
    parser.add_argument('--eq_time', type=int, default=20, help='Equilibration time in ns')
    parser.add_argument('--mab', required=True, help='Molecule name / system identifier')
    parser.add_argument('--block_length', type=int, default=10, help='Block length for S² averaging')

    args = parser.parse_args()

    temps = [t.strip() for t in args.temp.split(',')]
    analyze_trajectories(args.dir, temps, args.eq_time, args.mab, args.block_length)
