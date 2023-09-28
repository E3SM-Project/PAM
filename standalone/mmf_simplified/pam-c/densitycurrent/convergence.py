import yaml
import os
import subprocess
import numpy as np
from sys import argv

if __name__ == "__main__":
    nlevels = 2
    fpath = os.path.dirname(__file__) + "/../../inputs/pamc_idealized/pamc_input_extruded_densitycurrent.yaml"

    inputfile = yaml.safe_load(open(fpath))
    
    if(len(argv) > 1):
        rundir = argv[1]
    else:
        rundir = "dc_convergence"
    
    if(len(argv) > 2):
        base_dt = float(argv[2])
    else:
        base_dt = 2
    
    if(len(argv) > 3):
        diffusion_coeff = float(argv[3])
    else:
        diffusion_coeff = 75
    
    if(len(argv) > 4):
        driver = argv[4]
    else:
        driver = "driver"

    os.mkdir(rundir)
    os.chdir(rundir)

    timeend = 15 * 60
    outfreq = 1 * 60

    L = 51.2e3
    H = 6.4e3

    base_nz = 32

    errs = []
    dts = []
    dxs = []
    for l in range(nlevels):
        nz = base_nz * 2 ** l + 1
        nx = 8 * (nz - 1)

        dx = L / nx
        dz = H / (nz - 1)

        dt = base_dt / (2 ** l)
        steps = int(np.ceil(timeend / dt))

        ofname = f"output_{nx}_{nz}_"

        inputfile["crm_nx"] = nx
        inputfile["crm_nz"] = nz
        inputfile["dt_crm_phys"] = dt
        inputfile["sim_time"] = timeend
        inputfile["out_freq"] = outfreq
        inputfile["stat_freq"] = dt
        inputfile["dycore_out_prefix"] = ofname
        inputfile["entropicvar_diffusion_coeff"] = diffusion_coeff
        inputfile["velocity_diffusion_coeff"] = diffusion_coeff

        ifname = f"input_{nx}_{nz}"
        yaml.dump(inputfile, open(ifname, "w"))

        print(f"running (nx, nz, dx, dz, dt, steps) = ({nx}, {nz}, {dx}, {dz}, {dt}, {steps})")
        run = subprocess.run([f"../{driver}", ifname], capture_output=True, text=True)
        open(ofname + "stdout.txt", "w").write(run.stdout)
        open(ofname + "stderr.txt", "w").write(run.stderr)
