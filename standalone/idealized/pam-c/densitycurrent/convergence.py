import yaml
import os
import subprocess
import numpy as np
from sys import argv

if __name__ == "__main__":
    nlevels = 2
    fpath = os.path.dirname(__file__) + "/../../inputs/pamc_input_extruded_densitycurrent.yaml"

    inputfile = yaml.safe_load(open(fpath))
    
    if(len(argv) > 1):
        rundir = argv[1]
    else:
        rundir = "straka_convergence"
    
    if(len(argv) > 2):
        base_dt = float(argv[2])
    else:
        base_dt = 2

    os.mkdir(rundir)
    os.chdir(rundir)

    timeend = 15 * 60
    outfreq = 1 * 60

    L = 56.2e3
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
        outsteps = int(np.ceil(outfreq / dt))

        ofname = f"output_{nx}_{nz}_"

        inputfile["crm_nx"] = nx
        inputfile["crm_nz"] = nz
        inputfile["dtphys"] = dt
        inputfile["simSteps"] = steps
        inputfile["outSteps"] = outsteps
        inputfile["statSteps"] = 1
        inputfile["dycore_out_prefix"] = ofname

        ifname = f"input_{nx}_{nz}"
        yaml.dump(inputfile, open(ifname, "w"))

        print(f"running (nx, nz, dx, dz, dt, steps) = ({nx}, {nz}, {dx}, {dz}, {dt}, {steps})")
        run = subprocess.run(["../driver", ifname], capture_output=True, text=True)
        open(ofname + "stdout.txt", "w").write(run.stdout)
        open(ofname + "stderr.txt", "w").write(run.stderr)
