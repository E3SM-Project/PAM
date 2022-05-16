
import subprocess
import os

import xarray as xr
import numpy as np


input_template = """
simSteps  : {nsteps}
crm_nx       : {nx}
crm_ny       : {ny}
crm_nz       : 1

initData : {init}

dtphys: {dt}
crm_per_phys: 1
# Output filename
dycore_out_prefix  : {output}

# Output frequency in number of CRM time steps
outSteps: {nsteps}
statSteps: {nsteps}

tstype: ssprk3

# Number of ensembles
nens : 1

# Parallel decomposition
nprocx : 1
nprocy : 1
"""

if __name__ == "__main__":
    nlevels = 4

    errs0 = np.empty(nlevels, dtype=np.float64)
    errsc = np.empty(nlevels, dtype=np.float64)
    errsv = np.empty(nlevels, dtype=np.float64)

    #init = "gaussian1d"
    init = "sol2d"
    for l in range(nlevels):

        nx = 10 * 2 ** l
        if init == "gaussian1d":
            ny = 2
            timeend = 0.5
        else:
            ny = nx
            timeend = 1.0

        dx = 100 / nx
        dy = 100 / ny

        scale = 10 / nx

        dt = scale * 0.1
        nsteps = int(np.ceil(timeend / dt))

        print(f"nsteps = {nsteps}, dt = {dt}")
        input = "input_conv_linswe"

        output = f"output_{l}_"
        open(input, "w").write(input_template.format(init=init, nsteps=nsteps, nx=nx,ny=ny,dt=dt,output=output))
        subprocess.run(["./driver", input], capture_output=True)

        DS = xr.open_dataset(f"{output}0.nc")
        DS.load()

        it = -1
        HC = DS.h[it, 0, 0, :, :, 0]
        HEC = DS["HEC diag"][it, 0, 0, :, :, 0]
        H0 = DS["H0 diag"][it, 0, 0, :, :, 0]
        HE0 = DS["HE0 diag"][it, 0, 0, :, :, 0]
        V = DS.v[it, 0, 0, :, :, 0]
        VE = DS["VE diag"][it, 0, 0, :, :, 0]

        errsc[l] = np.linalg.norm(HC - HEC)
        errsv[l] = np.linalg.norm(V - VE)
        errs0[l] = np.sqrt(np.sum(dx * dy * (H0 - HE0) ** 2))
        print(f"l = {l}, errc = {errsc[l]}, err0 = {errs0[l]}, errv = {errsv[l]}")
        
        DS.close()
        os.remove(input)
        os.remove(f"{output}0.nc")
    
    ratesc = np.log2(errsc[0:nlevels-1] / errsc[1:nlevels])
    ratesv = np.log2(errsv[0:nlevels-1] / errsv[1:nlevels])
    rates0 = np.log2(errs0[0:nlevels-1] / errs0[1:nlevels])

    print(ratesc)
    print(rates0)
    print(ratesv)


