import yaml
import os
import subprocess
import numpy as np
from sys import argv
from netCDF4 import Dataset

def create_vcoords(name, H, nz):
    dz = H / (nz - 1)
    zint = [k * dz + dz / 2 for k in range(0, nz - 1)]
    zint = [0] + zint + [H]
    
    nc = Dataset(name, "w", format="NETCDF4")
    nc.z0 = 0
    nc.ztop = H
    nc.nlev = len(zint)
    nc.createDimension("num_interfaces", len(zint))
    nc_zint = nc.createVariable("vertical_interfaces", "f8", ("num_interfaces",))
    nc_zint[:] = zint

    nc.close()


def compute_Ediss_and_Edisp(a, b):
    cov_M = np.cov(np.vstack((np.reshape(a, np.size(a)), np.reshape(b, np.size(b)))))
    sigma_a = np.sqrt(cov_M[0, 0])
    sigma_b = np.sqrt(cov_M[1, 1])
    cov_ab = cov_M[0, 1]
    mean_a = np.mean(a)
    mean_b = np.mean(b)
    Ediss = (sigma_a - sigma_b) ** 2 + (mean_a - mean_b) ** 2
    Edisp = 2 * sigma_a * sigma_b - 2 * cov_ab
    return (Ediss, Edisp)

def compute_errors(dataset, nx, nz, dx, dz):
    _rho = 0
    _S = 1
    ti = -1

    rho = dataset["dens"][ti, _rho, :, 0 ,:, 0]
    S = dataset["dens"][ti, _S, :, 0 ,:, 0]
    w = dataset["w"][ti, 0, :, 0 ,:, 0]
    T = dataset["T"][ti, 0, :, 0 ,:, 0]

    rho_exact = dataset["dense"][ti, _rho, :, 0 ,:, 0]
    S_exact = dataset["dense"][ti, _S, :, 0 ,:, 0]
    w_exact = dataset["we"][ti, 0, :, 0 ,:, 0]
    T_exact = dataset["Te"][ti, 0, :, 0 ,:, 0]

    rho_b = dataset["densb"][ti, _rho, :, 0 ,:, 0]
    S_b = dataset["densb"][ti, _S, :, 0 ,:, 0]

    w /= dz
    w_exact /= dz

    rho /= (dx * dz)
    rho_exact /= (dx * dz)
    S /= (dx * dz)
    S_exact /= (dx * dz)
    
    Linf_rho = np.max(np.abs((rho - rho_exact)))
    L2_rho = np.sqrt(np.sum((rho - rho_exact) ** 2) / (nx * nz))
    Ediss_rho, Edisp_rho = compute_Ediss_and_Edisp(rho, rho_exact)
    
    Linf_S = np.max(np.abs((S - S_exact)))
    L2_S = np.sqrt(np.sum((S - S_exact) ** 2) / (nx * nz))
    Ediss_S, Edisp_S = compute_Ediss_and_Edisp(S, S_exact)

    Linf_T = np.max(np.abs((T - T_exact)))
    L2_T = np.sqrt(np.sum((T - T_exact) ** 2) / (nx * nz))
    Ediss_T, Edisp_T = compute_Ediss_and_Edisp(T, T_exact)

    Linf_w = np.max(np.abs((w - w_exact)))
    L2_w = np.sqrt(np.sum((w - w_exact) ** 2) / (nx * nz))
    Ediss_w, Edisp_w = compute_Ediss_and_Edisp(w, w_exact)

    ret = {"T" : (Linf_T, L2_T, Ediss_T, Edisp_T),
           "w" : (Linf_w, L2_w, Ediss_w, Edisp_w),
           "rho" : (Linf_rho, L2_rho, Ediss_rho, Edisp_rho),
           "S" : (Linf_S, L2_S, Ediss_S, Edisp_S)}
    return ret

if __name__ == "__main__":
    nlevels = 2
    fpath = os.path.dirname(__file__) + "/../../inputs/pamc_idealized/pamc_input_extruded_gravitywave.yaml"

    inputfile = yaml.safe_load(open(fpath))
    
    if(len(argv) > 1):
        rundir = argv[1]
    else:
        rundir = "gw_convergence"
    
    if(len(argv) > 2):
        base_dt = float(argv[2])
    else:
        base_dt = 20

    os.mkdir(rundir)
    os.chdir(rundir)

    timeend = 30 * 60

    L = 300e3
    H = 10e3

    base_nz = 20

    errs = []
    dts = []
    dxs = []
    for l in range(nlevels):
        nz = base_nz * 2 ** l + 1
        nx = 15 * (nz - 1)

        dx = L / nx
        dz = H / (nz - 1)

        dt = base_dt / (2 ** l)
        steps = int(np.ceil(timeend / dt))

        ofname = f"output_{nx}_{nz}_"

        vcoords = "uniform"
        #vcoords = "uniform_variable"
        #vcoords = "variable"

        if vcoords == "uniform_variable":
            vcoords = f"vcoords_{nz}.nc"
            create_vcoords(vcoords, H, nz)
            
        if vcoords == "variable":
            vcoords = f"vcoords_{nz}.nc"
            vert_levels_script = "../../../../utils/generate_vertical_levels.py"
            vert_func = "exp"
            subprocess.run(["python3", vert_levels_script, f"--nlev={nz-1}", f"--function={vert_func}",
                            f"--ztop={H}", f"--output={vcoords}"],
                           capture_output=True, text=True)


        inputfile["crm_nx"] = nx
        inputfile["crm_nz"] = nz
        inputfile["vcoords"] = vcoords
        inputfile["dt_crm_phys"] = dt
        inputfile["sim_time"] = timeend
        inputfile["out_freq"] = timeend
        inputfile["stat_freq"] = dt
        inputfile["dycore_out_prefix"] = ofname

        ifname = f"input_{nx}_{nz}"
        yaml.dump(inputfile, open(ifname, "w"))

        print(f"running (nx, nz, dx, dz, dt, steps) = ({nx}, {nz}, {dx}, {dz}, {dt}, {steps})")
        run = subprocess.run(["../driver", ifname], capture_output=True, text=True)
        open(ofname + "stdout.txt", "w").write(run.stdout)
        open(ofname + "stderr.txt", "w").write(run.stderr)

        dataset = Dataset(ofname + "0.nc", "r", format="NETCDF4")
        err = compute_errors(dataset, nx, nz, dx, dz)
        dts.append(dt)
        dxs.append(dx)
        errs.append(err)

    variables = ("T", "w", "rho", "S")
    outfiles  = {var : open(f"errors_{var}.txt", "w") for var in variables}
    header = "{:5} {:8} {:8} {:10} {:10} {:10} {:10} {:10} {:10}\n".format("lev", "dx", "dt", "Linf", "Linf_r", "L2", "L2_r", "Ediss", "Edisp")
    for var in variables:
        outfiles[var].write(header)
        for l in range(nlevels):
            err = errs[l][var]
            if l > 0:
                rate_Linf = np.log2(errs[l-1][var][0] / errs[l][var][0])
                rate_L2 = np.log2(errs[l-1][var][1] / errs[l][var][1])
            else:
                rate_Linf = 0
                rate_L2 = 0
            line = f"{l:<5} {dxs[l]:<8.2f} {dts[l]:<8.2f} {err[0]:<10.2e} {rate_Linf:<10.1e} {err[1]:<10.2e} {rate_L2:<10.1e} {err[2]:<10.2e} {err[3]:<10.2e}\n"

            outfiles[var].write(line)

