from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from sys import argv

plot_difference = False
fnames = argv[1:]

for fname in fnames:
    dset = Dataset(fname, "r", format="NETCDF4")

    ti = -1
    _rho = 0
    _S = 1

    rho = dset["dens"][ti, _rho, :, 0 ,:, 0]
    S = dset["dens"][ti, _S, :, 0 ,:, 0]
    w = dset["w"][ti, 0, :, 0 ,:, 0]
    T = dset["T"][ti, 0, :, 0 ,:, 0]

    rho_exact = dset["dense"][ti, _rho, :, 0 ,:, 0]
    S_exact = dset["dense"][ti, _S, :, 0 ,:, 0]
    w_exact = dset["we"][ti, 0, :, 0 ,:, 0]
    T_exact = dset["Te"][ti, 0, :, 0 ,:, 0]

    rho_b = dset["densb"][ti, _rho, :, 0 ,:, 0]
    S_b = dset["densb"][ti, _S, :, 0 ,:, 0]

    nz, nx = rho.shape
    dz = 10e3 / (nz - 1)
    dx = 300e3 / nx
    
    T = T - 250
    T_exact = T_exact - 250

    w /= dz
    w_exact /= dz

    rho -= rho_b
    rho_exact -= rho_b
    
    S -= S_b
    S_exact -= S_b

    k5000 = int(5000 / dz)
    w5000 = w[k5000, :]
    w5000_exact = w_exact[k5000, :]


    fig, ax = plt.subplots(1, 1, figsize=(16, 10), sharex = "col")
    ax.set_title("w at 5km", fontsize=22)
    ax.plot(w5000, "k")
    ax.plot(w5000_exact, "k--")
    plt.savefig(fname.rstrip(".nc") + "_w5km.pdf")
    plt.clf()

    plotdata = {"T" : (T, T_exact), "w" : (w, w_exact), "rho" : (rho, rho_exact), "S" : (S, S_exact)}

    for vpair in (("T", "w"), ("rho", "S")):
        v1 = vpair[0]
        v2 = vpair[1]

        if vpair == ("T", "w") and not plot_difference:
            sl = 0.0006
            levels = np.array([sl * (i + 1) for i in range(0, 6)])
            levels = np.concatenate([-levels[::-1], levels])
        else:
            levels = 10

        cmap = "PuOr"
        fig, ax = plt.subplots(2, 1, figsize=(16, 10), sharex = "col")

        ax[0].set_title(v1, fontsize=22)
        if plot_difference:
            ax[0].contour(plotdata[v1][0] - plotdata[v1][1], colors=("k",), levels=levels)
            cset = ax[0].contourf(plotdata[v1][0] - plotdata[v1][1], cmap=cmap, levels=levels)
        else:
            ax[0].contour(plotdata[v1][1], colors=("k",), levels=levels)
            cset = ax[0].contourf(plotdata[v1][0], cmap=cmap, levels=levels)
        cbar = plt.colorbar(cset, ax = ax[0], shrink = 1.0)

        ax[1].set_title(v2, fontsize=22)
        if plot_difference:
            ax[1].contour(plotdata[v2][0] - plotdata[v2][1], colors=("k",), levels=levels)
            cset = ax[1].contourf(plotdata[v2][0] - plotdata[v2][1], cmap=cmap, levels=levels)
        else:
            ax[1].contour(plotdata[v2][1], colors=("k",), levels=levels)
            cset = ax[1].contourf(plotdata[v2][0], cmap=cmap, levels=levels)
        cbar = plt.colorbar(cset, ax = ax[1], shrink = 1.0)

        plt.savefig(fname.rstrip(".nc") + f"_contours_{v1}{v2}.pdf")
        plt.clf()
