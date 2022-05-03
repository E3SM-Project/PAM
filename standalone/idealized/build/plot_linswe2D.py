import xarray as xr
import numpy as np
import sys
import matplotlib.pyplot as plt

def plot_scalar_comp(plotname,v1, v2, i):
    fig, ax = plt.subplots(3, 1, figsize=(10,8))
    p0 = ax[0].contourf(v1)
    ax[0].contour(v1)
    p1 = ax[1].contourf(v2)
    ax[1].contour(v2)
    p2 = ax[2].contourf(v2 - v1)
    ax[2].contour(v2 - v1)
    
    plt.colorbar(p0, ax=ax[0])
    plt.colorbar(p1, ax=ax[1])
    plt.colorbar(p2, ax=ax[2])
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig(plotname + '.' + str(i) + '.png')
    plt.close('all')

DS = xr.open_dataset('test_dycore.nc')
DS.load()

nt = DS.dims['t']

hc = DS.h
hec = DS["HEC diag"]
v = DS.v

Nlist = np.arange(0,nt)
for i in Nlist:
    hc_i = hc[i, 0, 0, :,:, 0]
    hec_i = hec[i, 0, 0, :,:, 0]
    plot_scalar_comp('h_comp', hc_i, hec_i, i)
