import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

DS = xr.open_dataset('output.nc')
DS.load()

h = DS.h
v = DS.v
q = DS.q
sl = DS.sl

M = DS.mass
B = DS.bouyancy
KE = DS.kinetic_energy
PE = DS.potential_energy
TE = DS.total_energy
PV = DS.potential_vorticity

g = 9.80616

def plot_stat(statname, data):
    plt.figure(figsize=(10,8))
    plt.plot( (data - data[0])/data[0]*100. )
    plt.xlabel('Nsteps')
    plt.ylabel('Fractional Change in ' + statname)
    plt.tight_layout()
    plt.savefig(statname + '.png')

def plot_rawstat(statname, data):
    plt.figure(figsize=(10,8))
    plt.plot(data)
    plt.xlabel('Nsteps')
    plt.ylabel(statname)
    plt.tight_layout()
    plt.savefig(statname + 'raw.png')

plot_stat('mass', M)
plot_stat('total_energy', TE)
plot_rawstat('kinetic_energy', KE)
plot_rawstat('potential_energy', PE)
plot_stat('pv', PV)
plot_stat('bouyancy', B)

Nlist = np.arange(21)

# FIX THIS STUFF UP...
# Probably control what to plot based on the test case run...

for i in Nlist:
    plt.figure(figsize=(10,8))
    plt.contourf(q.isel(t=i,q_ndofs=0,ncells_z=0))
    plt.colorbar()
    plt.contour(q.isel(t=i,q_ndofs=0,ncells_z=0))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig('q' + str(i) + '.png')
    plt.close('all')

    plt.figure(figsize=(10,8))
    plt.contourf(h.isel(t=i,h_ndofs=0,ncells_z=0))
    plt.colorbar()
    plt.contour(h.isel(t=i,h_ndofs=0,ncells_z=0))
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig('h' + str(i) + '.png')
    plt.close('all')

    plt.figure(figsize=(10,8))
    plt.contourf(sl.isel(t=i,sl_ndofs=0,ncells_z=0) - g)
    plt.colorbar()
    plt.contour(sl.isel(t=i,sl_ndofs=0,ncells_z=0) - g)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig('sl' + str(i) + '.png')
    plt.close('all')

    # plt.figure(figsize=(10,8))
    # plt.contourf(v.isel(t=i,v_ndofs=0,ncells_z=0))
    # plt.colorbar()
    # plt.contour(v.isel(t=i,v_ndofs=0,ncells_z=0))
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.savefig('u' + str(i) + '.png')
    # plt.close('all')
    #
    # plt.figure(figsize=(10,8))
    # plt.contourf(v.isel(t=i,v_ndofs=1,ncells_z=0))
    # plt.colorbar()
    # plt.contour(v.isel(t=i,v_ndofs=1,ncells_z=0))
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.savefig('v' + str(i) + '.png')
    # plt.close('all')
    #
    # plt.figure(figsize=(10,8))
    # plt.quiver(v.isel(t=i,v_ndofs=0,ncells_z=0), v.isel(t=i,v_ndofs=1,ncells_z=0))
    # # FIX THIS TO BE BASED ON WIND MAGNITUDE
    # plt.colorbar()
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.savefig('vvec' + str(i) + '.png')
    # plt.close('all')
