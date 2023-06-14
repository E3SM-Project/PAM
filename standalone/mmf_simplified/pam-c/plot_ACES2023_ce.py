import xarray as xr
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Circle


matplotlib.rcParams.update({'font.size': 30})
core = sys.argv[1]
test = sys.argv[2]
DS = xr.open_dataset(core+'-'+test+'/test_dycore0.nc')
DS.load()



# STATS PLOTTING

energy = DS.energy
pv = DS.pv
mass = DS.mass
densmax = DS.densmax
densmin = DS.densmin

plt.figure(figsize=(10,8))

plt.plot( (energy.isel(energy_ndofs=0) - energy.isel(energy_ndofs=0)[0])/energy.isel(energy_ndofs=0)[0]*100.)
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in Energy')
plt.tight_layout()
plt.savefig('total_energy.png')


plt.figure(figsize=(10,8))
if core == 'CE':
	plt.plot( (mass.isel(mass_ndofs=0) - mass.isel(mass_ndofs=0)[0])/mass.isel(mass_ndofs=0)[0]*100., label='rho')
	plt.plot( (mass.isel(mass_ndofs=1) - mass.isel(mass_ndofs=1)[0])/mass.isel(mass_ndofs=1)[0]*100., label='S')
elif core == 'AN':
	plt.plot( (mass.isel(mass_ndofs=0) - mass.isel(mass_ndofs=0)[0])/mass.isel(mass_ndofs=0)[0]*100., label='S')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in Densities')
plt.tight_layout()
plt.savefig('total_masses.png')

# VARIABLE PLOTTING

dens = DS.dens
densl = DS.densl
QXZl = DS.QXZl
total_dens = DS.total_dens

nt = DS.dims['t']
Nlist = [0,nt-1]
for n in Nlist:


	plt.figure(figsize=(10,8))
	if core == 'CE':
		plt.contourf(dens.isel(t=n,dens_ndofs=1,dual_ncells_y=0,nens=0) / dens.isel(t=n,dens_ndofs=0,dual_ncells_y=0,nens=0))
	if core == 'AN':
		plt.contourf(dens.isel(t=n,dens_ndofs=0,dual_ncells_y=0,nens=0) / total_dens.isel(t=n,total_dens_ndofs=0,dual_ncells_y=0,nens=0))	
	plt.colorbar()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.tight_layout()
	plt.savefig('Sc' + str(n) + '.png',transparent=True)
	plt.close('all')

