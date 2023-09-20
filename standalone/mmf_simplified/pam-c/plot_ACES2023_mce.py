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

masslabels = []
massnames = []
nmass = 50

if core == 'MCE':
	nmass = 3
	masslabels = [r'$\rho$', 'S', r'$\rho_v$']
	massnames = ['rho','S','rhov']
if core == 'MAN':
	nmass = 2
	masslabels = ['S', r'$\rho_v$']
	massnames = ['S','rhov']
	
for i,label,name in zip(range(nmass),masslabels,massnames):
	plt.figure(figsize=(10,8))
	plt.plot( (mass.isel(mass_ndofs=i) - mass.isel(mass_ndofs=i)[0])/mass.isel(mass_ndofs=i)[0]*100., label=label)
	plt.legend()
	plt.xlabel('Nsteps')
	plt.ylabel('Fractional Change in '+label)
	plt.tight_layout()
	plt.savefig('total_'+name+'.png')	

plt.figure(figsize=(10,8))
if core == 'MCE':
	plt.plot(densmin.isel(densmin_ndofs=2))
if core == 'MAN':
	plt.plot(densmin.isel(densmin_ndofs=1))
plt.xlabel('Nsteps')
plt.ylabel('Water Vapor Minima')
plt.tight_layout()
plt.savefig('vapor_minima.png')

plt.figure(figsize=(10,8))
if core == 'MCE':
	plt.plot(densmax.isel(densmax_ndofs=2))
if core == 'MAN':
	plt.plot(densmax.isel(densmax_ndofs=1))
plt.xlabel('Nsteps')
plt.ylabel('Water Vapor Maxima')
plt.tight_layout()
plt.savefig('vapor_maxima.png')

plt.figure(figsize=(10,8))
plt.plot( (pv.isel(pv_ndofs=0) - pv.isel(pv_ndofs=0)[0])/pv.isel(pv_ndofs=0)[0]*100.)
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in PV')
plt.tight_layout()
plt.savefig('total_PV.png')

# VARIABLE PLOTTING

dens = DS.dens
densl = DS.densl
QXZl = DS.QXZl
total_dens = DS.total_dens

nt = DS.dims['t']
Nlist = [0,nt-1]
for n in Nlist:

	plt.figure(figsize=(10,8))
	if core == 'MCE':
		plt.contourf(densl.isel(t=n,densl_ndofs=2,primal_ncells_y=0,nens=0))
	if core == 'MAN':
		plt.contourf(densl.isel(t=n,densl_ndofs=1,primal_ncells_y=0,nens=0))
	plt.colorbar()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.tight_layout()
	plt.savefig('rhovl' + str(n) + '.png',transparent=True)
	plt.close('all')


	plt.figure(figsize=(10,8))
	if core == 'MCE':
		plt.contourf(dens.isel(t=n,dens_ndofs=1,dual_ncells_y=0,nens=0) / dens.isel(t=n,dens_ndofs=0,dual_ncells_y=0,nens=0))
	if core == 'MAN':
		plt.contourf(dens.isel(t=n,dens_ndofs=0,dual_ncells_y=0,nens=0) / total_dens.isel(t=n,total_dens_ndofs=0,dual_ncells_y=0,nens=0))
	plt.colorbar()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.tight_layout()
	plt.savefig('Sc' + str(n) + '.png',transparent=True)
	plt.close('all')


