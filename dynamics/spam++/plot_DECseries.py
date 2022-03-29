import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Circle

DS_weno = xr.open_dataset('SWE/WENOFUNC9-DualWENOFUNC9-100/output.nc')
DS_cfv = xr.open_dataset('SWE/CFV10-DualCFV10-100/output.nc')
DS_low = xr.open_dataset('SWE/CFV2-DualCFV2-100/output.nc')
DS_weno.load()
DS_cfv.load()
DS_low.load()

matplotlib.rcParams.update({'font.size': 30})


# STATS PLOTTING


energy_weno = DS_weno.energy
energy_cfv = DS_cfv.energy
energy_low = DS_low.energy

pv_weno = DS_weno.pv
pv_cfv = DS_cfv.pv
pv_low = DS_low.pv

mass_weno = DS_weno.mass
densmax_weno = DS_weno.densmax
densmin_weno = DS_weno.densmin

mass_cfv = DS_cfv.mass
densmax_cfv = DS_cfv.densmax
densmin_cfv = DS_cfv.densmin

mass_low = DS_low.mass
densmax_low = DS_low.densmax
densmin_low = DS_low.densmin

massfct_weno = DS_weno.massfct
densfctmax_weno = DS_weno.densfctmax
densfctmin_weno = DS_weno.densfctmin

massfct_cfv = DS_cfv.massfct
densfctmax_cfv = DS_cfv.densfctmax
densfctmin_cfv = DS_cfv.densfctmin

plt.figure(figsize=(10,8))
plt.plot( (energy_cfv.isel(energy_ndofs=2) - energy_cfv.isel(energy_ndofs=2)[0])/energy_cfv.isel(energy_ndofs=2)[0]*100., label='cfv')
plt.plot( (energy_weno.isel(energy_ndofs=2) - energy_weno.isel(energy_ndofs=2)[0])/energy_weno.isel(energy_ndofs=2)[0]*100., label='weno')
plt.plot( (energy_low.isel(energy_ndofs=2) - energy_low.isel(energy_ndofs=2)[0])/energy_low.isel(energy_ndofs=2)[0]*100., label='low')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in Energy')
plt.tight_layout()
plt.savefig('total_energy.png')

plt.figure(figsize=(10,8))
plt.plot( (mass_cfv.isel(mass_ndofs=0) - mass_cfv.isel(mass_ndofs=0)[0])/mass_cfv.isel(mass_ndofs=0)[0]*100., label='cfv')
plt.plot( (mass_weno.isel(mass_ndofs=0) - mass_weno.isel(mass_ndofs=0)[0])/mass_weno.isel(mass_ndofs=0)[0]*100., label='weno')
plt.plot( (mass_low.isel(mass_ndofs=0) - mass_low.isel(mass_ndofs=0)[0])/mass_low.isel(mass_ndofs=0)[0]*100., label='low')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in Mass')
plt.tight_layout()
plt.savefig('total_M.png')

plt.figure(figsize=(10,8))
plt.plot( (mass_cfv.isel(mass_ndofs=1) - mass_cfv.isel(mass_ndofs=1)[0])/mass_cfv.isel(mass_ndofs=1)[0]*100., label='cfv')
plt.plot( (mass_weno.isel(mass_ndofs=1) - mass_weno.isel(mass_ndofs=1)[0])/mass_weno.isel(mass_ndofs=1)[0]*100., label='weno')
plt.plot( (mass_low.isel(mass_ndofs=1) - mass_low.isel(mass_ndofs=1)[0])/mass_low.isel(mass_ndofs=1)[0]*100., label='low')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in S')
plt.tight_layout()
plt.savefig('total_S.png')

plt.figure(figsize=(10,8))
plt.plot( (massfct_cfv.isel(massfct_ndofs=0) - massfct_cfv.isel(massfct_ndofs=0)[0])/massfct_cfv.isel(massfct_ndofs=0)[0]*100., label='cfv')
plt.plot( (massfct_weno.isel(massfct_ndofs=0) - massfct_weno.isel(massfct_ndofs=0)[0])/massfct_weno.isel(massfct_ndofs=0)[0]*100., label='weno')
plt.plot( (mass_low.isel(mass_ndofs=2) - mass_low.isel(mass_ndofs=2)[0])/mass_low.isel(mass_ndofs=2)[0]*100., label='low')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in T')
plt.tight_layout()
plt.savefig('total_T.png')

plt.figure(figsize=(10,8))
plt.plot( (pv_cfv.isel(pv_ndofs=0) - pv_cfv.isel(pv_ndofs=0)[0])/pv_cfv.isel(pv_ndofs=0)[0]*100., label='cfv')
plt.plot( (pv_weno.isel(pv_ndofs=0) - pv_weno.isel(pv_ndofs=0)[0])/pv_weno.isel(pv_ndofs=0)[0]*100., label='weno')
plt.plot( (pv_low.isel(pv_ndofs=0) - pv_weno.isel(pv_ndofs=0)[0])/pv_weno.isel(pv_ndofs=0)[0]*100., label='low')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in PV')
plt.tight_layout()
plt.savefig('total_PV.png')



plt.figure(figsize=(10,8))
plt.plot( densmax_cfv.isel(densmax_ndofs=1), label='cfv')
plt.plot( densmax_weno.isel(densmax_ndofs=1), label='weno')
plt.plot( densmax_low.isel(densmax_ndofs=1), label='low')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Max')
plt.tight_layout()
plt.savefig('S_max.png')

plt.figure(figsize=(10,8))
plt.plot( densmin_cfv.isel(densmin_ndofs=1), label='cfv')
plt.plot( densmin_weno.isel(densmin_ndofs=1), label='weno')
plt.plot( densmin_low.isel(densmin_ndofs=1), label='low')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Min')
plt.tight_layout()
plt.savefig('S_min.png')

plt.figure(figsize=(10,8))
plt.plot( densfctmax_cfv.isel(densfctmax_ndofs=0), label='cfv')
plt.plot( densfctmax_weno.isel(densfctmax_ndofs=0), label='weno')
plt.plot( densmax_low.isel(densmax_ndofs=2), label='low')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Max')
plt.tight_layout()
plt.savefig('T_max.png')

plt.figure(figsize=(10,8))
plt.plot( densfctmin_cfv.isel(densfctmin_ndofs=0), label='cfv')
plt.plot( densfctmin_weno.isel(densfctmin_ndofs=0), label='weno')
plt.plot( densmin_low.isel(densmin_ndofs=2), label='low')
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Min')
plt.tight_layout()
plt.savefig('T_min.png')























# VARIABLE PLOTTING

dens_weno = DS_weno.dens
dens_cfv = DS_cfv.dens
dens_low = DS_low.dens

densfct_weno = DS_weno.densfct
densfct_cfv = DS_cfv.densfct

densl_weno = DS_weno.densl
densl_cfv = DS_cfv.densl
densl_low = DS_low.densl

densfctl_weno = DS_weno.densfctl
densfctl_cfv = DS_cfv.densfctl

q_weno = DS_weno.q
q_cfv = DS_cfv.q
q_low = DS_low.q

# at t=0

plt.figure(figsize=(10,8))
plt.contourf(dens_cfv.isel(t=0,dens_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(dens_cfv.isel(t=0,dens_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('h0.png',transparent=True)
plt.close('all')


plt.figure(figsize=(10,8))
plt.contourf(dens_cfv.isel(t=0,dens_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(dens_cfv.isel(t=0,dens_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('S0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densl_cfv.isel(t=0,densl_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(densl_cfv.isel(t=0,densl_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('sl0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densfct_cfv.isel(t=0,densfct_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(densfct_cfv.isel(t=0,densfct_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('T0.png',transparent=True)
plt.close('all')


plt.figure(figsize=(10,8))
plt.contourf(densfctl_cfv.isel(t=0,densfctl_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(densfctl_cfv.isel(t=0,densfctl_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('tl0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(q_cfv.isel(t=0,q_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(q_cfv.isel(t=0,q_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('q0.png',transparent=True)
plt.close('all')



# at t=20
plt.figure(figsize=(10,8))
plt.contourf(dens_cfv.isel(t=20,dens_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(dens_cfv.isel(t=20,dens_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('h20-cfv.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(dens_weno.isel(t=20,dens_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(dens_weno.isel(t=20,dens_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('h20-weno.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(dens_low.isel(t=20,dens_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(dens_low.isel(t=20,dens_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('h20-low.png',transparent=True)
plt.close('all')






plt.figure(figsize=(10,8))
plt.contourf(dens_cfv.isel(t=20,dens_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(dens_cfv.isel(t=20,dens_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('S20-cfv.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(dens_weno.isel(t=20,dens_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(dens_weno.isel(t=20,dens_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('S20-weno.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(dens_low.isel(t=20,dens_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(dens_low.isel(t=20,dens_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('S20-low.png',transparent=True)
plt.close('all')






plt.figure(figsize=(10,8))
plt.contourf(densl_cfv.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(densl_cfv.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('sl20-cfv.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densl_weno.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(densl_weno.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('sl20-weno.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densl_low.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.colorbar()
#plt.contour(densl_low.isel(t=20,densl_ndofs=1,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('sl20-low.png',transparent=True)
plt.close('all')






plt.figure(figsize=(10,8))
plt.contourf(densfct_cfv.isel(t=20,densfct_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(densfct_cfv.isel(t=20,densfct_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('T20-cfv.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densfct_weno.isel(t=20,densfct_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(densfct_weno.isel(t=20,densfct_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('T20-weno.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(dens_low.isel(t=20,dens_ndofs=2,ncells_z=0))
plt.colorbar()
#plt.contour(dens_low.isel(t=20,dens_ndofs=2,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('T20-low.png',transparent=True)
plt.close('all')





plt.figure(figsize=(10,8))
plt.contourf(densfctl_cfv.isel(t=20,densfctl_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(densfctl_cfv.isel(t=20,densfctl_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('tl20-cfv.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densfctl_weno.isel(t=20,densfctl_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(densfctl_weno.isel(t=20,densfctl_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('tl20-weno.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densl_low.isel(t=20,densl_ndofs=2,ncells_z=0))
plt.colorbar()
#plt.contour(densl_low.isel(t=20,densl_ndofs=2,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('tl20-low.png',transparent=True)
plt.close('all')





plt.figure(figsize=(10,8))
plt.contourf(q_cfv.isel(t=20,q_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(q_cfv.isel(t=20,q_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('q20-cfv.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(q_weno.isel(t=20,q_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(q_weno.isel(t=20,q_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('q20-weno.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(q_low.isel(t=20,q_ndofs=0,ncells_z=0))
plt.colorbar()
#plt.contour(q_low.isel(t=20,q_ndofs=0,ncells_z=0))
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('q20-low.png',transparent=True)
plt.close('all')



