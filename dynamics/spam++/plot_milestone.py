import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Circle


matplotlib.rcParams.update({'font.size': 30})

########   Thermal Shallow Water #######
recon_types = ['WENO','CFV']
recon_orders = [5,]

Nlist = [20,]

DSdict = {}
energy_dict = {}
mass_dict = {}
massfct_dict = {}
densmax_dict = {}
densmin_dict = {}
densfctmax_dict = {}
densfctmin_dict = {}
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        DSdict[simname] = xr.open_dataset('TSWE/' + simname + '/output.nc')
        DSdict[simname].load()
        
        energy_dict[simname] = DSdict[simname].energy
        mass_dict[simname] = DSdict[simname].mass
        massfct_dict[simname] = DSdict[simname].massfct
        densmax_dict[simname] = DSdict[simname].densmax
        densmin_dict[simname] = DSdict[simname].densmin
        densfctmax_dict[simname] = DSdict[simname].densfctmax
        densfctmin_dict[simname] = DSdict[simname].densfctmin  
                      
# Stats Plotting
plt.figure(figsize=(20,16))
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        energy = energy_dict[simname]
        plt.plot( (energy.isel(energy_ndofs=0) - energy.isel(energy_ndofs=0)[0])/energy.isel(energy_ndofs=0)[0]*100., label=simname)        
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in Energy')
plt.tight_layout()
plt.savefig('total_energy-tswe.png')
plt.close('all')


plt.figure(figsize=(20,16))
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        mass = mass_dict[simname]
        plt.plot( (mass.isel(mass_ndofs=0) - mass.isel(mass_ndofs=0)[0])/mass.isel(mass_ndofs=0)[0]*100., label=simname)        
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in Mass')
plt.tight_layout()
plt.savefig('total_mass-tswe.png')
plt.close('all')

plt.figure(figsize=(20,16))
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        mass = mass_dict[simname]
        plt.plot( (mass.isel(mass_ndofs=1) - mass.isel(mass_ndofs=1)[0])/mass.isel(mass_ndofs=1)[0]*100., label=simname)          
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in Bouyancy')
plt.tight_layout()
plt.savefig('total_bouyancy-tswe.png')
plt.close('all')

plt.figure(figsize=(20,16))
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        mass = mass_dict[simname]
        plt.plot( (mass.isel(mass_ndofs=2) - mass.isel(mass_ndofs=2)[0])/mass.isel(mass_ndofs=2)[0]*100., label=simname)        
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in T1')
plt.tight_layout()
plt.savefig('total_T1-tswe.png')
plt.close('all')

plt.figure(figsize=(20,16))
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        massfct = massfct_dict[simname]
        plt.plot( (massfct.isel(massfct_ndofs=0) - massfct.isel(massfct_ndofs=0)[0])/massfct.isel(massfct_ndofs=0)[0]*100., label=simname)        
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Fractional Change in T2')
plt.tight_layout()
plt.savefig('total_T2-tswe.png')
plt.close('all')

plt.figure(figsize=(20,16))
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        densmax = densmax_dict[simname]
        plt.plot(densmax.isel(densmax_ndofs=2), label=simname)        
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Maximum')
plt.tight_layout()
plt.savefig('max_T1-tswe.png')
plt.close('all')

plt.figure(figsize=(20,16))
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        densmin = densmin_dict[simname]
        plt.plot(densmin.isel(densmin_ndofs=2), label=simname)        
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Minimum')
plt.tight_layout()
plt.savefig('min_T1-tswe.png')
plt.close('all')


plt.figure(figsize=(20,16))
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        densfctmax = densfctmax_dict[simname]
        plt.plot(densfctmax.isel(densfctmax_ndofs=0), label=simname)        
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Maximum')
plt.tight_layout()
plt.savefig('max_T2-tswe.png')
plt.close('all')

plt.figure(figsize=(20,16))
for recon_type in recon_types:
    for recon_order in recon_orders:
        simname = recon_type + str(recon_order)
        densfctmin = densfctmin_dict[simname]
        plt.plot(densfctmin.isel(densfctmin_ndofs=0), label=simname)        
plt.legend()
plt.xlabel('Nsteps')
plt.ylabel('Minimum')
plt.tight_layout()
plt.savefig('min_T2-tswe.png')
plt.close('all')


simname = 'WENO' + str(5)
q = DSdict[simname].q
dens = DSdict[simname].dens
densfct = DSdict[simname].densfct
densl = DSdict[simname].densl
densfctl = DSdict[simname].densfctl


plt.figure(figsize=(10,8))
plt.contourf(dens.isel(t=0,dens_ndofs=0,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('h-tswe-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(dens.isel(t=0,dens_ndofs=1,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('S-tswe-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(dens.isel(t=0,dens_ndofs=2,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('T1-tswe-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densfct.isel(t=0,densfct_ndofs=0,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('T2-tswe-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(q.isel(t=0,q_ndofs=0,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('q-tswe-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densl.isel(t=0,densl_ndofs=1,primal_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('s-tswe-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densl.isel(t=0,densl_ndofs=2,primal_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('t1-tswe-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
plt.contourf(densfctl.isel(t=0,densfctl_ndofs=0,primal_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('t2-tswe-0.png',transparent=True)
plt.close('all')



nhoriz = len(recon_orders)
nvert = len(recon_types)
for time in Nlist:
    l = 1
    plt.figure(figsize=(10*nhoriz,8*nvert))
    for recon_type in recon_types:
        for recon_order in recon_orders:
            simname = recon_type + str(recon_order)
            dens = DSdict[simname].dens
            plt.subplot(nvert,nhoriz,l)
            plt.contourf(dens.isel(t=time,dens_ndofs=0,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
            plt.colorbar()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(simname)
            l = l+1
    plt.tight_layout()
    plt.savefig('h-tswe-' + str(time) + '.png',transparent=True)
    plt.close('all')

    l = 1
    plt.figure(figsize=(10*nhoriz,8*nvert))
    for recon_type in recon_types:
        for recon_order in recon_orders:
            simname = recon_type + str(recon_order)
            dens = DSdict[simname].dens
            plt.subplot(nvert,nhoriz,l)
            plt.contourf(dens.isel(t=time,dens_ndofs=1,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
            plt.colorbar()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(simname)
            l = l+1
    plt.tight_layout()
    plt.savefig('S-tswe-' + str(time) + '.png',transparent=True)
    plt.close('all')

    l = 1
    plt.figure(figsize=(10*nhoriz,8*nvert))
    for recon_type in recon_types:
        for recon_order in recon_orders:
            simname = recon_type + str(recon_order)
            dens = DSdict[simname].dens
            plt.subplot(nvert,nhoriz,l)
            plt.contourf(dens.isel(t=time,dens_ndofs=2,dual_ncells_z=0), np.linspace(-0.15,1.05,25) * 1e10,  cmap=plt.get_cmap('viridis'))
            plt.colorbar()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(simname)
            l = l+1
    plt.tight_layout()
    plt.savefig('T1-tswe-' + str(time) + '.png',transparent=True)
    plt.close('all')

    l = 1
    plt.figure(figsize=(10*nhoriz,8*nvert))
    for recon_type in recon_types:
        for recon_order in recon_orders:
            simname = recon_type + str(recon_order)
            densfct = DSdict[simname].densfct
            plt.subplot(nvert,nhoriz,l)
            plt.contourf(densfct.isel(t=time,densfct_ndofs=0,dual_ncells_z=0), np.linspace(0,1.05,22) * 1e10, cmap=plt.get_cmap('viridis'))
            plt.colorbar()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(simname)
            l = l+1
    plt.tight_layout()
    plt.savefig('T2-tswe-' + str(time) + '.png',transparent=True)
    plt.close('all')

    l = 1
    plt.figure(figsize=(10*nhoriz,8*nvert))
    for recon_type in recon_types:
        for recon_order in recon_orders:
            simname = recon_type + str(recon_order)
            q = DSdict[simname].q
            plt.subplot(nvert,nhoriz,l)
            plt.contourf(q.isel(t=time,q_ndofs=0,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
            plt.colorbar()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(simname)
            l = l+1
    plt.tight_layout()
    plt.savefig('q-tswe-' + str(time) + '.png',transparent=True)
    plt.close('all')
    
    l = 1
    plt.figure(figsize=(10*nhoriz,8*nvert))
    for recon_type in recon_types:
        for recon_order in recon_orders:
            simname = recon_type + str(recon_order)
            densl = DSdict[simname].densl
            plt.subplot(nvert,nhoriz,l)
            plt.contourf(densl.isel(t=time,densl_ndofs=1,primal_ncells_z=0), np.linspace(9.76,10.32,15), cmap=plt.get_cmap('viridis'))
            plt.colorbar()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(simname)
            l = l+1
    plt.tight_layout()
    plt.savefig('s-tswe-' + str(time) + '.png',transparent=True)
    plt.close('all')

    l = 1
    plt.figure(figsize=(10*nhoriz,8*nvert))
    for recon_type in recon_types:
        for recon_order in recon_orders:
            simname = recon_type + str(recon_order)
            densl = DSdict[simname].densl
            plt.subplot(nvert,nhoriz,l)
            plt.contourf(densl.isel(t=time,densl_ndofs=2,primal_ncells_z=0), np.linspace(-0.0008,0.0056,17), cmap=plt.get_cmap('viridis'))
            plt.colorbar()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(simname)
            l = l+1
    plt.tight_layout()
    plt.savefig('t1-tswe-' + str(time) + '.png',transparent=True)
    plt.close('all')

    l = 1
    plt.figure(figsize=(10*nhoriz,8*nvert))
    for recon_type in recon_types:
        for recon_order in recon_orders:
            simname = recon_type + str(recon_order)
            densfctl = DSdict[simname].densfctl
            plt.subplot(nvert,nhoriz,l)
            plt.contourf(densfctl.isel(t=time,densfctl_ndofs=0,primal_ncells_z=0),  np.linspace(0,0.0056,15), cmap=plt.get_cmap('viridis'))
            plt.colorbar()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(simname)
            l = l+1
    plt.tight_layout()
    plt.savefig('t2-tswe-' + str(time) + '.png',transparent=True)
    plt.close('all')
    
########   Compressible Euler Plotting #######

recon_types = ['WENOFUNC','CFV']
recon_orders = [3,5,7]
diff_orders = [2,]
sizes = [100,]

Nlist = [9,]

DSdict = {}
energy_dict = {}
mass_dict = {}
for recon_type in recon_types:
    for recon_order in recon_orders:
        for diff_order in diff_orders:
            for size in sizes:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                DSdict[simname] = xr.open_dataset('RB/RB' + simname + '/output.nc')
                DSdict[simname].load()
                
                energy_dict[simname] = DSdict[simname].energy
                mass_dict[simname] = DSdict[simname].mass

# Stats Plotting
for size in sizes:
    for diff_order in diff_orders:
        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                shortname = recon_type + str(recon_order)
                energy = energy_dict[simname]
                plt.plot( (energy.isel(energy_ndofs=0) - energy.isel(energy_ndofs=0)[0])/energy.isel(energy_ndofs=0)[0]*100., label=shortname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Energy')
        plt.tight_layout()
        plt.savefig('total_energy-CE.png')
        plt.close('all')


        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                shortname = recon_type + str(recon_order)
                mass = mass_dict[simname]
                plt.plot( (mass.isel(mass_ndofs=0) - mass.isel(mass_ndofs=0)[0])/mass.isel(mass_ndofs=0)[0]*100., label=shortname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Mass')
        plt.tight_layout()
        plt.savefig('total_mass-CE.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                shortname = recon_type + str(recon_order)
                mass = mass_dict[simname]
                plt.plot( (mass.isel(mass_ndofs=1) - mass.isel(mass_ndofs=1)[0])/mass.isel(mass_ndofs=1)[0]*100., label=shortname)          
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Entropic Density')
        plt.tight_layout()
        plt.savefig('total_entropic_var_density-CE.png')
        plt.close('all')
        

plt.figure(figsize=(10,8))
simname = str(100) + '-' + 'WENOFUNC' + str(7) + '-' + 'HODGE' + str(2)
q = DSdict[simname].q
plt.contourf(q.isel(t=0,q_ndofs=0,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('q-CE-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
simname = str(100) + '-' + 'WENOFUNC' + str(7) + '-' + 'HODGE' + str(2)
densl = DSdict[simname].densl
plt.contourf(densl.isel(t=0,densl_ndofs=1,primal_ncells_z=0) - 300.0, cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('thetaprime-CE-0.png',transparent=True)
plt.close('all')

            
nhoriz = len(recon_orders)
nvert = len(recon_types)
for time in Nlist:
    for size in sizes:
        for diff_order in diff_orders:
            l = 1
            plt.figure(figsize=(10*nhoriz,8*nvert))
            for recon_type in recon_types:
                for recon_order in recon_orders:
                    simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    figname = recon_type + str(recon_order)
                    q = DSdict[simname].q
                    plt.subplot(nvert,nhoriz,l)
                    plt.contourf(q.isel(t=time,q_ndofs=0,dual_ncells_z=0), np.linspace(-0.2,0.2,41), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('q-CE-' + str(time) + '.png',transparent=True)
            plt.close('all')

# THIS ASSUMES THAT ENTROPIC VARIABLE IS POTTEMP DENSITY, AND DOES A VERY ROUGH SUBTRACTION OF 300K TO GET PERTURBATION
            l = 1
            plt.figure(figsize=(10*nhoriz,8*nvert))
            for recon_type in recon_types:
                for recon_order in recon_orders:
                    simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    figname = recon_type + str(recon_order)
                    densl = DSdict[simname].densl
                    plt.subplot(nvert,nhoriz,l)
                    plt.contourf(densl.isel(t=time,densl_ndofs=1,primal_ncells_z=0) - 300.0, np.linspace(-0.08,0.64,19), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('thetaprime-CE-' + str(time) + '.png',transparent=True)
            plt.close('all')
            
            
########   Moist Compressible Euler Plotting #######

recon_types = ['WENOFUNC', 'CFV']
recon_orders = [3,5,7]
diff_orders = [2,] 
sizes = [100,] 

Nlist = [9,]

DSdict = {}
energy_dict = {}
mass_dict = {}
massfct_dict = {}
densmax_dict = {}
densmin_dict = {}
densfctmax_dict = {}
densfctmin_dict = {}
for recon_type in recon_types:
    for recon_order in recon_orders:
        for diff_order in diff_orders:
            for size in sizes:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                DSdict[simname] = xr.open_dataset('MRB/MRB' + simname + '/output.nc')
                DSdict[simname].load()
                
                energy_dict[simname] = DSdict[simname].energy
                mass_dict[simname] = DSdict[simname].mass
                massfct_dict[simname] = DSdict[simname].massfct
                densmax_dict[simname] = DSdict[simname].densmax
                densmin_dict[simname] = DSdict[simname].densmin
                densfctmax_dict[simname] = DSdict[simname].densfctmax
                densfctmin_dict[simname] = DSdict[simname].densfctmin

# Stats Plotting
for size in sizes:
    for diff_order in diff_orders:
        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                shortname = recon_type + str(recon_order)
                energy = energy_dict[simname]
                plt.plot( (energy.isel(energy_ndofs=0) - energy.isel(energy_ndofs=0)[0])/energy.isel(energy_ndofs=0)[0]*100., label=shortname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Energy')
        plt.tight_layout()
        plt.savefig('total_energy-MCE.png')
        plt.close('all')


        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                shortname = recon_type + str(recon_order)
                mass = mass_dict[simname]
                plt.plot( (mass.isel(mass_ndofs=0) - mass.isel(mass_ndofs=0)[0])/mass.isel(mass_ndofs=0)[0]*100., label=shortname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Mass')
        plt.tight_layout()
        plt.savefig('total_mass-MCE.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                shortname = recon_type + str(recon_order)
                mass = mass_dict[simname]
                plt.plot( (mass.isel(mass_ndofs=1) - mass.isel(mass_ndofs=1)[0])/mass.isel(mass_ndofs=1)[0]*100., label=shortname)          
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Entropic Density')
        plt.tight_layout()
        plt.savefig('total_entropic_var_density-MCE.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                shortname = recon_type + str(recon_order)
                massfct = massfct_dict[simname]
                plt.plot( (massfct.isel(massfct_ndofs=0) - massfct.isel(massfct_ndofs=0)[0])/massfct.isel(massfct_ndofs=0)[0]*100., label=shortname)          
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Vapor Mass')
        plt.tight_layout()
        plt.savefig('total_vapor-MCE.png')
        plt.close('all')        

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                shortname = recon_type + str(recon_order)
                densfctmax = densfctmax_dict[simname]
                plt.plot(densfctmax.isel(densfctmax_ndofs=0), label=shortname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Maximum')
        plt.tight_layout()
        plt.savefig('max_vapor-MCE.png')
        plt.close('all')
        
        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                shortname = recon_type + str(recon_order)                
                densfctmin = densfctmin_dict[simname]
                plt.plot(densfctmin.isel(densfctmin_ndofs=0), label=shortname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Minimum')
        plt.tight_layout()
        plt.savefig('min_vapor-MCE.png')
        plt.close('all')
        
plt.figure(figsize=(10,8))
simname = str(100) + '-' + 'WENOFUNC' + str(7) + '-' + 'HODGE' + str(2)
q = DSdict[simname].q
plt.contourf(q.isel(t=0,q_ndofs=0,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('q-MCE-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
simname = str(100) + '-' + 'WENOFUNC' + str(7) + '-' + 'HODGE' + str(2)
densl = DSdict[simname].densl
plt.contourf(densl.isel(t=0,densl_ndofs=1,primal_ncells_z=0) - 300.0, cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('thetaprime-MCE-0.png',transparent=True)
plt.close('all')

plt.figure(figsize=(10,8))
simname = str(100) + '-' + 'WENOFUNC' + str(7) + '-' + 'HODGE' + str(2)
densfctl = DSdict[simname].densfctl
plt.contourf(densfctl.isel(t=0,densfctl_ndofs=0,primal_ncells_z=0), cmap=plt.get_cmap('viridis'))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('qv-MCE-0.png',transparent=True)
plt.close('all')
            
nhoriz = len(recon_orders)
nvert = len(recon_types)
for time in Nlist:
    for size in sizes:
        for diff_order in diff_orders:
            l = 1
            plt.figure(figsize=(10*nhoriz,8*nvert))
            for recon_type in recon_types:
                for recon_order in recon_orders:
                    simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    figname = recon_type + str(recon_order)
                    q = DSdict[simname].q
                    plt.subplot(nvert,nhoriz,l)
                    plt.contourf(q.isel(t=time,q_ndofs=0,dual_ncells_z=0), np.linspace(-0.2,0.2,41), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('q-MCE-' + str(time) + '.png',transparent=True)
            plt.close('all')

# THIS ASSUMES THAT ENTROPIC VARIABLE IS POTTEMP DENSITY, AND DOES A VERY ROUGH SUBTRACTION OF 300K TO GET PERTURBATION
            l = 1
            plt.figure(figsize=(10*nhoriz,8*nvert))
            for recon_type in recon_types:
                for recon_order in recon_orders:
                    simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    figname = recon_type + str(recon_order)
                    densl = DSdict[simname].densl
                    plt.subplot(nvert,nhoriz,l)
                    plt.contourf(densl.isel(t=time,densl_ndofs=1,primal_ncells_z=0) - 300.0, np.linspace(-0.08,0.64,19), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('thetaprime-MCE-' + str(time) + '.png',transparent=True)
            plt.close('all')
            
            l = 1
            plt.figure(figsize=(10*nhoriz,8*nvert))
            for recon_type in recon_types:
                for recon_order in recon_orders:
                    simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    figname = recon_type + str(recon_order)
                    densfctl = DSdict[simname].densfctl
                    plt.subplot(nvert,nhoriz,l)
                    plt.contourf(densfctl.isel(t=time,densfctl_ndofs=0,primal_ncells_z=0), np.linspace(0,0.028,15), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('qv-MCE-' + str(time) + '.png',transparent=True)
            plt.close('all')