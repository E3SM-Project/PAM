import xarray as xr
from plot_helpers import plotvar_scalar2D, plotvar_vector2D, plot_stat, plot_rawstat
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Circle

recon_types = ['WENOFUNC',] #'CFV'
recon_orders = [1,3,5,7,9]
diff_orders = [2,] #4
sizes = [100,] #200

Nlist = [0,6,9]


DSdict = {}
energy_dict = {}
pv_dict = {}
mass_dict = {}
densmax_dict = {}
densmin_dict = {}
for recon_type in recon_types:
    for recon_order in recon_orders:
        for diff_order in diff_orders:
            for size in sizes:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                DSdict[simname] = xr.open_dataset('MRB/MRB' + simname + '/output.nc')
                DSdict[simname].load()
                
                energy_dict[simname] = DSdict[simname].energy
                pv_dict[simname] = DSdict[simname].pv
                mass_dict[simname] = DSdict[simname].mass
                densmax_dict[simname] = DSdict[simname].densmax
                densmin_dict[simname] = DSdict[simname].densmin


matplotlib.rcParams.update({'font.size': 30})


# Stats Plotting

for size in sizes:
    for diff_order in diff_orders:
        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                energy = energy_dict[simname]
                plt.plot( (energy.isel(energy_ndofs=0) - energy.isel(energy_ndofs=0)[0])/energy.isel(energy_ndofs=0)[0]*100., label=simname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Energy')
        plt.tight_layout()
        plt.savefig('total_energy.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                pv = pv_dict[simname]
                plt.plot( (pv.isel(pv_ndofs=0) - pv.isel(pv_ndofs=0)[1])/pv.isel(pv_ndofs=0)[1]*100., label=simname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in PV')
        plt.tight_layout()
        plt.savefig('total_pv.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                mass = mass_dict[simname]
                plt.plot( (mass.isel(mass_ndofs=0) - mass.isel(mass_ndofs=0)[0])/mass.isel(mass_ndofs=0)[0]*100., label=simname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Mass')
        plt.tight_layout()
        plt.savefig('total_mass.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                mass = mass_dict[simname]
                plt.plot( (mass.isel(mass_ndofs=1) - mass.isel(mass_ndofs=1)[0])/mass.isel(mass_ndofs=1)[0]*100., label=simname)          
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Entropic Density')
        plt.tight_layout()
        plt.savefig('total_entropic_var_density.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                mass = mass_dict[simname]
                plt.plot( (mass.isel(mass_ndofs=2) - mass.isel(mass_ndofs=2)[0])/mass.isel(mass_ndofs=2)[0]*100., label=simname)          
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Vapor Mass')
        plt.tight_layout()
        plt.savefig('total_vapor.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')


        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                mass = mass_dict[simname]
                plt.plot( (mass.isel(mass_ndofs=3) - mass.isel(mass_ndofs=3)[0])/mass.isel(mass_ndofs=3)[0]*100., label=simname)          
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Liquid Mass')
        plt.tight_layout()
        plt.savefig('total_liquid.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                mass = mass_dict[simname]
                plt.plot( (mass.isel(mass_ndofs=4) - mass.isel(mass_ndofs=4)[0])/mass.isel(mass_ndofs=4)[0]*100., label=simname)          
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Fractional Change in Ice Mass')
        plt.tight_layout()
        plt.savefig('total_ice.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                densmax = densmax_dict[simname]
                plt.plot(densmax.isel(densmax_ndofs=2), label=simname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Maximum')
        plt.tight_layout()
        plt.savefig('max_vapor.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')
        
        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                densmin = densmin_dict[simname]
                plt.plot(densmin.isel(densmin_ndofs=2), label=simname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Minimum')
        plt.tight_layout()
        plt.savefig('min_vapor.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                densmax = densmax_dict[simname]
                plt.plot(densmax.isel(densmax_ndofs=3), label=simname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Maximum')
        plt.tight_layout()
        plt.savefig('max_liquid.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')
        
        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                densmin = densmin_dict[simname]
                plt.plot(densmin.isel(densmin_ndofs=3), label=simname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Minimum')
        plt.tight_layout()
        plt.savefig('min_liquid.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')

        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                densmax = densmax_dict[simname]
                plt.plot(densmax.isel(densmax_ndofs=4), label=simname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Maximum')
        plt.tight_layout()
        plt.savefig('max_ice.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
        plt.close('all')
        
        plt.figure(figsize=(20,16))
        for recon_type in recon_types:
            for recon_order in recon_orders:
                simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                densmin = densmin_dict[simname]
                plt.plot(densmin.isel(densmin_ndofs=4), label=simname)        
        plt.legend()
        plt.xlabel('Nsteps')
        plt.ylabel('Minimum')
        plt.tight_layout()
        plt.savefig('min_ice.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png')
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
                    figname = recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    #q = DSdict[simname].q
                    QXZl = DSdict[simname].QXZl
                    plt.subplot(nvert,nhoriz,l)
                    ##plt.contourf(q.isel(t=time,q_ndofs=0,dual_ncells_z=0), np.linspace(-0.2,0.2,41), cmap=plt.get_cmap('viridis'))
                    #plt.contourf(q.isel(t=time,q_ndofs=0,dual_ncells_z=0), cmap=plt.get_cmap('viridis'))
                    #plt.contourf(QXZl.isel(t=time,QXZl_ndofs=0,dual_ncells_y=0), np.linspace(-0.2,0.2,41), cmap=plt.get_cmap('viridis'))
                    plt.contourf(QXZl.isel(t=time,QXZl_ndofs=0,dual_ncells_y=0), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('q' + str(time) + '.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png',transparent=True)
            plt.close('all')

# THIS ASSUMES THAT ENTROPIC VARIABLE IS POTTEMP DENSITY, AND DOES A VERY ROUGH SUBTRACTION OF 300K TO GET PERTURBATION
            l = 1
            plt.figure(figsize=(10*nhoriz,8*nvert))
            for recon_type in recon_types:
                for recon_order in recon_orders:
                    simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    figname = recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    densl = DSdict[simname].densl
                    plt.subplot(nvert,nhoriz,l)
                    specific_entropic_var = densl.isel(t=time,densl_ndofs=1,primal_ncells_y=0) / densl.isel(t=time,densl_ndofs=0,primal_ncells_y=0)
                    plt.contourf(specific_entropic_var - 300.0, cmap=plt.get_cmap('viridis'))
                    #plt.contourf(specific_entropic_var - 300.0, np.linspace(-0.08,0.64,19), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('thetaprime' + str(time) + '.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png',transparent=True)
            plt.close('all')

            l = 1
            plt.figure(figsize=(10*nhoriz,8*nvert))
            for recon_type in recon_types:
                for recon_order in recon_orders:
                    simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    figname = recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    densl = DSdict[simname].densl
                    plt.subplot(nvert,nhoriz,l)
                    qv = densl.isel(t=time,densl_ndofs=2,primal_ncells_y=0) / densl.isel(t=time,densl_ndofs=0,primal_ncells_y=0)
                    plt.contourf(qv, cmap=plt.get_cmap('viridis'))
                    #plt.contourf(qv, np.linspace(-0.08,0.64,19), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('qv' + str(time) + '.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png',transparent=True)
            plt.close('all')

            l = 1
            plt.figure(figsize=(10*nhoriz,8*nvert))
            for recon_type in recon_types:
                for recon_order in recon_orders:
                    simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    figname = recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    densl = DSdict[simname].densl
                    plt.subplot(nvert,nhoriz,l)
                    ql = densl.isel(t=time,densl_ndofs=3,primal_ncells_y=0) / densl.isel(t=time,densl_ndofs=0,primal_ncells_y=0)
                    plt.contourf(ql, cmap=plt.get_cmap('viridis'))
                    #plt.contourf(ql, np.linspace(-0.08,0.64,19), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('ql' + str(time) + '.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png',transparent=True)
            plt.close('all')

            l = 1
            plt.figure(figsize=(10*nhoriz,8*nvert))
            for recon_type in recon_types:
                for recon_order in recon_orders:
                    simname = str(size) + '-' + recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    figname = recon_type + str(recon_order) + '-' + 'HODGE' + str(diff_order)
                    densl = DSdict[simname].densl
                    plt.subplot(nvert,nhoriz,l)
                    qi = densl.isel(t=time,densl_ndofs=4,primal_ncells_y=0) / densl.isel(t=time,densl_ndofs=0,primal_ncells_y=0)
                    plt.contourf(qi, cmap=plt.get_cmap('viridis'))
                    #plt.contourf(qi, np.linspace(-0.08,0.64,19), cmap=plt.get_cmap('viridis'))
                    plt.colorbar()
                    plt.xlabel('x')
                    plt.ylabel('y')
                    plt.title(figname)
                    l = l+1
            plt.tight_layout()
            plt.savefig('qi' + str(time) + '.' + str(size) + '.' + 'HODGE' + str(diff_order) + '.png',transparent=True)
            plt.close('all')
