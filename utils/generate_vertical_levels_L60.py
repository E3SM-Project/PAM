'''
This script is meant to generate PAM vertical coordinate data that roughly 
matches the default vertical grid used by E3SM-MMF, which is currently L60.
'''
#-------------------------------------------------------------------------------
import os,sys, getopt,numpy as np
from netCDF4 import Dataset

p0 = 1000e2
ps = 1000e2

# Set up terminal colors
class tcolor: ENDC,RED,GREEN,CYAN = '\033[0m','\033[31m','\033[32m','\033[36m'
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def main():

    #---------------------------------------------------------------------------
    # simple recipe using a list of height thicknesses
    #---------------------------------------------------------------------------
    grid_name = 'L60'
    num_smooth = 20
    dk_list = [ 12,  8,  8,  8,  8,  8,  8]
    dz_list = [100,200,400,500,1000,2e3,4e3]

    #---------------------------------------------------------------------------
    # Turn dk_list and dz_list into arrays of dz and z
    #---------------------------------------------------------------------------
    dz,zint = [],np.zeros(np.sum(dk_list)+1)
    kk = 1
    for d,dk in enumerate(dk_list):
        for k in range(dk):
            zint[kk] = zint[kk-1] + dz_list[d]
            dz.append(dz_list[d])
            kk += 1
    
    num_ilev = len(zint)
    num_mlev = num_ilev-1

    #---------------------------------------------------------------------------
    # Smoothing
    #---------------------------------------------------------------------------
    dz_smoothed   = np.zeros(num_mlev)
    zint_smoothed = np.copy(zint)

    for s in range(num_smooth):
        zs_tmp = np.copy(zint_smoothed)
        for k in range(1,num_mlev): 
            zint_smoothed[k] = ( 0.25*zs_tmp[k-1] + 0.5*zs_tmp[k] + 0.25*zs_tmp[k+1] )
    
    for k in range(0,num_mlev): dz_smoothed[k] = zint_smoothed[k+1] - zint_smoothed[k]

    zint = zint_smoothed

    ### print interface levels
    # for k in range(num_ilev): 
    #     k2 = num_ilev-k-1
    #     print(f'{k:02}  ({k2:02})    {zint[k2]:8.1f}')

    ### calculate mid-level height
    zmid = np.zeros(num_mlev)
    for k in range(num_mlev): 
        k2 = num_mlev-k-1
        zmid[k] = ( zint[k2+1] + zint[k2] ) / 2.
        print(f'{k:02}  ({k2:02})  {zmid[k]:8.1f}   {zint[k2]:8.1f} ')
    
    # exit()

    #---------------------------------------------------------------------------
    # Write to file
    #---------------------------------------------------------------------------
    ofile = f'PAM_vcoord_MMF_{grid_name}.nc'

    nc = Dataset(ofile,"w")
    nc.description = "Vertical height coordinates for PAM corresponding to default MMF L60 grid (as of 2023)"
    nc.z0 = 0
    nc.ztop = np.max(zmid)
    nc.nlev = len(zmid)
    
    nc.createDimension("num_interfaces",len(zint))
    nc_zint = nc.createVariable("vertical_interfaces","f8",("num_interfaces",))
    nc_zint.units = "m"
    nc_zint.description = "Height of vertical cell interfaces"
    nc_zint[:] = zint
    nc.close()

    print(f'\n{ofile}\n')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------