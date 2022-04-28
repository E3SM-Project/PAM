from netCDF4 import Dataset
import numpy as np

def get_col_avg( data ) :
    nx = data.shape[2]
    i1 = int(   nx /4)
    i2 = int((3*nx)/4)
    return np.apply_over_axes( np.mean , data[:,:,i1:i2+1] , [1,2] )[:,0,0]

nc = Dataset("supercell_200_100_200.nc","r")
density_dry_col  = get_col_avg( nc.variables["density_dry" ] )
uvel_col         = get_col_avg( nc.variables["uvel"        ] )
vvel_col         = get_col_avg( nc.variables["vvel"        ] )
wvel_col         = get_col_avg( nc.variables["wvel"        ] )
temp_col         = get_col_avg( nc.variables["temp"        ] )
water_vapor_col  = get_col_avg( nc.variables["water_vapor" ] )
cloud_liquid_col = get_col_avg( nc.variables["cloud_liquid"] )

nc2 = Dataset("supercell_avg_column.nc","w")
nc2.createDimension("z",len(nc.dimensions["z"]))
vn = "z"           ;  nc2.createVariable(vn,nc.variables[vn].datatype,("z",));  nc2.variables[vn][:] = nc.variables[vn][:]
vn = "density_dry" ;  nc2.createVariable(vn,nc.variables[vn].datatype,("z",));  nc2.variables[vn][:] = density_dry_col 
vn = "uvel"        ;  nc2.createVariable(vn,nc.variables[vn].datatype,("z",));  nc2.variables[vn][:] = uvel_col        
vn = "vvel"        ;  nc2.createVariable(vn,nc.variables[vn].datatype,("z",));  nc2.variables[vn][:] = vvel_col        
vn = "wvel"        ;  nc2.createVariable(vn,nc.variables[vn].datatype,("z",));  nc2.variables[vn][:] = wvel_col        
vn = "temp"        ;  nc2.createVariable(vn,nc.variables[vn].datatype,("z",));  nc2.variables[vn][:] = temp_col        
vn = "water_vapor" ;  nc2.createVariable(vn,nc.variables[vn].datatype,("z",));  nc2.variables[vn][:] = water_vapor_col 
vn = "cloud_liquid";  nc2.createVariable(vn,nc.variables[vn].datatype,("z",));  nc2.variables[vn][:] = cloud_liquid_col
nc2.close()
