from netCDF4 import Dataset

src = Dataset("output.nc","r")
dst = Dataset("output4.nc","w",format="NETCDF4")

# copy dimensions
for name, dimension in src.dimensions.items():
    dst.createDimension(
        name, (len(dimension) if not dimension.isunlimited() else None))
# copy all file data except for the excluded
for name, variable in src.variables.items():
    x = dst.createVariable(name, variable.datatype, variable.dimensions)
    dst[name][:] = src[name][:]

src.close()
dst.close()

