from netCDF4 import Dataset
import numpy as np

nc = Dataset("test.nc","r")
print(np.min(nc.variables["tracer_block"][:]))
