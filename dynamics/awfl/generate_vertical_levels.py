
import numpy as np
import sys, getopt
import matplotlib.pyplot as plt
from netCDF4 import Dataset

avail = ["help","function=","exp-base=","tanh-inflect=","tanh-steep=","tanh-scale=","z0=","ztop=","nlev=","show-plots","output="]

function="tanh"
tanh_inflect=0.4
tanh_steep=5.
tanh_scale=10.
exp_base=10.
z0=0.
ztop=10000.
nlev=64
showplots=False
output="vcoords.nc"

def print_help() :
  print("Usage: python generate_vertical_levels.py [...]")
  print("Options:")
  print("--function=[equal|exp|tanh]")
  print("  equal: self-explanatory")
  print("  exp:")
  print("    * This gives you a consistent coarsening as you increase in height")
  print("    * dz(k) = exp( pow( base , k/(nlev-1) ) )  for k in {0,...,nlev-1}")
  print("      [This is scaled afterward to match z0 and ztop]")
  print("    --exp-base=10")
  print("      * Top cell grid spacing will be 10x larger than bottom cell")
  print("  tanh:")
  print("    * This allows you to concentrate layers at the surface,")
  print("      then keep grid spacing fairly constant after that")
  print("    * dz(k) = ( ( tanh( (k/(nlev-1)-inflect)*steep ) + 1 ) / 2 * (scale-1) ) + 1")
  print("            for k in {0,...,nlev-1}")
  print("      [This is scaled afterward to match z0 and ztop]")
  print("    --tanh-inflect=0.4")
  print("      * tanh inflection point is 40% through the vertical indices")
  print("    --tanh-steep=5")
  print("      * How quickly grid spacing changes around tanh inflection point")
  print("        (generally between 2 and 10)")
  print("    --tanh-scale=10")
  print("      * The maximum bound on the grid spacing ratio of top cell to bottom cell")
  print("--help")
  print("  * Print this message")
  print("--z0=0        [Bottom of the vertical grid]")
  print("--ztop=10000  [Top of the vertical grid]")
  print("--nlev=64     [Number of vertical cells]")
  print("--show-plots")
  print("  * Show matplotlib plots over X11")
  print("--output=vcoords.nc  [NetCDF filename to save the coordinates to]")

print_help()

try:
  opts , args = getopt.getopt(sys.argv[1:] , "h" , avail)
except:
  print("ERROR: Invalid argument")
  sys.exit(-1)
for opt,arg in opts:
  if (opt == "-h" or opt == "--help") :
    sys.exit()
  elif opt == "--function" :
    function = arg
  elif opt == "--exp-base" :
    exp_base = float(arg)
  elif opt == "--tanh-inflect" :
    tanh_inflect = float(arg)
  elif opt == "--tanh-steep" :
    tanh_steep = float(arg)
  elif opt == "--tanh-scale" :
    tanh_scale = float(arg)
  elif opt == "--z0" :
    z0 = float(arg)
  elif opt == "--ztop" :
    ztop = float(arg)
  elif opt == "--nlev" :
    nlev = int(arg)
  elif opt == "--show-plots" :
    showplots = True
  elif opt == "--output" :
    output = arg
  else :
    print("ERROR: Invalid argument")
    sys.exit(-1)
    

print("\n***** USING THE FOLLOWING PARAMETERS *****")
print("nlev: "+str(nlev))
print("z0: "+str(z0))
print("ztop: "+str(ztop))
print("function: "+function)
if function == "exp" :
  print("  exp_base: "+str(exp_base))
if function == "tanh" :
  print("  tanh_inflect: "+str(tanh_inflect))
  print("  tanh_steep: "+str(tanh_steep))
  print("  tanh_scale: "+str(tanh_scale))
print("output: "+output)
print("show plots: "+str(showplots))

zthick = ztop - z0

dz   = np.array([0. for i in range(nlev)])
zint = np.array([0. for i in range(nlev+1)])
zmid = np.array([0. for i in range(nlev)])

zint[0] = 0;
  
for i in range(nlev) :
  zloc = float(i)/float(nlev-1)
  if (function == "exp") :
    dz[i] = pow(exp_base,zloc)
  elif (function == "tanh") :
    dz[i] = ((np.tanh((zloc-tanh_inflect)*tanh_steep)+1)/2*(tanh_scale-1))+1
  elif (function == "equal") :
    dz[i] = 1
  zint[i+1] = zint[i] + dz[i]

scale = zthick / zint[nlev]
for i in range(nlev+1) :
  zint[i] = zint[i] * scale + z0

for i in range(nlev) :
  dz[i] = zint[i+1] - zint[i]
  zmid[i] = ( zint[i+1] + zint[i] ) / 2

plt.plot(range(nlev),dz)
plt.xlabel("vertical level index")
plt.ylabel("grid spacing")
plt.savefig("grid_spacing_vs_level_index.png")
if showplots :
  plt.show()

plt.close()

skip = int( nlev / 50 )
plt.hlines(zint[::skip],0,1,linewidth=0.2)
plt.ylabel("height")
plt.savefig("grid_spacing_visual.png",dpi=600)
if showplots :
  plt.show()

print("\n**************** grid spacing ****************")
print(dz)
print("\n**************** interface locations ****************")
print(zint)

print("\nFigures saved in grid_spacing_vs_level_index.png and grid_spacing_visual.png")


nc = Dataset(output,"w")
nc.createDimension("num_interfaces",nlev+1)
nc_zint = nc.createVariable("vertical_interfaces","f8",("num_interfaces",))
nc_zint[:] = zint
nc.close()


