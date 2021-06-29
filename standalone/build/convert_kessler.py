import xarray as xr
import numpy as np

nc  = xr.open_dataset("kessler_ml_data.nc")
nc2 = xr.Dataset()

ncol = nc["ncol"].size
nz   = nc["nz"  ].size
nt   = nc["t"   ].size

nc2.expand_dims({'nsamples':nt*nz*ncol})
nc2.expand_dims({'three':3})


print( "nc: ml_in_theta"   , np.mean( nc["ml_in_theta"  ].to_masked_array() ) )
print( "nc: ml_in_qc"      , np.mean( nc["ml_in_qc"     ].to_masked_array() ) )
print( "nc: ml_in_qv"      , np.mean( nc["ml_in_qv"     ].to_masked_array() ) )
print( "nc: ml_in_qr"      , np.mean( nc["ml_in_qr"     ].to_masked_array() ) )
print( "nc: ml_in_rho_dry" , np.mean( nc["ml_in_rho_dry"].to_masked_array() ) )
print( "nc: ml_in_z"       , np.mean( nc["ml_in_z"      ].to_masked_array() ) )
print( "nc: ml_in_exner"   , np.mean( nc["ml_in_exner"  ].to_masked_array() ) )
print( "nc: ml_out_theta"  , np.mean( nc["ml_out_theta" ].to_masked_array() ) )
print( "nc: ml_out_qv"     , np.mean( nc["ml_out_qv"    ].to_masked_array() ) )
print( "nc: ml_out_qc"     , np.mean( nc["ml_out_qc"    ].to_masked_array() ) )
print( "nc: ml_out_qr"     , np.mean( nc["ml_out_qr"    ].to_masked_array() ) )


nc2["ml_in_theta"  ] = xr.DataArray( nc["ml_in_theta"  ].to_masked_array().reshape(nt*nz*ncol,3) , dims=("nsamples","three") );
nc2["ml_in_qc"     ] = xr.DataArray( nc["ml_in_qc"     ].to_masked_array().reshape(nt*nz*ncol,3) , dims=("nsamples","three") );
nc2["ml_in_qv"     ] = xr.DataArray( nc["ml_in_qv"     ].to_masked_array().reshape(nt*nz*ncol,3) , dims=("nsamples","three") );
nc2["ml_in_qr"     ] = xr.DataArray( nc["ml_in_qr"     ].to_masked_array().reshape(nt*nz*ncol,3) , dims=("nsamples","three") );
nc2["ml_in_rho_dry"] = xr.DataArray( nc["ml_in_rho_dry"].to_masked_array().reshape(nt*nz*ncol,3) , dims=("nsamples","three") );
nc2["ml_in_z"      ] = xr.DataArray( nc["ml_in_z"      ].to_masked_array().reshape(nt*nz*ncol  ) , dims=("nsamples"        ) );
nc2["ml_in_exner"  ] = xr.DataArray( nc["ml_in_exner"  ].to_masked_array().reshape(nt*nz*ncol,3) , dims=("nsamples","three") );
nc2["ml_out_theta" ] = xr.DataArray( nc["ml_out_theta" ].to_masked_array().reshape(nt*nz*ncol  ) , dims=("nsamples"        ) );
nc2["ml_out_qv"    ] = xr.DataArray( nc["ml_out_qv"    ].to_masked_array().reshape(nt*nz*ncol  ) , dims=("nsamples"        ) );
nc2["ml_out_qc"    ] = xr.DataArray( nc["ml_out_qc"    ].to_masked_array().reshape(nt*nz*ncol  ) , dims=("nsamples"        ) );
nc2["ml_out_qr"    ] = xr.DataArray( nc["ml_out_qr"    ].to_masked_array().reshape(nt*nz*ncol  ) , dims=("nsamples"        ) );


print( "nc2: ml_in_theta"   , np.mean( nc2["ml_in_theta"  ].to_masked_array() ) )
print( "nc2: ml_in_qc"      , np.mean( nc2["ml_in_qc"     ].to_masked_array() ) )
print( "nc2: ml_in_qv"      , np.mean( nc2["ml_in_qv"     ].to_masked_array() ) )
print( "nc2: ml_in_qr"      , np.mean( nc2["ml_in_qr"     ].to_masked_array() ) )
print( "nc2: ml_in_rho_dry" , np.mean( nc2["ml_in_rho_dry"].to_masked_array() ) )
print( "nc2: ml_in_z"       , np.mean( nc2["ml_in_z"      ].to_masked_array() ) )
print( "nc2: ml_in_exner"   , np.mean( nc2["ml_in_exner"  ].to_masked_array() ) )
print( "nc2: ml_out_theta"  , np.mean( nc2["ml_out_theta" ].to_masked_array() ) )
print( "nc2: ml_out_qv"     , np.mean( nc2["ml_out_qv"    ].to_masked_array() ) )
print( "nc2: ml_out_qc"     , np.mean( nc2["ml_out_qc"    ].to_masked_array() ) )
print( "nc2: ml_out_qr"     , np.mean( nc2["ml_out_qr"    ].to_masked_array() ) )

nc2.to_netcdf("kessler_ml_data_flat.nc")

