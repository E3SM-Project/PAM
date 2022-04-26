#pragma once

#include "common.h"
#include "thermo.h"
#include "DataManager.h"
#include "MultipleFields.h"
#include "pam_coupler.h" //Has DataManager and pam_const
using pam::PamCoupler;

class VariableSet
{
public: 
  
  std::string dens_name[ndensity];      // Name of each density
  std::string dens_desc[ndensity];      // Description of each density
  //bool1d                   dens_pos;       // Whether each density is positive-definite
  bool                   dens_pos[ndensity];       // Whether each density is positive-definite
  
  int dm_id_vap = -1;
  int dm_id_liq = -1;
  int dm_id_ice = -1;

  bool ice_found = false;
  bool liquid_found = false;

  ThermoPotential *thermo;
  Geometry *primal_geometry;
  Geometry *dual_geometry;
  Topology *primal_topology;
  Topology *dual_topology;
  
  void initialize(PamCoupler &coupler, ModelParameters &params, ThermoPotential &thermodynamics, Topology &primal_topo, Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom)
  {

    this->thermo = &thermodynamics;
    this->primal_geometry = &primal_geom;
    this->dual_geometry = &dual_geom;
    this->primal_topology = &primal_topo;
    this->dual_topology = &dual_topo;
  
//dens_pos IS NOT BEING PROPERLY DEALLOCATED AT THE END OF THE RUN IE WHEN THE POOL IS DESTROYED
//THIS IS REALLY WEIRD
    // Allocate device arrays for whether densities are positive-definite
    //this->dens_pos       = bool1d("dens_pos"      ,ndensity);
    //boolHost1d dens_pos_host      ("dens_pos_host"      ,ndensity);
    
    for (int l=ndensity_dycore; l<ndensity_nophysics; l++)
    {
      //dens_pos_host(l) = params.dycore_tracerpos[l-ndensity_dycore];      
      dens_pos[l] = params.dycore_tracerpos[l-ndensity_dycore];      
    }
    
    std::vector<std::string> tracer_names_loc = coupler.get_tracer_names();
    bool water_vapor_found = false;
    for (int tr=0; tr < ntracers_physics; tr++) {
      bool found, positive, adds_mass;
      std::string desc;
      coupler.get_tracer_info( tracer_names_loc[tr] , desc , found , positive , adds_mass );
      dens_name[tr+ndensity_nophysics] = tracer_names_loc[tr];
      dens_desc[tr+ndensity_nophysics] = desc;
      //dens_pos_host      (tr+ndensity_nophysics) = positive ;
      dens_pos      [tr+ndensity_nophysics] = positive ;
      if (tracer_names_loc[tr] == std::string("water_vapor")) {
        dm_id_vap = tr;
        water_vapor_found = true;
      }
      if (tracer_names_loc[tr] == std::string("cloud_liquid") || tracer_names_loc[tr] == std::string("cloud_water")) {
        dm_id_liq = tr;
        liquid_found = true;
      }
      if (tracer_names_loc[tr] == std::string("cloud_ice")) {
        dm_id_ice = tr;
        ice_found = true;
      }
    }
    if (ntracers_physics>0)
    {
    if (! water_vapor_found) {endrun("ERROR: processed registered tracers, and water_vapor was not found");}
    //if (! liquid_found) {endrun("ERROR: processed registered tracers, and cloud_liquid was not found");}
    //if (! ice_found) {endrun("ERROR: processed registered tracers, and cloud_ice was not found");}
    }
    
    for (int i=ndensity_dycore; i<ndensity_nophysics; i++)
    {
      dens_name[i] = "Tracer" + std::to_string(i-ndensity_dycore);
      dens_desc[i] = "Dycore Tracer" + std::to_string(i-ndensity_dycore);
    }
    
    for (int i=0;i < ndensity; i++)
    {
      serial_print("Density" + std::to_string(i) + " Name: " + dens_name[i]  + " Desc: " + dens_desc[i] + " Pos: " +  std::to_string(dens_pos[i]), params.masterproc);
    }

    //dens_pos_host      .deep_copy_to(dens_pos      );
    yakl::fence();
  }
  
  virtual real YAKL_INLINE get_total_density(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};
  virtual real YAKL_INLINE get_dry_density(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};
  virtual real YAKL_INLINE get_entropic_var(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};
  virtual void YAKL_INLINE set_density(real dens, real dry_dens, real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};
  virtual void YAKL_INLINE set_entropic_density(real entropic_var_density, real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};
  virtual real YAKL_INLINE get_alpha(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};
  virtual real YAKL_INLINE get_qd(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};
  virtual real YAKL_INLINE get_qv(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};
  virtual real YAKL_INLINE get_ql(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};
  virtual real YAKL_INLINE get_qi(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {};

  void convert_dynamics_to_coupler_state(PamCoupler &coupler, const FieldSet<nprognostic> &prog_vars, const FieldSet<nconstant> &const_vars) {};
  void convert_coupler_to_dynamics_state(PamCoupler &coupler, FieldSet<nprognostic> &prog_vars, const FieldSet<nconstant> &const_vars) {};
};


class VariableSet_Couple : public VariableSet
{
public: 
  void convert_dynamics_to_coupler_state(PamCoupler &coupler, const FieldSet<nprognostic> &prog_vars, const FieldSet<nconstant> &const_vars)
  {

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;
    
    real4d dm_dens_dry = coupler.dm.get<real,4>( "density_dry"      );
    real4d dm_uvel     = coupler.dm.get<real,4>( "uvel"             );
    real4d dm_vvel     = coupler.dm.get<real,4>( "vvel"             );
    real4d dm_wvel     = coupler.dm.get<real,4>( "wvel"             );
    real4d dm_temp     = coupler.dm.get<real,4>( "temp"             );

    pam::MultipleFields<ntracers_physics,real4d> dm_tracers;
    for (int tr = 0; tr < ntracers_physics; tr++) {
      auto trac = coupler.dm.get<real,4>( dens_name[tr+ndensity_nophysics] );
      dm_tracers.add_field( trac );
    }

    parallel_for( "Dynamics to Coupler State" , Bounds<4>(dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA (int k, int j, int i, int n) {
      

  // IN 3D THIS IS MORE COMPLICATED
  
       //real uvel_l  = prog_vars.fields_arr[VVAR].data(0,k+pks,j+pjs,i+pis,n) / pgeom.get_area_10entity(k+pks,j+pjs,i+pis);
       //real uvel_r  = prog_vars.fields_arr[VVAR].data(0,k+pks,j+pjs,i+pis+1,n) / pgeom.get_area_10entity(k+pks,j+pjs,i+pis+1);
//DO SOMETHING DIFFERENT AT K=0 ....
       //real wvel_u  = prog_vars.fields_arr[WVAR].data(0,k+pks,j+pjs,i+pis,n) / pgeom.get_area_01entity(k+pks,j+pjs,i+pis);
       //real wvel_d  = prog_vars.fields_arr[WVAR].data(0,k+pks-1,j+pjs,i+pis,n) / pgeom.get_area_01entity(k+pks-1,j+pjs,i+pis);
   //EVENTUALLY FIX THIS FOR 3D...
       //real vvel  = 0.0_fp;
       
      real qd = get_qd(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);
      real qv = get_qv(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);
      real alpha = get_alpha(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);
      real entropic_var = get_entropic_var(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);   
      
      real ql = 0.0_fp;
      if (liquid_found) {ql = get_qv(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);}

      real qi = 0.0_fp;
      if (ice_found) {qi = get_qi(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);}
      
      real temp = thermo->compute_T(alpha, entropic_var, qd, qv, ql, qi);
        
      dm_dens_dry(k,j,i,n) = get_dry_density(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n) / dual_geometry->get_area_11entity(k+dks,j+djs,i+dis);
      // dm_uvel    (k,j,i,n) = (uvel_l + uvel_r) * 0.5_fp;
      // dm_vvel    (k,j,i,n) = vvel;
      // dm_wvel    (k,j,i,n) = (wvel_u + wvel_d) * 0.5_fp;
      dm_temp    (k,j,i,n) = temp;
      
      for (int tr=ndensity_nophysics; tr < ndensity; tr++) {
        dm_tracers(tr-ndensity_nophysics,k,j,i,n) = prog_vars.fields_arr[DENSVAR].data(tr,k+dks,j+djs,i+dis,n) / dual_geometry->get_area_11entity(k+dks,j+djs,i+dis);
      }
    });
    
  }


  void convert_coupler_to_dynamics_state(PamCoupler &coupler, FieldSet<nprognostic> &prog_vars, const FieldSet<nconstant> &const_vars)
  {

    int dis = dual_topology->is;
    int djs = dual_topology->js;
    int dks = dual_topology->ks;
        
    auto dm_dens_dry = coupler.dm.get<real const,4>( "density_dry"      );
    auto dm_uvel     = coupler.dm.get<real const,4>( "uvel"             );
    auto dm_vvel     = coupler.dm.get<real const,4>( "vvel"             );
    auto dm_wvel     = coupler.dm.get<real const,4>( "wvel"             );
    auto dm_temp     = coupler.dm.get<real const,4>( "temp"             );

    pam::MultipleFields<ntracers_physics,realConst4d> dm_tracers;
    for (int tr = 0; tr < ntracers_physics; tr++) {
      auto trac = coupler.dm.get<real const,4>( dens_name[tr+ndensity_nophysics] );
      dm_tracers.add_field( trac );
    }

    parallel_for( "Coupler to Dynamics State" , Bounds<4>(dual_topology->nl, dual_topology->n_cells_y, dual_topology->n_cells_x, dual_topology->nens) , YAKL_LAMBDA (int k, int j, int i, int n) {

    real temp = dm_temp(k,j,i,n);
    
    //ADD VELOCITY COUPLING IN HERE PROPERLY
    //SHOULD REALLY BE LOOPING OVER PRIMAL EDGES, PROBABLY...
    real dens_dry = dm_dens_dry(k,j,i,n);
    real dens_vap = dm_tracers(dm_id_vap,k,j,i,n);
    real dens_liq = 0.0_fp;
    real dens_ice = 0.0_fp;
    if (liquid_found) {dens_liq = dm_tracers(dm_id_liq,k,j,i,n);}
    if (ice_found) {dens_ice = dm_tracers(dm_id_ice,k,j,i,n);}
    real dens = dens_dry + dens_ice + dens_liq + dens_vap;

    real qd = dens_dry / dens;
    real qv = dens_vap / dens;
    real ql = dens_liq / dens;
    real qi = dens_ice / dens;

    real alpha = 1.0_fp / dens;
    real entropic_var = thermo->compute_entropic_var_from_T(alpha, temp, qd, qv, ql, qi);
    
    set_density(dens * dual_geometry->get_area_11entity(k+dks,j+djs,i+dis), dens_dry * dual_geometry->get_area_11entity(k+dks,j+djs,i+dis), prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);
    set_entropic_density(entropic_var * dens * dual_geometry->get_area_11entity(k+dks,j+djs,i+dis), prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);

    for (int tr=ndensity_nophysics; tr < ndensity; tr++) {
      prog_vars.fields_arr[DENSVAR].data(tr,k+dks,j+djs,i+dis,n) = dm_tracers(tr-ndensity_nophysics,k,j,i,n) * dual_geometry->get_area_11entity(k+dks,j+djs,i+dis);
    }

    });

  }
};

//ADD ANELASTIC AND PRESSURE VERSIONS HERE ALSO
class VariableSet_SWE : public VariableSet {
public:
  void initialize(PamCoupler &coupler, ModelParameters &params, ThermoPotential &thermodynamics, Topology &primal_topo, Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom)
  {
    dens_name[0] = "h";
    dens_desc[0] = "fluid height";
    dens_pos[0] = false;
    VariableSet::initialize(coupler, params, thermodynamics, primal_topo, dual_topo, primal_geom, dual_geom);
  }
};

class VariableSet_TSWE : public VariableSet {
public:
  void initialize(PamCoupler &coupler, ModelParameters &params, ThermoPotential &thermodynamics, Topology &primal_topo, Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom)
  {
    dens_name[0] = "h";
    dens_name[1] = "S";
    dens_desc[0] = "fluid height";
    dens_desc[1] = "bouyancy density";
    dens_pos[0] = false;
    dens_pos[1] = false;
    VariableSet::initialize(coupler, params, thermodynamics, primal_topo, dual_topo, primal_geom, dual_geom);
  }
};

class VariableSet_CE : public VariableSet {
public:
  void initialize(PamCoupler &coupler, ModelParameters &params, ThermoPotential &thermodynamics, Topology &primal_topo, Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom)
  {
    dens_name[0] = "rho";
    dens_name[1] = "S";
    dens_desc[0] = "fluid density";
    dens_desc[1] = "entropic variable density";
    dens_pos[0] = false;
    dens_pos[1] = false;
    VariableSet::initialize(coupler, params, thermodynamics, primal_topo, dual_topo, primal_geom, dual_geom);
  }
  real YAKL_INLINE get_entropic_var(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(1,k+ks,j+js,i+is,n) / densvar(0,k+ks,j+js,i+is,n);
  };
  real YAKL_INLINE get_alpha(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return dual_geometry->get_area_11entity(k+ks,j+js,i+is) / densvar(0,k+ks,j+js,i+is,n);
  };
};
  
  
//We rely on physics packages ie micro to provide water species- must at least have vapor and cloud liquid

class VariableSet_MCE_rho : public VariableSet_Couple {
public:
  void initialize(PamCoupler &coupler, ModelParameters &params, ThermoPotential &thermodynamics, Topology &primal_topo, Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom)
  {
    dens_name[0] = "rho";
    dens_name[1] = "S";
    dens_desc[0] = "fluid density";
    dens_desc[1] = "entropic variable density";
    dens_pos[0] = false;
    dens_pos[1] = false;
    VariableSet::initialize(coupler, params, thermodynamics, primal_topo, dual_topo, primal_geom, dual_geom);
  }
  

  real YAKL_INLINE get_total_density(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(0,k+ks,j+js,i+is,n);
  }

  real YAKL_INLINE get_entropic_var(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(1,k+ks,j+js,i+is,n) / densvar(0,k+ks,j+js,i+is,n);
  }

  real YAKL_INLINE get_alpha(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return dual_geometry->get_area_11entity(k+ks,j+js,i+is) / densvar(0,k+ks,j+js,i+is,n);
  }

  real YAKL_INLINE get_qv(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(dm_id_vap + ndensity_nophysics,k+ks,j+js,i+is,n) / densvar(0,k+ks,j+js,i+is,n);
  }
  real YAKL_INLINE get_ql(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(dm_id_liq + ndensity_nophysics,k+ks,j+js,i+is,n) / densvar(0,k+ks,j+js,i+is,n);
  }
  real YAKL_INLINE get_qi(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(dm_id_ice + ndensity_nophysics,k+ks,j+js,i+is,n) / densvar(0,k+ks,j+js,i+is,n);
  }
  
  real YAKL_INLINE _water_dens(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    real vap_dens = densvar(dm_id_vap + ndensity_nophysics,k+ks,j+js,i+is,n);
    real liq_dens = 0.0_fp;
    real ice_dens = 0.0_fp;
    if (liquid_found) {liq_dens = densvar(dm_id_liq + ndensity_nophysics,k+ks,j+js,i+is,n);}
    if (ice_found) {ice_dens = densvar(dm_id_ice + ndensity_nophysics,k+ks,j+js,i+is,n);}
    return vap_dens + liq_dens + ice_dens;
  }   
  
  real YAKL_INLINE get_dry_density(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return (densvar(0,k+ks,j+js,i+is,n) - _water_dens(densvar, k, j, i, ks, js, is, n));
  }
  real YAKL_INLINE get_qd(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return (densvar(0,k+ks,j+js,i+is,n) - _water_dens(densvar, k, j, i, ks, js, is, n)) / densvar(0,k+ks,j+js,i+is,n);
  }

  void YAKL_INLINE set_density(real dens, real dryden, real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    densvar(0,k+ks,j+js,i+is,n) = dens;
  }
  void YAKL_INLINE set_entropic_density(real entropic_var_density, real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    densvar(1,k+ks,j+js,i+is,n) = entropic_var_density;
  }
  
  
};

class VariableSet_MCE_rhod : public VariableSet_Couple {
public:
  void initialize(PamCoupler &coupler, ModelParameters &params, ThermoPotential &thermodynamics, Topology &primal_topo, Topology &dual_topo, Geometry &primal_geom, Geometry &dual_geom)
  {
    dens_name[0] = "rho_d";
    dens_name[1] = "S";
    dens_desc[0] = "fluid dry density";
    dens_desc[1] = "entropic variable density";
    dens_pos[0] = false;
    dens_pos[1] = false;
    VariableSet::initialize(coupler, params, thermodynamics, primal_topo, dual_topo, primal_geom, dual_geom);
  }

real YAKL_INLINE get_total_density(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
  real dry_dens = densvar(0,k+ks,j+js,i+is,n);
  real vap_dens = densvar(dm_id_vap + ndensity_nophysics,k+ks,j+js,i+is,n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  if (liquid_found) {liq_dens = densvar(dm_id_liq + ndensity_nophysics,k+ks,j+js,i+is,n);}
  if (ice_found) {ice_dens = densvar(dm_id_ice + ndensity_nophysics,k+ks,j+js,i+is,n);}
  return dry_dens + vap_dens + liq_dens + ice_dens;
}
  
  real YAKL_INLINE get_dry_density(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(0,k+ks,j+js,i+is,n);
  }
  real YAKL_INLINE get_entropic_var(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(1,k+ks,j+js,i+is,n) / get_total_density(densvar, k, j, i, ks, js, is, n);
  }
  
  real YAKL_INLINE get_alpha(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return dual_geometry->get_area_11entity(k+ks,j+js,i+is) / get_total_density(densvar, k, j, i, ks, js, is, n);
  }
  real YAKL_INLINE get_qd(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(0,k+ks,j+js,i+is,n) / get_total_density(densvar, k, j, i, ks, js, is, n);
  }
  
  real YAKL_INLINE get_qv(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(dm_id_vap + ndensity_nophysics,k+ks,j+js,i+is,n) /  get_total_density(densvar, k, j, i, ks, js, is, n);  
  }
  real YAKL_INLINE get_ql(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(dm_id_liq + ndensity_nophysics,k+ks,j+js,i+is,n) /  get_total_density(densvar, k, j, i, ks, js, is, n);  
  }
  real YAKL_INLINE get_qi(const real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    return densvar(dm_id_ice + ndensity_nophysics,k+ks,j+js,i+is,n) /  get_total_density(densvar, k, j, i, ks, js, is, n);  
  }
  void YAKL_INLINE set_density(real dens, real drydens, real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    densvar(0,k+ks,j+js,i+is,n) = drydens;
  }
  void YAKL_INLINE set_entropic_density(real entropic_var_density, real5d densvar, int k, int j, int i, int ks, int js, int is, int n) {
    densvar(1,k+ks,j+js,i+is,n) = entropic_var_density;
  }
};