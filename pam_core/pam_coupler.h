
#pragma once

#include "pam_const.h"
#include "DataManager.h"
#include "vertical_interp.h"
//#include "YAKL_netcdf.h"
#include "Options.h"


namespace pam {



  YAKL_INLINE real hydrostatic_pressure( realConst3d hy_params , real z_in , real z0 , real dz ,
                                         int k, int iens             ) {
    real z = ( z_in - z0 ) / dz;
    real a0 = hy_params(k,0,iens);
    real a1 = hy_params(k,1,iens);
    real a2 = hy_params(k,2,iens);
    real a3 = hy_params(k,3,iens);
    real a4 = hy_params(k,4,iens);
    real lnp = a0 + ( a1 + ( a2 + ( a3 + ( a4)*z)*z)*z)*z;
    return exp(lnp);
  }



  YAKL_INLINE real hydrostatic_density( realConst3d hy_params , real z_in , real z0 , real dz ,
                                        int k, int iens , real grav ) {
    real z = ( z_in - z0 ) / dz;
    real a1 = hy_params(k,1,iens);
    real a2 = hy_params(k,2,iens);
    real a3 = hy_params(k,3,iens);
    real a4 = hy_params(k,4,iens);
    real p = hydrostatic_pressure( hy_params , z_in , z0 , dz , k , iens );
    real mult = a1 + (2*a2 + (3*a3 + (4*a4)*z)*z)*z;
    real dpdz = mult*p/dz;
    return -dpdz/grav;
  }



  YAKL_INLINE real compute_pressure( real rho_d, real rho_v, real T, real R_d, real R_v ) {
    return rho_d*R_d*T + rho_v*R_v*T;
  }



  class PamCoupler {
  protected:

    std::thread::id thread_id;

    Options options;

    real R_d;    // Dry air gas constant
    real R_v;    // Water vapor gas constant
    real cp_d;   // Dry air specific heat at constant pressure
    real cp_v;   // Water vapor specific heat at constant pressure
    real grav;   // Acceleration due to gravity (m s^-2): typically 9.81
    real p0;     // Reference pressure (Pa): typically 10^5
    real xlen;   // Domain length in the x-direction in meters
    real ylen;   // Domain length in the y-direction in meters
    real dt_gcm; // Time step of the GCM for this MMF invocation

    DataManager     dm;
    DataManagerHost dm_host;

    struct Tracer {
      std::string name;
      std::string desc;
      bool        positive;
      bool        adds_mass;
    };
    std::vector<Tracer> tracers;

    struct DycoreFunction {
      std::string                                   name;
      std::function< void ( PamCoupler & , real ) > func;
    };
    std::vector< DycoreFunction > dycore_functions;

    struct MMFFunction {
      std::string                                   name;
      std::function< void ( PamCoupler & , real ) > func;
    };
    std::vector< MMFFunction > pam_functions;


  public:

    PamCoupler() {
      this->R_d    = 287 ;
      this->R_v    = 461 ;
      this->cp_d   = 1004;
      this->cp_v   = 1859;
      this->grav   = 9.81;
      this->p0     = 1.e5;
      this->xlen   = -1;
      this->ylen   = -1;
      this->dt_gcm = -1;
      this->thread_id = std::this_thread::get_id();
    }


    PamCoupler(PamCoupler &&) = default;
    PamCoupler &operator=(PamCoupler &&) = default;
    PamCoupler(PamCoupler const &) = delete;
    PamCoupler &operator=(PamCoupler const &) = delete;


    ~PamCoupler() {
      dm.finalize();
      options.finalize();
      tracers = std::vector<Tracer>();
      this->R_d    = 287 ;
      this->R_v    = 461 ;
      this->cp_d   = 1004;
      this->cp_v   = 1859;
      this->grav   = 9.81;
      this->p0     = 1.e5;
      this->xlen   = -1;
      this->ylen   = -1;
      this->dt_gcm = -1;
    }


    void set_dt_gcm(real dt_gcm) { this->dt_gcm = dt_gcm; }


    std::thread::id         get_thread_id                  () const { return this->thread_id    ; }
    real                    get_R_d                        () const { return this->R_d          ; }
    real                    get_R_v                        () const { return this->R_v          ; }
    real                    get_cp_d                       () const { return this->cp_d         ; }
    real                    get_cp_v                       () const { return this->cp_v         ; }
    real                    get_grav                       () const { return this->grav         ; }
    real                    get_p0                         () const { return this->p0           ; }
    real                    get_xlen                       () const { return this->xlen         ; }
    real                    get_ylen                       () const { return this->ylen         ; }
    real                    get_dx                         () const { return get_xlen()/get_nx(); }
    real                    get_dy                         () const { return get_ylen()/get_ny(); }
    real                    get_dt_gcm                     () const { return this->dt_gcm       ; }
    DataManager const &     get_data_manager_readonly      () const { return this->dm           ; }
    DataManager       &     get_data_manager_readwrite     ()       { return this->dm           ; }
    DataManagerHost const & get_data_manager_host_readonly () const { return this->dm_host      ; }
    DataManagerHost       & get_data_manager_host_readwrite()       { return this->dm_host      ; }


    int get_nx() const {
      if (dm.find_dimension("x") == -1) return -1;
      return dm.get_dimension_size("x");
    }


    int get_ny() const {
      if (dm.find_dimension("y") == -1) return -1;
      return dm.get_dimension_size("y");
    }


    int get_nz() const {
      if (dm.find_dimension("z") == -1) return -1;
      return dm.get_dimension_size("z");
    }


    int get_nens() const {
      if (dm.find_dimension("nens") == -1) return -1;
      return dm.get_dimension_size("nens");
    }


    int get_ncrms() const {
      return get_nens();
    }


    template <class T>
    void add_option( std::string key , T value ) {
      options.add_option<T>(key,value);
    }


    template <class T>
    void set_option( std::string key , T value ) {
      options.set_option<T>(key,value);
    }


    template <class T>
    T get_option( std::string key ) const {
      return options.get_option<T>(key);
    }


    bool option_exists( std::string key ) const {
      return options.option_exists(key);
    }


    void delete_option( std::string key ) {
      options.delete_option(key);
    }


    void add_dycore_function( std::string name , std::function< void ( PamCoupler & , real ) > func ) {
      dycore_functions.push_back( { name , func } );
    }


    void add_pam_function( std::string name , std::function< void ( PamCoupler & , real ) > func ) {
      pam_functions.push_back( { name , func } );
    }


    int get_num_pam_functions() const { return pam_functions.size(); }


    std::vector<std::string> get_pam_function_names() const {
      std::vector<std::string> names;
      for (auto & func : pam_functions) {
        names.push_back(func.name);
      }
      return names;
    }


    int get_num_dycore_functions() const { return dycore_functions.size(); }


    std::vector<std::string> get_dycore_function_names() const {
      std::vector<std::string> names;
      for (auto & func : dycore_functions) {
        names.push_back(func.name);
      }
      return names;
    }


    void run_pam_function( std::string name , real dt_crm ) {
      for (int i=0; i < pam_functions.size(); i++) {
        if (name == pam_functions[i].name) {
          #ifdef PAM_FUNCTION_TRACE
            dm.clean_all_entries();
          #endif
          #ifdef PAM_FUNCTION_TIMERS
            yakl::timer_start( pam_functions[i].name.c_str() );
          #endif
          pam_functions[i].func( *this , dt_crm );
          #ifdef PAM_FUNCTION_TIMERS
            yakl::timer_stop ( pam_functions[i].name.c_str() );
          #endif
          #ifdef PAM_FUNCTION_TRACE
            auto dirty_entry_names = dm.get_dirty_entries();
            std::cout << "MMF Function " << pam_functions[i].name << " ran with a time step of "
                      << dt_crm << " seconds and wrote to the following coupler entries: ";
            for (int e=0; e < dirty_entry_names.size(); e++) {
              std::cout << dirty_entry_names[e];
              if (e < dirty_entry_names.size()-1) std::cout << ", ";
            }
            std::cout << "\n\n";
          #endif
          return;
        }
      }
      endrun("ERROR: run_pam_function called with invalid function name: " + name);
    }


    void run_pam_functions(real dt_crm) {
      for (int i=0; i < pam_functions.size(); i++) {
        #ifdef PAM_FUNCTION_TRACE
          dm.clean_all_entries();
        #endif
        #ifdef PAM_FUNCTION_TIMERS
          yakl::timer_start( pam_functions[i].name.c_str() );
        #endif
        pam_functions[i].func( *this , dt_crm );
        #ifdef PAM_FUNCTION_TIMERS
          yakl::timer_stop ( pam_functions[i].name.c_str() );
        #endif
        #ifdef PAM_FUNCTION_TRACE
          auto dirty_entry_names = dm.get_dirty_entries();
          std::cout << "MMF Function " << pam_functions[i].name << " ran with a time step of "
                    << dt_crm << " seconds and wrote to the following coupler entries: ";
          for (int e=0; e < dirty_entry_names.size(); e++) {
            std::cout << dirty_entry_names[e];
            if (e < dirty_entry_names.size()-1) std::cout << ", ";
          }
          std::cout << "\n\n";
        #endif
      }
    }


    void run_dycore_function( std::string name , real dt_dycore ) {
      for (int i=0; i < pam_functions.size(); i++) {
        if (name == pam_functions[i].name) {
          #ifdef PAM_FUNCTION_TRACE
            dm.clean_all_entries();
          #endif
          #ifdef PAM_FUNCTION_TIMERS
            yakl::timer_start( dycore_functions[i].name.c_str() );
          #endif
          dycore_functions[i].func( *this , dt_dycore );
          #ifdef PAM_FUNCTION_TIMERS
            yakl::timer_stop ( dycore_functions[i].name.c_str() );
          #endif
          #ifdef PAM_FUNCTION_TRACE
            auto dirty_entry_names = dm.get_dirty_entries();
            std::cout << "Dycore Function " << dycore_functions[i].name << " ran with a time step of "
                      << dt_dycore << " seconds and wrote to the following coupler entries: ";
            for (int e=0; e < dirty_entry_names.size(); e++) {
              std::cout << dirty_entry_names[e];
              if (e < dirty_entry_names.size()-1) std::cout << ", ";
            }
            std::cout << "\n\n";
          #endif
          return;
        }
      }
      endrun("ERROR: run_dycore_function called with invalid function name: " + name);
    }


    void run_dycore_functions(real dt_dycore) {
      for (int i=0; i < dycore_functions.size(); i++) {
        #ifdef PAM_FUNCTION_TRACE
          dm.clean_all_entries();
        #endif
        #ifdef PAM_FUNCTION_TIMERS
          yakl::timer_start( dycore_functions[i].name.c_str() );
        #endif
        dycore_functions[i].func( *this , dt_dycore );
        #ifdef PAM_FUNCTION_TIMERS
          yakl::timer_stop ( dycore_functions[i].name.c_str() );
        #endif
        #ifdef PAM_FUNCTION_TRACE
          auto dirty_entry_names = dm.get_dirty_entries();
          std::cout << "Dycore Function " << dycore_functions[i].name << " ran with a time step of "
                    << dt_dycore << " seconds and wrote to the following coupler entries: ";
          for (int e=0; e < dirty_entry_names.size(); e++) {
            std::cout << dirty_entry_names[e];
            if (e < dirty_entry_names.size()-1) std::cout << ", ";
          }
          std::cout << "\n\n";
        #endif
      }
    }


    void set_phys_constants(real R_d, real R_v, real cp_d, real cp_v, real grav=9.81, real p0=1.e5) {
      this->R_d  = R_d ;
      this->R_v  = R_v ;
      this->cp_d = cp_d;
      this->cp_v = cp_v;
      this->grav = grav;
      this->p0   = p0  ;
    }


    void set_grid(real xlen, real ylen, realConst2d zint_in) {
      using yakl::c::parallel_for;
      using yakl::c::SimpleBounds;

      int nz    = get_nz();
      int nens  = get_nens();
      this->xlen = xlen;
      this->ylen = ylen;
      auto zint = dm.get<real,2>("vertical_interface_height");
      auto dz   = dm.get<real,2>("vertical_cell_dz");
      auto zmid = dm.get<real,2>("vertical_midpoint_height" );
      parallel_for( "vert grid 1" , SimpleBounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
        zint(k,iens) = zint_in(k,iens);
        if (k < nz) {
          zmid(k,iens) = 0.5_fp * (zint_in(k,iens) + zint_in(k+1,iens));
          dz  (k,iens) = zint_in(k+1,iens) - zint_in(k,iens);
        }
      });
    }


    void set_grid(real xlen, real ylen, realConst1d zint_in) {
      using yakl::c::parallel_for;
      using yakl::c::SimpleBounds;

      int nz    = get_nz();
      int nens  = get_nens();
      this->xlen = xlen;
      this->ylen = ylen;
      auto zint = dm.get<real,2>("vertical_interface_height");
      auto dz   = dm.get<real,2>("vertical_cell_dz");
      auto zmid = dm.get<real,2>("vertical_midpoint_height" );
      parallel_for( "vert grid 2" , SimpleBounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
        zint(k,iens) = zint_in(k);
        if (k < nz) {
          zmid(k,iens) = 0.5_fp * (zint_in(k) + zint_in(k+1));
          dz  (k,iens) = zint_in(k+1) - zint_in(k);
        }
      });
    }


    
    void add_tracer( std::string tracer_name , std::string tracer_desc , bool positive , bool adds_mass ) {
      int nz   = get_nz  ();
      int ny   = get_ny  ();
      int nx   = get_nx  ();
      int nens = get_nens();
      dm.register_and_allocate<real>( tracer_name , tracer_desc , {nz,ny,nx,nens} , {"z","y","x","nens"} );
      tracers.push_back( { tracer_name , tracer_desc , positive , adds_mass } );
    }



    int get_num_tracers() const { return tracers.size(); }


    
    std::vector<std::string> get_tracer_names() const {
      std::vector<std::string> ret;
      for (int i=0; i < tracers.size(); i++) { ret.push_back( tracers[i].name ); }
      return ret;
    }


    
    void get_tracer_info(std::string tracer_name , std::string &tracer_desc, bool &tracer_found ,
                         bool &positive , bool &adds_mass) const {
      std::vector<std::string> ret;
      for (int i=0; i < tracers.size(); i++) {
        if (tracer_name == tracers[i].name) {
          positive     = tracers[i].positive ;
          tracer_desc  = tracers[i].desc     ;
          adds_mass    = tracers[i].adds_mass;
          tracer_found = true;
          return;
        }
      }
      tracer_found = false;
    }


    
    bool tracer_exists( std::string tracer_name ) const {
      for (int i=0; i < tracers.size(); i++) {
        if (tracer_name == tracers[i].name) return true;
      }
      return false;
    }



    void allocate_coupler_state( int nz, int ny, int nx, int nens ) {
      using yakl::c::parallel_for;
      using yakl::c::SimpleBounds;

      dm.register_and_allocate<real>("density_dry"              ,"dry density"                ,{nz,ny,nx,nens},{"z","y","x","nens"});
      dm.register_and_allocate<real>("uvel"                     ,"x-direction velocity"       ,{nz,ny,nx,nens},{"z","y","x","nens"});
      dm.register_and_allocate<real>("vvel"                     ,"y-direction velocity"       ,{nz,ny,nx,nens},{"z","y","x","nens"});
      dm.register_and_allocate<real>("wvel"                     ,"z-direction velocity"       ,{nz,ny,nx,nens},{"z","y","x","nens"});
      dm.register_and_allocate<real>("temp"                     ,"temperature"                ,{nz,ny,nx,nens},{"z","y","x","nens"});
      dm.register_and_allocate<real>("vertical_interface_height","vertical interface height"  ,{nz+1    ,nens},{"zp1"      ,"nens"});
      dm.register_and_allocate<real>("vertical_cell_dz"         ,"vertical grid spacing"      ,{nz      ,nens},{"z"        ,"nens"});
      dm.register_and_allocate<real>("vertical_midpoint_height" ,"vertical midpoint height"   ,{nz      ,nens},{"z"        ,"nens"});
      dm.register_and_allocate<real>("hydrostasis_parameters"   ,"hydrostasis parameters"     ,{nz,5    ,nens},{"z","nhy"  ,"nens"});
      dm.register_and_allocate<real>("gcm_density_dry"          ,"GCM column dry density"     ,{nz      ,nens},{"z"        ,"nens"});
      dm.register_and_allocate<real>("gcm_uvel"                 ,"GCM column u-velocity"      ,{nz      ,nens},{"z"        ,"nens"});
      dm.register_and_allocate<real>("gcm_vvel"                 ,"GCM column v-velocity"      ,{nz      ,nens},{"z"        ,"nens"});
      dm.register_and_allocate<real>("gcm_wvel"                 ,"GCM column w-velocity"      ,{nz      ,nens},{"z"        ,"nens"});
      dm.register_and_allocate<real>("gcm_temp"                 ,"GCM column temperature"     ,{nz      ,nens},{"z"        ,"nens"});
      dm.register_and_allocate<real>("gcm_water_vapor"          ,"GCM column water vapor mass",{nz      ,nens},{"z"        ,"nens"});

      auto density_dry  = dm.get_collapsed<real>("density_dry"              );
      auto uvel         = dm.get_collapsed<real>("uvel"                     );
      auto vvel         = dm.get_collapsed<real>("vvel"                     );
      auto wvel         = dm.get_collapsed<real>("wvel"                     );
      auto temp         = dm.get_collapsed<real>("temp"                     );
      auto zint         = dm.get_collapsed<real>("vertical_interface_height");
      auto dz           = dm.get_collapsed<real>("vertical_cell_dz"         );
      auto zmid         = dm.get_collapsed<real>("vertical_midpoint_height" );
      auto hy_params    = dm.get_collapsed<real>("hydrostasis_parameters"   );
      auto gcm_rho_d    = dm.get_collapsed<real>("gcm_density_dry"          );
      auto gcm_uvel     = dm.get_collapsed<real>("gcm_uvel"                 );
      auto gcm_vvel     = dm.get_collapsed<real>("gcm_vvel"                 );
      auto gcm_wvel     = dm.get_collapsed<real>("gcm_wvel"                 );
      auto gcm_temp     = dm.get_collapsed<real>("gcm_temp"                 );
      auto gcm_rho_v    = dm.get_collapsed<real>("gcm_water_vapor"          );

      parallel_for( "coupler zero" , SimpleBounds<1>(nz*ny*nx*nens) , YAKL_LAMBDA (int i) {
        density_dry (i) = 0;
        uvel        (i) = 0;
        vvel        (i) = 0;
        wvel        (i) = 0;
        temp        (i) = 0;
        if (i < (nz+1)*nens) zint(i) = 0;
        if (i < (nz  )*nens) {
          zmid     (i) = 0;
          dz       (i) = 0;
          gcm_rho_d(i) = 0;
          gcm_uvel (i) = 0;
          gcm_vvel (i) = 0;
          gcm_wvel (i) = 0;
          gcm_temp (i) = 0;
          gcm_rho_v(i) = 0;
        }
        if (i < nz*3  *nens) hy_params(i) = 0;
      });
    }



    void update_hydrostasis() {
      using yakl::c::parallel_for;
      using yakl::c::SimpleBounds;
      using yakl::intrinsics::matmul_cr;
      using yakl::intrinsics::matinv_ge;
      using yakl::atomicAdd;

      auto zint      = dm.get<real const,2>("vertical_interface_height");
      auto zmid      = dm.get<real const,2>("vertical_midpoint_height" );
      auto dz        = dm.get<real const,2>("vertical_cell_dz");
      auto hy_params = dm.get<real,3>("hydrostasis_parameters"   );

      auto dens_dry = dm.get<real const,4>("density_dry");
      auto dens_wv  = dm.get<real const,4>("water_vapor");
      auto temp     = dm.get<real const,4>("temp");

      int nz   = get_nz();
      int ny   = get_ny();
      int nx   = get_nx();
      int nens = get_nens();

      YAKL_SCOPE( R_d , this->R_d );
      YAKL_SCOPE( R_v , this->R_v );

      // Compute average column of pressure for each ensemble
      real2d pressure_col("pressure_col",nz,nens);
      memset( pressure_col , 0._fp );
      real r_nx_ny = 1._fp / (nx*ny);
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        real rho_d = dens_dry(k,j,i,iens);
        real rho_v = dens_wv (k,j,i,iens);
        real T     = temp    (k,j,i,iens);
        atomicAdd( pressure_col(k,iens) , compute_pressure( rho_d , rho_v , T , R_d , R_v )*r_nx_ny );
      });

      parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
        int constexpr npts = 5;
        int kmid = k;
        int kbot = k-2;
        int ktop = k+2;
        while (kbot < 0   ) { kbot++; ktop++; }
        while (ktop > nz-1) { kbot--; ktop--; }

        SArray<double,1,npts> z;
        real z0 = zmid(kmid,iens);
        for (int i=0; i < npts; i++) { z(i) = ( zmid(kbot+i,iens) - z0 ) / dz(k,iens); }

        SArray<double,2,npts,npts> vand;
        for (int j=0; j < npts; j++) {
          for (int i=0; i < npts; i++) {
            vand(j,i) = pow( z(i) , (double) j );
          }
        }

        auto vand_inv = matinv_ge( vand );

        SArray<double,1,npts> logp;
        for (int i=0; i < npts; i++) {
          logp(i) = log(pressure_col(kbot+i,iens));
        }

        auto params = matmul_cr( vand_inv , logp );

        for (int i=0; i < npts; i++) { hy_params(k,i,iens) = params(i); }
      });

    }



    real4d compute_pressure_array() const {
      using yakl::c::parallel_for;
      using yakl::c::SimpleBounds;

      auto dens_dry = dm.get<real const,4>("density_dry");
      auto dens_wv  = dm.get<real const,4>("water_vapor");
      auto temp     = dm.get<real const,4>("temp");

      int nz   = get_nz();
      int ny   = get_ny();
      int nx   = get_nx();
      int nens = get_nens();

      real4d pressure("pressure",nz,ny,nx,nens);

      YAKL_SCOPE( R_d , this->R_d );
      YAKL_SCOPE( R_v , this->R_v );

      parallel_for( "coupler pressure" , SimpleBounds<4>(nz,ny,nx,nens) ,
                    YAKL_LAMBDA (int k, int j, int i, int iens) {
        real rho_d = dens_dry(k,j,i,iens);
        real rho_v = dens_wv (k,j,i,iens);
        real T     = temp    (k,j,i,iens);
        pressure(k,j,i,iens) = compute_pressure( rho_d , rho_v , T , R_d , R_v );
      });

      return pressure;
    }


  };

}


