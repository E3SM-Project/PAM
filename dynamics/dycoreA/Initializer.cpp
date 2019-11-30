
#include "Initializer.h"


void initialize_mpi( int *argc , char ***argv , Parallel &par ) {
  int ierr = MPI_Init( argc , argv );
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);

  //Determine if I'm the master process
  if (par.myrank == 0) {
    par.masterproc = 1;
  } else {
    par.masterproc = 0;
  }
}

void initialize(realArr &state, Domain &dom, Parallel &par, Exchange &exch, TimeIntegrator &tint) {
  int ierr;
  SArray<real,ord> gllOrdPoints;
  SArray<real,ord> gllOrdWeights;
  SArray<real,tord> gllTordPoints;
  SArray<real,tord> gllTordWeights;

  if (dom.ny_glob == 1) { dom.run2d = 1; } else { dom.run2d = 0; }

  //Get GLL points and weights
  TransformMatrices<real> trans;
  trans.get_gll_points(gllOrdPoints);
  trans.get_gll_weights(gllOrdWeights);
  trans.get_gll_points(gllTordPoints);
  trans.get_gll_weights(gllTordWeights);
  if (par.nranks != par.nproc_x*par.nproc_y) {
    std::cerr << "ERROR: nproc_x*nproc_y != nranks\n";
    std::cerr << par.nproc_x << " " << par.nproc_y << " " << par.nranks << "\n";
    exit(-1);
  }

  //Get my x and y process grid ID
  par.px = par.myrank % par.nproc_x;
  par.py = par.myrank / par.nproc_x;

  //Get my beginning and ending global indices
  double nper;
  nper = ((double) dom.nx_glob)/par.nproc_x;
  par.i_beg = (long) round( nper* par.px    );
  par.i_end = (long) round( nper*(par.px+1) )-1;
  nper = ((double) dom.ny_glob)/par.nproc_y;
  par.j_beg = (long) round( nper* par.py    );
  par.j_end = (long) round( nper*(par.py+1) )-1;
  //Determine my number of grid cells
  dom.nx = par.i_end - par.i_beg + 1;
  dom.ny = par.j_end - par.j_beg + 1;
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      int pxloc = par.px+i-1;
      if (pxloc < 0            ) pxloc = pxloc + par.nproc_x;
      if (pxloc > par.nproc_x-1) pxloc = pxloc - par.nproc_x;
      int pyloc = par.py+j-1;
      if (pyloc < 0            ) pyloc = pyloc + par.nproc_y;
      if (pyloc > par.nproc_y-1) pyloc = pyloc - par.nproc_y;
      par.neigh(j,i) = pyloc * par.nproc_x + pxloc;
    }
  }

  // Debug output for the parallel decomposition
  if (0) {
    for (int rr=0; rr < par.nranks; rr++) {
      if (rr == par.myrank) {
        std::cout << "Hello! My Rank is what, my rank is who, my rank is: " << par.myrank << "\n";
        std::cout << "My proc grid ID is: " << par.px << " , " << par.py << "\n";
        std::cout << "I have: " << dom.nx << " x " << dom.ny << " columns." << "\n";
        std::cout << "I start at index: " << par.i_beg << " x " << par.j_beg << "\n";
        std::cout << "My neighbor matrix is:\n";
        for (int j = 2; j >= 0; j--) {
          for (int i = 0; i < 3; i++) {
            std::cout << std::setw(6) << par.neigh(j,i) << " ";
          }
          printf("\n");
        }
        printf("\n");
      }
      ierr = MPI_Barrier(MPI_COMM_WORLD);
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }

  // Initialize the grid
  dom.etime = 0;
  dom.nz = dom.nz_glob;
  dom.dx = dom.xlen / dom.nx_glob;
  dom.dy = dom.ylen / dom.ny_glob;
  dom.dz = dom.zlen / dom.nz_glob;

  tint.initialize(dom);

  // Allocate the MPI exchange buffers
  exch.allocate(dom);

  // Allocate the fluid state variable
  state = realArr( "state" , numState , dom.nz+2*hs , dom.ny+2*hs , dom.nx+2*hs );

  dom.hyDensCells      = realArr( "hyCellsR"  , dom.nz+2*hs );
  dom.hyDensThetaCells = realArr( "hyCellsRT" , dom.nz+2*hs );
  dom.hyThetaCells     = realArr( "hyCellsT"  , dom.nz+2*hs );
  dom.hyPressureCells  = realArr( "hyCellsp"  , dom.nz+2*hs );

  dom.hyDensGLL      = realArr( "hyGLLR"  , dom.nz , tord );
  dom.hyDensThetaGLL = realArr( "hyGLLRT" , dom.nz , tord );
  dom.hyThetaGLL     = realArr( "hyGLLT"  , dom.nz , tord );
  dom.hyPressureGLL  = realArr( "hyGLLp"  , dom.nz , tord );

  // Initialize the hydrostatic background state for cell averages
  // for (int k=0; k<dom.nz; k++) {
  yakl::parallel_for( "initHydroCells" , dom.nz , YAKL_LAMBDA (int const k) {
    dom.hyDensCells     (hs+k) = 0;
    dom.hyDensThetaCells(hs+k) = 0;
    dom.hyThetaCells    (hs+k) = 0;
    dom.hyPressureCells (hs+k) = 0;
    // Perform ord-point GLL quadrature for the cell averages
    for (int kk=0; kk<ord; kk++) {
      real zloc = (k + 0.5_fp)*dom.dz + gllOrdPoints(kk)*dom.dz;
      real r0, t0;

      if (dom.dataInit == DATA_INIT_THERMAL || dom.dataInit == DATA_INIT_COLLISION || dom.dataInit == DATA_INIT_STRAKA) {
        t0 = 300._fp;
        hydro::hydroConstTheta( t0 , zloc , r0 );
      }

      dom.hyDensCells     (hs+k) += gllOrdWeights(kk) * r0;
      dom.hyDensThetaCells(hs+k) += gllOrdWeights(kk) * r0*t0;
      dom.hyThetaCells    (hs+k) += gllOrdWeights(kk) * t0;
      dom.hyPressureCells (hs+k) += gllOrdWeights(kk) * C0*pow( r0*t0 , GAMMA );
    }
  });

  // Enforce vertical boundaries
  // for (int ii=0; ii<hs; ii++) {
  yakl::parallel_for( "boundariesHydroCells" , hs , YAKL_LAMBDA (int const ii) {
    dom.hyDensCells     (ii) = dom.hyDensCells     (hs);
    dom.hyDensThetaCells(ii) = dom.hyDensThetaCells(hs);
    dom.hyThetaCells    (ii) = dom.hyThetaCells    (hs);
    dom.hyPressureCells (ii) = dom.hyPressureCells (hs);
    dom.hyDensCells     (dom.nz+hs+ii) = dom.hyDensCells     (dom.nz+hs-1);
    dom.hyDensThetaCells(dom.nz+hs+ii) = dom.hyDensThetaCells(dom.nz+hs-1);
    dom.hyThetaCells    (dom.nz+hs+ii) = dom.hyThetaCells    (dom.nz+hs-1);
    dom.hyPressureCells (dom.nz+hs+ii) = dom.hyPressureCells (dom.nz+hs-1);
  });

  // Initialize the hydrostatic background state for GLL points
  // Perform ord-point GLL quadrature for the cell averages
  // for (int k=0; k<dom.nz; k++) {
  //   for (int kk=0; kk<tord; kk++) {
  yakl::parallel_for( "initHydroGLL" , dom.nz,tord , YAKL_LAMBDA (int k, int kk) {
    real zloc = (k + 0.5_fp)*dom.dz + gllTordPoints(kk)*dom.dz;
    real r0, t0;

    if (dom.dataInit == DATA_INIT_THERMAL || dom.dataInit == DATA_INIT_COLLISION || dom.dataInit == DATA_INIT_STRAKA) {
      t0 = 300._fp;
      hydro::hydroConstTheta( t0 , zloc , r0 );
    }

    dom.hyDensGLL     (k,kk) = r0;
    dom.hyDensThetaGLL(k,kk) = r0*t0;
    dom.hyThetaGLL    (k,kk) = t0;
    dom.hyPressureGLL (k,kk) = C0*pow( r0*t0 , GAMMA );
  });

  // Initialize the state
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( "InitFluidState" , dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    // Initialize the state to zero
    for (int l=0; l<numState; l++) {
      state(l,hs+k,hs+j,hs+i) = 0;
    }
    // Perform ord-point GLL quadrature for the cell averages
    for (int kk=0; kk<ord; kk++) {
      for (int jj=0; jj<ord; jj++) {
        for (int ii=0; ii<ord; ii++) {
          real xloc = (par.i_beg + i + 0.5_fp)*dom.dx + gllOrdPoints(ii)*dom.dx;
          real yloc = (par.j_beg + j + 0.5_fp)*dom.dy + gllOrdPoints(jj)*dom.dy;
          real zloc = (k + 0.5_fp)*dom.dz + gllOrdPoints(kk)*dom.dz;
          real r0, t0, r, t;

          if (dom.run2d) yloc = dom.ylen/2;

          if (dom.dataInit == DATA_INIT_THERMAL || dom.dataInit == DATA_INIT_COLLISION || dom.dataInit == DATA_INIT_STRAKA) {
            t0 = 300._fp;
            hydro::hydroConstTheta( t0 , zloc , r0 );
          }

          if        (dom.dataInit == DATA_INIT_COLLISION) {
            t  = ellipsoid_linear(xloc, yloc, zloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 2000, 2000,  20);
            t += ellipsoid_linear(xloc, yloc, zloc, dom.xlen/2, dom.ylen/2, 8000, 2000, 2000, 2000, -20);
          } else if (dom.dataInit == DATA_INIT_THERMAL  ) {
            t  = ellipsoid_linear(xloc, yloc, zloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 2000, 2000, 2  );
          } else if (dom.dataInit == DATA_INIT_STRAKA   ) {
            t  = ellipsoid_cosine(xloc, yloc, zloc, dom.xlen/2, dom.ylen/2, 3000, 4000, 4000, 2000, -15 , 1 );
          }

          // Set the initial density such that pressure is constant (to get rid of distracting acoustic waves)
          // r = r0*(t0/(t0+t)-1._fp);
          r = 0;

          real wt = gllOrdWeights(ii)*gllOrdWeights(jj)*gllOrdWeights(kk);

          state(idR ,hs+k,hs+j,hs+i) += wt * r  ;
          state(idRT,hs+k,hs+j,hs+i) += wt * ( (r0+r)*(t0+t) - r0*t0 );
        }
      }
    }
  });

  realArr dt3d("dt3d",dom.nz,dom.ny,dom.nx);
  // Compute the time step based on the CFL value
  // for (int k=0; k<dom.nz; k++) {
  //   for (int j=0; j<dom.ny; j++) {
  //     for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( "Compute_dt3d" , dom.nz,dom.ny,dom.nx , YAKL_LAMBDA (int k, int j, int i) {
    // Grab state variables
    real r, u, v, w, t, p, cs;
    r = state(idR ,hs+k,hs+j,hs+i) + dom.hyDensCells(hs+k);
    u = state(idRU,hs+k,hs+j,hs+i);
    v = state(idRV,hs+k,hs+j,hs+i);
    w = state(idRW,hs+k,hs+j,hs+i);
    t = ( state(idRT,hs+k,hs+j,hs+i) + dom.hyDensThetaCells(hs+k) ) / r;
    p = C0 * pow( r*t , GAMMA );
    cs = sqrt( GAMMA * p / r );

    // Compute the max wave
    real maxWave = mymax( mymax( fabs(u) , fabs(v)) , fabs(w)) + cs;

    // Compute the time step
    real mindx = mymin( mymin(dom.dx,dom.dy) , dom.dz);
    dt3d(k,j,i) = dom.cfl * mindx / maxWave;
  });

  yakl::ParallelMin<real,yakl::memDevice> pmin( dom.nx*dom.ny*dom.nz );
  dom.dt = pmin( dt3d.data() );

  if (par.nranks > 1) {
    real dtloc = dom.dt;
    ierr = MPI_Allreduce(&dtloc, &dom.dt, 1, MPI_REAL , MPI_MIN, MPI_COMM_WORLD);
  }

  if (par.masterproc) {
    std::cout << "dx: " << dom.dx << "\n";
    std::cout << "dy: " << dom.dy << "\n";
    std::cout << "dz: " << dom.dz << "\n";
    std::cout << "dt: " << dom.dt << "\n";
  }

}


