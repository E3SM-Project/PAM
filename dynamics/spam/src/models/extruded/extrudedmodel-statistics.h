#pragma once

namespace pamc {

class ModelStats : public Stats {
  using VS = VariableSet;

public:
  real3d TEarr, KEarr, PEarr, IEarr, PENSarr, trimmed_density;
  real4d PVarr;

  void initialize(ModelParameters &params, Parallel &par,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom, Equations &equations) {
    Stats::initialize(params, par, primal_geom, dual_geom, equations);
    this->stats_arr[DENSSTAT].initialize("mass", VS::ndensity_prognostic,
                                         this->statsize, this->nens,
                                         this->masterproc);
    this->stats_arr[DENSMAXSTAT].initialize("densmax", VS::ndensity_prognostic,
                                            this->statsize, this->nens,
                                            this->masterproc);
    this->stats_arr[DENSMINSTAT].initialize("densmin", VS::ndensity_prognostic,
                                            this->statsize, this->nens,
                                            this->masterproc);
    this->stats_arr[ESTAT].initialize("energy", 4, this->statsize, this->nens,
                                      this->masterproc);
    this->stats_arr[PVSTAT].initialize("pv", ndims > 1 ? 3 : 1, this->statsize,
                                       this->nens, this->masterproc);
    this->stats_arr[PESTAT].initialize("pens", 1, this->statsize, this->nens,
                                       this->masterproc);

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    TEarr = real3d("TE", dual_topology.nl, dual_topology.n_cells_y,
                   dual_topology.n_cells_x);
    KEarr = real3d("KE", dual_topology.nl, dual_topology.n_cells_y,
                   dual_topology.n_cells_x);
    IEarr = real3d("IE", dual_topology.nl, dual_topology.n_cells_y,
                   dual_topology.n_cells_x);
    PEarr = real3d("PE", dual_topology.nl, dual_topology.n_cells_y,
                   dual_topology.n_cells_x);
    PVarr = real4d("PV", ndims > 1 ? 3 : 1, primal_topology.nl,
                   primal_topology.n_cells_y, primal_topology.n_cells_x);
    PENSarr = real3d("PENS", primal_topology.nl, primal_topology.n_cells_y,
                     primal_topology.n_cells_x);
    trimmed_density = real3d("trimmed_density", dual_topology.nl,
                             dual_topology.n_cells_y, dual_topology.n_cells_x);
  }

  void compute(FieldSet<nprognostic> &progvars, FieldSet<nconstant> &constvars,
               int tind) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    for (int n = 0; n < nens; n++) {

      SArray<real, 1, VS::ndensity_prognostic> masslocal, massglobal;
      SArray<real, 1, VS::ndensity_prognostic> densmaxlocal, densmaxglobal;
      SArray<real, 1, VS::ndensity_prognostic> densminlocal, densminglobal;
      SArray<real, 1, (ndims > 1 ? 3 : 1)> pvlocal, pvglobal;
      SArray<real, 1, 4> elocal, eglobal;
      SArray<real, 1, 1> pelocal, peglobal;

      pvlocal(0) = 0.;
      pvglobal(0) = 0.;
      pelocal(0) = 0.;
      peglobal(0) = 0.;
      elocal(0) = 0.;
      elocal(1) = 0.;
      elocal(2) = 0.;
      elocal(3) = 0.;
      eglobal(0) = 0.;
      eglobal(1) = 0.;
      eglobal(2) = 0.;
      eglobal(3) = 0.;
      for (int l = 0; l < VS::ndensity_prognostic; l++) {
        masslocal(l) = 0.;
        massglobal(l) = 0.;
      }
      for (int l = 0; l < VS::ndensity_prognostic; l++) {
        densmaxlocal(l) = 0.;
        densmaxglobal(l) = 0.;
      }
      for (int l = 0; l < VS::ndensity_prognostic; l++) {
        densminlocal(l) = 0.;
        densminglobal(l) = 0.;
      }

      int dis = dual_topology.is;
      int djs = dual_topology.js;
      int dks = dual_topology.ks;

      YAKL_SCOPE(Hk, equations->Hk);
      YAKL_SCOPE(Hs, equations->Hs);
      parallel_for(
          "Compute energetics stats",
          SimpleBounds<3>(dual_topology.nl - 2, dual_topology.n_cells_y,
                          dual_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int k, int j, int i) {
            real KE, PE, IE;
            KE = Hk.compute_KE(progvars.fields_arr[VVAR].data,
                               progvars.fields_arr[WVAR].data,
                               progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, k + 1, n);
            PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data,
                               constvars.fields_arr[HSVAR].data, dis, djs, dks,
                               i, j, k + 1, n);
            IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, k + 1, n);
            TEarr(k + 1, j, i) = KE + PE + IE;
            KEarr(k + 1, j, i) = KE;
            PEarr(k + 1, j, i) = PE;
            IEarr(k + 1, j, i) = IE;
          });
      parallel_for(
          "Compute energetics stats bnd",
          SimpleBounds<2>(dual_topology.n_cells_y, dual_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int j, int i) {
            real KE, PE, IE;
            KE = Hk.compute_KE_bottom(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data, dis, djs, dks, i, j, 0, n);
            PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data,
                               constvars.fields_arr[HSVAR].data, dis, djs, dks,
                               i, j, 0, n);
            IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, 0, n);
            TEarr(0, j, i) = KE + PE + IE;
            KEarr(0, j, i) = KE;
            PEarr(0, j, i) = PE;
            IEarr(0, j, i) = IE;
            KE = Hk.compute_KE_top(progvars.fields_arr[VVAR].data,
                                   progvars.fields_arr[WVAR].data,
                                   progvars.fields_arr[DENSVAR].data, dis, djs,
                                   dks, i, j, dual_topology.nl - 1, n);
            PE = Hs.compute_PE(progvars.fields_arr[DENSVAR].data,
                               constvars.fields_arr[HSVAR].data, dis, djs, dks,
                               i, j, dual_topology.nl - 1, n);
            IE = Hs.compute_IE(progvars.fields_arr[DENSVAR].data, dis, djs, dks,
                               i, j, dual_topology.nl - 1, n);
            TEarr(dual_topology.nl - 1, j, i) = KE + PE + IE;
            KEarr(dual_topology.nl - 1, j, i) = KE;
            PEarr(dual_topology.nl - 1, j, i) = PE;
            IEarr(dual_topology.nl - 1, j, i) = IE;
          });

      elocal(0) = yakl::intrinsics::sum(TEarr);
      elocal(1) = yakl::intrinsics::sum(KEarr);
      elocal(2) = yakl::intrinsics::sum(PEarr);
      elocal(3) = yakl::intrinsics::sum(IEarr);

      int pis = primal_topology.is;
      int pjs = primal_topology.js;
      int pks = primal_topology.ks;

      YAKL_SCOPE(PVPE, equations->PVPE);
      parallel_for(
          "Compute PV/PE stats",
          SimpleBounds<3>(primal_topology.nl - 2, primal_topology.n_cells_y,
                          primal_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int k, int j, int i) {
            pvpe_extruded vals_pvpe;
            vals_pvpe = PVPE.compute_PVPE(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data,
                constvars.fields_arr[CORIOLISHZVAR].data,
                constvars.fields_arr[CORIOLISXYVAR].data, pis, pjs, pks, i, j,
                k + 1, n);

            for (int d = 0; d < (ndims > 1 ? 3 : 1); ++d) {
              PVarr(d, k + 1, j, i) = vals_pvpe.pv(d);
            }
            PENSarr(k + 1, j, i) = vals_pvpe.pe;
          });
      parallel_for(
          "Compute PV/PE stats bnd",
          SimpleBounds<2>(primal_topology.n_cells_y, primal_topology.n_cells_x),
          YAKL_CLASS_LAMBDA(int j, int i) {
            pvpe_extruded vals_pvpe;
            vals_pvpe = PVPE.compute_PVPE_bottom(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data,
                constvars.fields_arr[CORIOLISHZVAR].data,
                constvars.fields_arr[CORIOLISXYVAR].data, pis, pjs, pks, i, j,
                0, n);
            for (int d = 0; d < (ndims > 1 ? 3 : 1); ++d) {
              PVarr(d, 0, j, i) = vals_pvpe.pv(d);
            }
            PENSarr(0, j, i) = vals_pvpe.pe;
            vals_pvpe = PVPE.compute_PVPE_top(
                progvars.fields_arr[VVAR].data, progvars.fields_arr[WVAR].data,
                progvars.fields_arr[DENSVAR].data,
                constvars.fields_arr[CORIOLISHZVAR].data,
                constvars.fields_arr[CORIOLISXYVAR].data, pis, pjs, pks, i, j,
                primal_topology.nl - 1, n);
            for (int d = 0; d < (ndims > 1 ? 3 : 1); ++d) {
              PVarr(d, primal_topology.nl - 1, j, i) = vals_pvpe.pv(d);
            }
            PENSarr(primal_topology.nl - 1, j, i) = vals_pvpe.pe;
          });

      pvlocal(0) = yakl::intrinsics::sum(
          PVarr.slice<3>(0, yakl::COLON, yakl::COLON, yakl::COLON));
      if (ndims > 1) {
        pvlocal(1) = yakl::intrinsics::sum(
            PVarr.slice<3>(1, yakl::COLON, yakl::COLON, yakl::COLON));
        pvlocal(2) = yakl::intrinsics::sum(
            PVarr.slice<3>(2, yakl::COLON, yakl::COLON, yakl::COLON));
      }
      pelocal(0) = yakl::intrinsics::sum(PENSarr);

      for (int l = 0; l < VS::ndensity_prognostic; l++) {
        parallel_for(
            "Compute trimmed density",
            SimpleBounds<3>(dual_topology.nl, dual_topology.n_cells_y,
                            dual_topology.n_cells_x),
            YAKL_CLASS_LAMBDA(int k, int j, int i) {
              trimmed_density(k, j, i) = progvars.fields_arr[DENSVAR].data(
                  l, k + dks, j + djs, i + dis, n);
            });

        masslocal(l) = yakl::intrinsics::sum(trimmed_density);
        densmaxlocal(l) = yakl::intrinsics::maxval(trimmed_density);
        densminlocal(l) = yakl::intrinsics::minval(trimmed_density);
      }

      // MPI sum/min/max
      this->ierr = MPI_Ireduce(&masslocal, &massglobal, VS::ndensity_prognostic,
                               PAMC_MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD,
                               &this->Req[DENSSTAT]);
      this->ierr = MPI_Ireduce(&densmaxlocal, &densmaxglobal,
                               VS::ndensity_prognostic, PAMC_MPI_REAL, MPI_MAX,
                               0, MPI_COMM_WORLD, &this->Req[DENSMAXSTAT]);
      this->ierr = MPI_Ireduce(&densminlocal, &densminglobal,
                               VS::ndensity_prognostic, PAMC_MPI_REAL, MPI_MIN,
                               0, MPI_COMM_WORLD, &this->Req[DENSMINSTAT]);
      this->ierr =
          MPI_Ireduce(&pvlocal, &pvglobal, ndims > 1 ? 3 : 1, PAMC_MPI_REAL,
                      MPI_SUM, 0, MPI_COMM_WORLD, &this->Req[PVSTAT]);
      this->ierr = MPI_Ireduce(&pelocal, &peglobal, 1, PAMC_MPI_REAL, MPI_SUM,
                               0, MPI_COMM_WORLD, &this->Req[PESTAT]);
      this->ierr = MPI_Ireduce(&elocal, &eglobal, 4, PAMC_MPI_REAL, MPI_SUM, 0,
                               MPI_COMM_WORLD, &this->Req[ESTAT]);

      this->ierr = MPI_Waitall(nstats, this->Req, this->Status);

      if (masterproc) {
        for (int l = 0; l < VS::ndensity_prognostic; l++) {
          this->stats_arr[DENSSTAT].data(l, tind, n) = massglobal(l);
          this->stats_arr[DENSMAXSTAT].data(l, tind, n) = densmaxglobal(l);
          this->stats_arr[DENSMINSTAT].data(l, tind, n) = densminglobal(l);
        }

        this->stats_arr[ESTAT].data(0, tind, n) = eglobal(0);
        this->stats_arr[ESTAT].data(1, tind, n) = eglobal(1);
        this->stats_arr[ESTAT].data(2, tind, n) = eglobal(2);
        this->stats_arr[ESTAT].data(3, tind, n) = eglobal(3);
        for (int d = 0; d < (ndims > 1 ? 3 : 1); ++d) {
          this->stats_arr[PVSTAT].data(d, tind, n) = pvglobal(d);
        }
        this->stats_arr[PESTAT].data(0, tind, n) = peglobal(0);
      }
    }
  }
};

} // namespace pamc
