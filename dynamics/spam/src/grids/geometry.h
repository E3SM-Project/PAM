#pragma once

#include "common.h"
#include "fields.h"
#include "params.h"
#include "profiles.h"
#include "topology.h"

enum class LINE_INTEGRAL_TYPE { TANGENT, NORMAL };

struct CoordsXYZ {
  real x = 0.0_fp;
  real y = 0.0_fp;
  real z = 0.0_fp;
};
struct VecXY {
  real u = 0, v = 0;
};
struct VecXYZ {
  real u = 0, v = 0, w = 0;
};

template <class T>
void set_ref_quad_pts_wts(std::vector<T> &pts, std::vector<T> &wts, int npts) {
  if (npts == 1) {
    pts[0] = 0.5_fp;
    wts[0] = 1.0_fp;
  }

  if (npts == 2) {
    pts[0] = -1.0_fp / (2.0_fp * sqrt(3.)) + 0.5_fp;
    pts[1] = 1.0_fp / (2.0_fp * sqrt(3.)) + 0.5_fp;

    wts[0] = 0.5_fp;
    wts[1] = 0.5_fp;
  }

  if (npts == 3) {
    pts[0] = -0.5_fp * sqrt(3.0_fp / 5.0_fp) + 0.5_fp;
    pts[1] = 0.5_fp;
    pts[2] = 0.5_fp * sqrt(3.0_fp / 5.0_fp) + 0.5_fp;

    wts[0] = 5.0_fp / 18.0_fp;
    wts[1] = 4.0_fp / 9.0_fp;
    wts[2] = 5.0_fp / 18.0_fp;
  }
  //   {
  //     pts(0) = 0.112701665379258311482073460022;
  //     pts(1) = 0.500000000000000000000000000000;
  //     pts(2) = 0.887298334620741688517926539980;
  //
  //     wts(0) = 0.277777777777777777777777777779;
  //     wts(1) = 0.444444444444444444444444444444;
  //     wts(2) = 0.277777777777777777777777777779;
  //   }

  if (npts == 4) {
    pts[0] = -0.5_fp * sqrt(3.0_fp / 7.0_fp +
                            2.0_fp / 7.0_fp * sqrt(6.0_fp / 5.0_fp)) +
             0.5_fp;
    pts[1] = -0.5_fp * sqrt(3.0_fp / 7.0_fp -
                            2.0_fp / 7.0_fp * sqrt(6.0_fp / 5.0_fp)) +
             0.5_fp;
    pts[2] = 0.5_fp * sqrt(3.0_fp / 7.0_fp -
                           2.0_fp / 7.0_fp * sqrt(6.0_fp / 5.0_fp)) +
             0.5_fp;
    pts[3] = 0.5_fp * sqrt(3.0_fp / 7.0_fp +
                           2.0_fp / 7.0_fp * sqrt(6.0_fp / 5.0_fp)) +
             0.5_fp;

    wts[0] = (18.0_fp - sqrt(30.0_fp)) / 72.0_fp;
    wts[1] = (18.0_fp + sqrt(30.0_fp)) / 72.0_fp;
    wts[2] = (18.0_fp + sqrt(30.0_fp)) / 72.0_fp;
    wts[3] = (18.0_fp - sqrt(30.0_fp)) / 72.0_fp;
  }

  if (npts == 5) {
    pts[0] = -1.0_fp / 6.0_fp * sqrt(5.0_fp + 2.0_fp * sqrt(10.0_fp / 7.0_fp)) +
             0.5_fp;
    pts[1] = -1.0_fp / 6.0_fp * sqrt(5.0_fp - 2.0_fp * sqrt(10.0_fp / 7.0_fp)) +
             0.5_fp;
    pts[2] = 0.5_fp;
    pts[3] = 1.0_fp / 6.0_fp * sqrt(5.0_fp - 2.0_fp * sqrt(10.0_fp / 7.0_fp)) +
             0.5_fp;
    pts[4] = 1.0_fp / 6.0_fp * sqrt(5.0_fp + 2.0_fp * sqrt(10.0_fp / 7.0_fp)) +
             0.5_fp;

    wts[0] = (322.0_fp - 13.0_fp * sqrt(70.0_fp)) / 1800.0_fp;
    wts[1] = (322.0_fp + 13.0_fp * sqrt(70.0_fp)) / 1800.0_fp;
    wts[2] = 64.0_fp / 225.0_fp;
    wts[3] = (322.0_fp + 13.0_fp * sqrt(70.0_fp)) / 1800.0_fp;
    wts[4] = (322.0_fp - 13.0_fp * sqrt(70.0_fp)) / 1800.0_fp;
  }
}

template <class T, uint npts>
void set_ref_quad_pts_wts(SArray<T, 1, npts> &pts, SArray<T, 1, npts> &wts) {
  std::vector<T> pts_v(npts);
  std::vector<T> wts_v(npts);
  set_ref_quad_pts_wts(pts_v, wts_v, npts);
  for (int m = 0; m < npts; ++m) {
    pts(m) = pts_v[m];
    wts(m) = wts_v[m];
  }
}

struct Straight {};
struct Twisted {};

template <class T> class Geometry {

public:
  Topology topology;
  bool is_initialized;
  static constexpr bool straight = std::is_same_v<T, Straight>;
  bool uniform_vertical;

  Geometry();

  SArray<real, 1, ic_quad_pts_x> x_quad_pts_ref;
  SArray<real, 1, ic_quad_pts_x> x_quad_wts_ref;

  SArray<real, 1, ic_quad_pts_y> y_quad_pts_ref;
  SArray<real, 1, ic_quad_pts_y> y_quad_wts_ref;

  SArray<real, 1, ic_quad_pts_z> z_quad_pts_ref;
  SArray<real, 1, ic_quad_pts_z> z_quad_wts_ref;

  void initialize(const Topology &topo, const ModelParameters &params);
  void printinfo() const;

  void YAKL_INLINE get_ll_corner(CoordsXYZ &llc, int k, int j, int i,
                                 int n) const;
  real YAKL_INLINE get_zint(int k, int n) const;
  real YAKL_INLINE get_dz(int k, int n) const;

  void YAKL_INLINE get_11face_normals(
      int k, int j, int i, int n,
      SArray<VecXYZ, 2, ic_quad_pts_x, ic_quad_pts_z> &xz_normals,
      SArray<VecXYZ, 2, ic_quad_pts_y, ic_quad_pts_z> &yz_normals) const;

  void YAKL_INLINE get_20face_normals(
      int k, int j, int i, int n,
      SArray<VecXYZ, 2, ic_quad_pts_x, ic_quad_pts_y> &xy_normals) const;

  void YAKL_INLINE get_10edge_tangents(
      int k, int j, int i, int n, SArray<VecXYZ, 1, ic_quad_pts_x> &x_tangent,
      SArray<VecXYZ, 1, ic_quad_pts_y> &y_tangent) const;

  void YAKL_INLINE
  get_01edge_tangents(int k, int j, int i, int n,
                      SArray<VecXYZ, 1, ic_quad_pts_z> &z_tangent) const;

  void YAKL_INLINE get_00form_quad_pts_wts(
      int k, int j, int i, int n, SArray<CoordsXYZ, 1, 1> &quad_pts_phys,
      SArray<real, 1, 1> &quad_wts_phys) const;

  void YAKL_INLINE get_10form_quad_pts_wts(
      int k, int j, int i, int n,
      SArray<CoordsXYZ, 1, ic_quad_pts_x> &x_quad_pts_phys,
      SArray<real, 1, ic_quad_pts_x> &x_quad_wts_phys,
      SArray<CoordsXYZ, 1, ic_quad_pts_y> &y_quad_pts_phys,
      SArray<real, 1, ic_quad_pts_y> &y_quad_wts_phys) const;

  void YAKL_INLINE get_01form_quad_pts_wts(
      int k, int j, int i, int n,
      SArray<CoordsXYZ, 1, ic_quad_pts_z> &z_quad_pts_phys,
      SArray<real, 1, ic_quad_pts_z> &z_quad_wts_phys) const;

  void YAKL_INLINE get_11form_quad_pts_wts(
      int k, int j, int i, int n,
      SArray<CoordsXYZ, 2, ic_quad_pts_x, ic_quad_pts_z> &xz_quad_pts_phys,
      SArray<real, 2, ic_quad_pts_x, ic_quad_pts_z> &xz_quad_wts_phys,
      SArray<CoordsXYZ, 2, ic_quad_pts_y, ic_quad_pts_z> &yz_quad_pts_phys,
      SArray<real, 2, ic_quad_pts_y, ic_quad_pts_z> &yz_quad_wts_phys) const;

  void YAKL_INLINE get_20form_quad_pts_wts(
      int k, int j, int i, int n,
      SArray<CoordsXYZ, 2, ic_quad_pts_x, ic_quad_pts_y> &xy_quad_pts_phys,
      SArray<real, 2, ic_quad_pts_x, ic_quad_pts_y> &xy_quad_wts_phys) const;

  void YAKL_INLINE get_21form_quad_pts_wts(
      int k, int j, int i, int n,
      SArray<CoordsXYZ, 3, ic_quad_pts_x, ic_quad_pts_y, ic_quad_pts_z>
          &quad_pts_phys,
      SArray<real, 3, ic_quad_pts_x, ic_quad_pts_y, ic_quad_pts_z>
          &quad_wts_phys) const;

  template <int l>
  real YAKL_INLINE get_area_lform(int d, int k, int j, int i, int n) const;

  real YAKL_INLINE get_area_00entity(int k, int j, int i, int n) const;
  real YAKL_INLINE get_area_10entity(int d, int k, int j, int i, int n) const;
  real YAKL_INLINE get_area_20entity(int k, int j, int i, int n) const;
  real YAKL_INLINE get_area_01entity(int k, int j, int i, int n) const;
  real YAKL_INLINE get_area_11entity(int d, int k, int j, int i, int n) const;
  real YAKL_INLINE get_area_21entity(int k, int j, int i, int n) const;

  real YAKL_INLINE get_area_n1entity(int k, int j, int i, int n) const;
  real YAKL_INLINE get_area_n0entity(int k, int j, int i, int n) const;
  real YAKL_INLINE get_area_nm11entity(int d, int k, int j, int i, int n) const;

  template <class F>
  void set_00form_values(F initial_value_function, Field &field,
                         int ndof) const;
  template <class F>
  void set_00form_values(F initial_value_function, Field &field, int ndof,
                         int nz) const;
  template <class F>
  void set_10form_values(F initial_value_function, Field &field,
                         int ndof) const;
  template <class F>
  void set_10form_values(F initial_value_function, Field &field, int ndof,
                         int nz) const;

  template <class F>
  void set_01form_values(F initial_value_function, Field &field,
                         int ndof) const;

  template <class F>
  void set_nm11form_values(F initial_value_function, Field &field,
                           int ndof) const;

  template <class F>
  void set_n0form_values(F initial_value_function, Field &field,
                         int ndof) const;

  template <class F>
  void set_n1form_values(F initial_value_function, Field &field,
                         int ndof) const;

  template <class F>
  void set_profile_00form_values(F initial_value_function, Profile &prof,
                                 int ndof) const;
  template <class F>
  void set_profile_n1form_values(F initial_value_function, Profile &prof,
                                 int ndof) const;

  template <class F>
  void set_0form_values(F initial_value_function, Field &field, int ndof) const;
  template <class F>
  void set_1form_values(F initial_value_function, Field &field, int ndof,
                        LINE_INTEGRAL_TYPE line_type) const;
  template <class F>
  void set_2form_values(F initial_value_function, Field &field, int ndof) const;

  real dx, dy;
  real Lx, Ly;
  real xc, yc;

  real2d dz;
  real2d zint;
};

template <class T> Geometry<T>::Geometry() { this->is_initialized = false; }

template <class T>
void Geometry<T>::initialize(const Topology &topo,
                             const ModelParameters &params) {

  this->topology = topo;

  set_ref_quad_pts_wts<real, ic_quad_pts_x>(this->x_quad_pts_ref,
                                            this->x_quad_wts_ref);
  set_ref_quad_pts_wts<real, ic_quad_pts_y>(this->y_quad_pts_ref,
                                            this->y_quad_wts_ref);
  set_ref_quad_pts_wts<real, ic_quad_pts_z>(this->z_quad_pts_ref,
                                            this->z_quad_wts_ref);

  this->Lx = params.xlen;
  this->Ly = params.ylen;

  this->xc = params.xc;
  this->yc = params.yc;

  this->dx = params.xlen / params.nx_glob;
  if (ndims > 1) {
    this->dy = params.ylen / params.ny_glob;
  } else {
    this->dy = 1;
  }

#ifdef PAMC_EXTRUDED

  if (straight) {
    this->dz = real2d("dz straight", topo.nl + 2 * topo.mirror_halo, topo.nens);
    this->zint =
        real2d("zint straight", topo.ni + 2 * topo.mirror_halo, topo.nens);

    int ks = topo.ks;

    YAKL_SCOPE(zint, this->zint);
    YAKL_SCOPE(dz, this->dz);
    // This code puts straight grid interfaces at the midpoint point between
    // twisted grid interfaces (other than the top and bottom levels)
    parallel_for(
        "Set zint straight", SimpleBounds<2>(topo.ni, topo.nens),
        YAKL_LAMBDA(int k, int n) {
          if (k == 0) {
            zint(k + topo.ks, n) = params.zint(k, n);
          } else if (k == topo.ni - 1) {
            zint(k + topo.ks, n) = params.zint(k + 1, n);
          } // need to add 1 here due to indexing between straight and twisted
            // grids
          else {
            zint(k + topo.ks, n) =
                (params.zint(k, n) + params.zint(k + 1, n)) / 2.0_fp;
          }
        });

    parallel_for(
        "Set zint straight halo", SimpleBounds<2>(topo.mirror_halo, topo.nens),
        YAKL_LAMBDA(int k, int n) {
          zint(-(k + 1) + ks, n) =
              zint(ks, n) - (zint(k + 1 + ks, n) - zint(ks, n));
          zint(k + ks + topo.ni, n) =
              zint(ks + topo.ni - 1, n) +
              (zint(ks + topo.ni - 1, n) - zint(-k + ks + topo.ni - 2, n));
        });

    parallel_for(
        "Set dz straight",
        SimpleBounds<2>(topo.nl + 2 * topo.mirror_halo, topo.nens),
        YAKL_LAMBDA(int k, int n) { dz(k, n) = zint(k + 1, n) - zint(k, n); });
  } else {
    this->dz = real2d("dz twisted", topo.nl + 2 * topo.mirror_halo, topo.nens);
    this->zint =
        real2d("zint twisted", topo.ni + 2 * topo.mirror_halo, topo.nens);

    int ks = topo.ks;

    YAKL_SCOPE(zint, this->zint);
    YAKL_SCOPE(dz, this->dz);
    parallel_for(
        "Set zint twisted", SimpleBounds<2>(topo.ni, topo.nens),
        YAKL_LAMBDA(int k, int n) {
          zint(k + topo.ks, n) = params.zint(k, n);
        });

    parallel_for(
        "Set zint twisted halo", SimpleBounds<2>(topo.mirror_halo, topo.nens),
        YAKL_LAMBDA(int k, int n) {
          zint(-(k + 1) + ks, n) =
              zint(ks, n) - (zint(k + 1 + ks, n) - zint(ks, n));
          zint(k + ks + topo.ni, n) =
              zint(ks + topo.ni - 1, n) +
              (zint(ks + topo.ni - 1, n) - zint(-k + ks + topo.ni - 2, n));
        });

    parallel_for(
        "Set dz twisted",
        SimpleBounds<2>(topo.nl + 2 * topo.mirror_halo, topo.nens),
        YAKL_LAMBDA(int k, int n) { dz(k, n) = zint(k + 1, n) - zint(k, n); });
  }

  this->uniform_vertical = params.uniform_vertical;
#endif

  this->is_initialized = true;
}

template <class T>
void YAKL_INLINE Geometry<T>::get_ll_corner(CoordsXYZ &llc, int k, int j, int i,
                                            int n) const {
  if (straight) {
    llc.x = (i + this->topology.i_beg) * this->dx + this->xc -
            this->Lx * 0.5_fp - this->dx * 0.5_fp;
    llc.y = (j + this->topology.j_beg) * this->dy + this->yc -
            this->Ly * 0.5_fp - this->dy * 0.5_fp;
  } else {
    llc.x =
        (i + this->topology.i_beg) * this->dx + this->xc - this->Lx * 0.5_fp;
    llc.y =
        (j + this->topology.j_beg) * this->dy + this->yc - this->Ly * 0.5_fp;
  }
  llc.z = this->get_zint(k + this->topology.ks, n);
}

template <class T> real YAKL_INLINE Geometry<T>::get_zint(int k, int n) const {
#ifdef PAMC_EXTRUDED
  return this->zint(k, n);
#else
  return 0;
#endif
}

template <class T> real YAKL_INLINE Geometry<T>::get_dz(int k, int n) const {
#ifdef PAMC_EXTRUDED
  return this->dz(k, n);
#else
  return 1;
#endif
}

template <class T>
real YAKL_INLINE Geometry<T>::get_area_00entity(int k, int j, int i,
                                                int n) const {
  return 1._fp;
}

template <class T>
real YAKL_INLINE Geometry<T>::get_area_01entity(int k, int j, int i,
                                                int n) const {
  return this->get_dz(k, n);
}

template <class T>
real YAKL_INLINE Geometry<T>::get_area_10entity(int d, int k, int j, int i,
                                                int n) const {
  if (ndims == 1) {
    return this->dx;
  } else {
    if (straight) {
      return (d == 0 ? this->dx : this->dy);
    } else {
      return (d == 0 ? this->dy : this->dx);
    }
  }
}

template <class T>
real YAKL_INLINE Geometry<T>::get_area_11entity(int d, int k, int j, int i,
                                                int n) const {
  return get_area_10entity(d, k, j, i, n) * this->get_dz(k, n);
}

template <class T>
real YAKL_INLINE Geometry<T>::get_area_20entity(int k, int j, int i,
                                                int n) const {
  return this->dx * this->dy;
}

template <class T>
real YAKL_INLINE Geometry<T>::get_area_21entity(int k, int j, int i,
                                                int n) const {
  return this->dx * this->dy * this->get_dz(k, n);
}

template <class T>
real YAKL_INLINE Geometry<T>::get_area_n1entity(int k, int j, int i,
                                                int n) const {
  return (ndims == 1) ? get_area_11entity(0, k, j, i, n)
                      : get_area_21entity(k, j, i, n);
}

template <class T>
real YAKL_INLINE Geometry<T>::get_area_n0entity(int k, int j, int i,
                                                int n) const {
  return (ndims == 1) ? get_area_10entity(0, k, j, i, n)
                      : get_area_20entity(k, j, i, n);
}

template <class T>
real YAKL_INLINE Geometry<T>::get_area_nm11entity(int d, int k, int j, int i,
                                                  int n) const {
  return (ndims == 1) ? get_area_01entity(k, j, i, n)
                      : get_area_11entity(d, k, j, i, n);
}

template <class T>
template <int l>
real YAKL_INLINE Geometry<T>::get_area_lform(int d, int k, int j, int i,
                                             int n) const {
  if constexpr (l == 0) {
    return get_area_00entity(k, j, i, n);
  }
  if constexpr (l == 1) {
    return get_area_10entity(d, k, j, i, n);
  }
  if constexpr (l == 2) {
    return get_area_20entity(k, j, i, n);
  }
}

template <class T>
void YAKL_INLINE Geometry<T>::get_11face_normals(
    int k, int j, int i, int n,
    SArray<VecXYZ, 2, ic_quad_pts_x, ic_quad_pts_z> &xz_normals,
    SArray<VecXYZ, 2, ic_quad_pts_y, ic_quad_pts_z> &yz_normals) const {

  if (straight) {
    for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
      for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
        xz_normals(nqx, nqz).u = 0;
        xz_normals(nqx, nqz).v = 1;
        xz_normals(nqx, nqz).w = 0;
      }
    }

    for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
      for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
        yz_normals(nqy, nqz).u = 1;
        yz_normals(nqy, nqz).v = 0;
        yz_normals(nqy, nqz).w = 0;
      }
    }
  } else {
    for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
      for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
        xz_normals(nqx, nqz).u = 0;
        xz_normals(nqx, nqz).v = 1;
        xz_normals(nqx, nqz).w = 0;
      }
    }

    for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
      for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
        yz_normals(nqy, nqz).u = 1;
        yz_normals(nqy, nqz).v = 0;
        yz_normals(nqy, nqz).w = 0;
      }
    }
  }
}

template <class T>
void YAKL_INLINE Geometry<T>::get_20face_normals(
    int k, int j, int i, int n,
    SArray<VecXYZ, 2, ic_quad_pts_x, ic_quad_pts_y> &xy_normals) const {

  if (straight) {
    for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
      for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
        xy_normals(nqx, nqy).u = 0;
        xy_normals(nqx, nqy).v = 0;
        if (ndims == 1) {
          xy_normals(nqx, nqy).w = -1;
        } else {
          xy_normals(nqx, nqy).w = 1;
        }
      }
    }
  } else {
    for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
      for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
        xy_normals(nqx, nqy).u = 0;
        xy_normals(nqx, nqy).v = 0;
        xy_normals(nqx, nqy).w = 1;
      }
    }
  }
}

template <class T>
void YAKL_INLINE Geometry<T>::get_10edge_tangents(
    int k, int j, int i, int n, SArray<VecXYZ, 1, ic_quad_pts_x> &x_tangent,
    SArray<VecXYZ, 1, ic_quad_pts_y> &y_tangent) const {

  if (straight) {
    for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
      x_tangent(nqx).u = 1;
      x_tangent(nqx).v = 0;
      x_tangent(nqx).w = 0;
    }
    for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
      y_tangent(nqy).u = 0;
      y_tangent(nqy).v = 1;
      y_tangent(nqy).w = 0;
    }
  } else {
    for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
#ifdef PAMC_LAYER
      x_tangent(nqx).u = -1;
#else
      x_tangent(nqx).u = 1;
#endif
      x_tangent(nqx).v = 0;
      x_tangent(nqx).w = 0;
    }
    for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
      y_tangent(nqy).u = 0;
      y_tangent(nqy).v = 1;
      y_tangent(nqy).w = 0;
    }
  }
}

template <class T>
void YAKL_INLINE Geometry<T>::get_01edge_tangents(
    int k, int j, int i, int n,
    SArray<VecXYZ, 1, ic_quad_pts_z> &z_tangent) const {

  if (straight) {
    for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
      z_tangent(nqz).u = 0;
      z_tangent(nqz).v = 0;
      z_tangent(nqz).w = 1;
    }
  } else {
    for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
      z_tangent(nqz).u = 0;
      z_tangent(nqz).v = 0;
      if (ndims == 1) {
        z_tangent(nqz).w = -1;
      } else {
        z_tangent(nqz).w = 1;
      }
    }
  }
}

template <class T>
void YAKL_INLINE Geometry<T>::get_00form_quad_pts_wts(
    int k, int j, int i, int n, SArray<CoordsXYZ, 1, 1> &quad_pts_phys,
    SArray<real, 1, 1> &quad_wts_phys) const {

  CoordsXYZ ll_corner;
  get_ll_corner(ll_corner, k, j, i, n);

  if (straight) {
    quad_pts_phys(0).x = ll_corner.x + this->dx;
    quad_pts_phys(0).y = ll_corner.y + this->dy;
    quad_pts_phys(0).z = ll_corner.z;
  } else {
    quad_pts_phys(0).x = ll_corner.x;
    quad_pts_phys(0).y = ll_corner.y;
    quad_pts_phys(0).z = ll_corner.z;
  }

  quad_wts_phys(0) = 1;
}

template <class T>
void YAKL_INLINE Geometry<T>::get_10form_quad_pts_wts(
    int k, int j, int i, int n,
    SArray<CoordsXYZ, 1, ic_quad_pts_x> &x_quad_pts_phys,
    SArray<real, 1, ic_quad_pts_x> &x_quad_wts_phys,
    SArray<CoordsXYZ, 1, ic_quad_pts_y> &y_quad_pts_phys,
    SArray<real, 1, ic_quad_pts_y> &y_quad_wts_phys) const {
  int ks = this->topology.ks;

  CoordsXYZ ll_corner;
  get_ll_corner(ll_corner, k, j, i, n);

  CoordsXYZ rr_corner;
  rr_corner.x = ll_corner.x + dx;
  rr_corner.y = ll_corner.y + dy;
  rr_corner.z = ll_corner.z + get_dz(k + ks, n);

  // x
  for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
    x_quad_pts_phys(nqx).x = ll_corner.x + this->x_quad_pts_ref(nqx) * this->dx;
    x_quad_pts_phys(nqx).y = straight ? rr_corner.y : ll_corner.y;
    x_quad_pts_phys(nqx).z = ll_corner.z;

    x_quad_wts_phys(nqx) = this->x_quad_wts_ref(nqx) * this->dx;
  }

  // y
  if (ndims > 1) {
    for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
      y_quad_pts_phys(nqy).x = straight ? rr_corner.x : ll_corner.x;
      y_quad_pts_phys(nqy).y =
          ll_corner.y + this->y_quad_pts_ref(nqy) * this->dy;
      y_quad_pts_phys(nqy).z = ll_corner.z;
      y_quad_wts_phys(nqy) = this->y_quad_wts_ref(nqy) * this->dy;
    }
  }
}

template <class T>
void YAKL_INLINE Geometry<T>::get_01form_quad_pts_wts(
    int k, int j, int i, int n,
    SArray<CoordsXYZ, 1, ic_quad_pts_z> &z_quad_pts_phys,
    SArray<real, 1, ic_quad_pts_z> &z_quad_wts_phys) const {
  int ks = this->topology.ks;

  CoordsXYZ ll_corner;
  get_ll_corner(ll_corner, k, j, i, n);

  for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
    z_quad_pts_phys(nqz).x = straight ? ll_corner.x + this->dx : ll_corner.x;
    z_quad_pts_phys(nqz).y = straight ? ll_corner.y + this->dy : ll_corner.y;
    z_quad_pts_phys(nqz).z =
        ll_corner.z + this->z_quad_pts_ref(nqz) * get_dz(k + ks, n);
    z_quad_wts_phys(nqz) = this->z_quad_wts_ref(nqz) * get_dz(k + ks, n);
  }
}

template <class T>
void YAKL_INLINE Geometry<T>::get_11form_quad_pts_wts(
    int k, int j, int i, int n,
    SArray<CoordsXYZ, 2, ic_quad_pts_x, ic_quad_pts_z> &xz_quad_pts_phys,
    SArray<real, 2, ic_quad_pts_x, ic_quad_pts_z> &xz_quad_wts_phys,
    SArray<CoordsXYZ, 2, ic_quad_pts_y, ic_quad_pts_z> &yz_quad_pts_phys,
    SArray<real, 2, ic_quad_pts_y, ic_quad_pts_z> &yz_quad_wts_phys) const {

  int ks = this->topology.ks;

  CoordsXYZ ll_corner;
  get_ll_corner(ll_corner, k, j, i, n);

  CoordsXYZ rr_corner;
  rr_corner.x = ll_corner.x + dx;
  rr_corner.y = ll_corner.y + dy;
  rr_corner.z = ll_corner.z + get_dz(k + ks, n);

  // xz
  for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
    for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
      xz_quad_pts_phys(nqx, nqz).x =
          ll_corner.x + this->x_quad_pts_ref(nqx) * this->dx;
      xz_quad_pts_phys(nqx, nqz).y = straight ? rr_corner.y : ll_corner.y;
      xz_quad_pts_phys(nqx, nqz).z =
          ll_corner.z + this->z_quad_pts_ref(nqz) * this->get_dz(k + ks, n);
      xz_quad_wts_phys(nqx, nqz) = this->x_quad_wts_ref(nqx) * this->dx *
                                   this->z_quad_wts_ref(nqz) *
                                   this->get_dz(k + ks, n);
    }
  }

  // yz
  for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
    for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
      yz_quad_pts_phys(nqy, nqz).x = straight ? rr_corner.x : ll_corner.x;
      yz_quad_pts_phys(nqy, nqz).y =
          ll_corner.y + this->y_quad_pts_ref(nqy) * this->dy;
      yz_quad_pts_phys(nqy, nqz).z =
          ll_corner.z + this->z_quad_pts_ref(nqz) * this->get_dz(k + ks, n);
      yz_quad_wts_phys(nqy, nqz) = this->y_quad_wts_ref(nqy) * this->dy *
                                   this->z_quad_wts_ref(nqz) *
                                   this->get_dz(k + ks, n);
    }
  }
}

template <class T>
void YAKL_INLINE Geometry<T>::get_20form_quad_pts_wts(
    int k, int j, int i, int n,
    SArray<CoordsXYZ, 2, ic_quad_pts_x, ic_quad_pts_y> &xy_quad_pts_phys,
    SArray<real, 2, ic_quad_pts_x, ic_quad_pts_y> &xy_quad_wts_phys) const {

  CoordsXYZ ll_corner;
  get_ll_corner(ll_corner, k, j, i, n);

  // xy
  for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
    for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
      xy_quad_pts_phys(nqx, nqy).x =
          ll_corner.x + this->x_quad_pts_ref(nqx) * this->dx;
      xy_quad_pts_phys(nqx, nqy).y =
          ll_corner.y + this->y_quad_pts_ref(nqy) * this->dy;
      xy_quad_pts_phys(nqx, nqy).z = ll_corner.z;
      xy_quad_wts_phys(nqx, nqy) = this->x_quad_wts_ref(nqx) * this->dx *
                                   this->y_quad_wts_ref(nqy) * this->dy;
    }
  }
}

template <class T>
void YAKL_INLINE Geometry<T>::get_21form_quad_pts_wts(
    int k, int j, int i, int n,
    SArray<CoordsXYZ, 3, ic_quad_pts_x, ic_quad_pts_y, ic_quad_pts_z>
        &quad_pts_phys,
    SArray<real, 3, ic_quad_pts_x, ic_quad_pts_y, ic_quad_pts_z> &quad_wts_phys)
    const {

  int ks = this->topology.ks;

  CoordsXYZ ll_corner;
  get_ll_corner(ll_corner, k, j, i, n);

  for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
    for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
      for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
        quad_pts_phys(nqx, nqy, nqz).x =
            ll_corner.x + this->x_quad_pts_ref(nqx) * this->dx;
        quad_pts_phys(nqx, nqy, nqz).y =
            ll_corner.y + this->y_quad_pts_ref(nqy) * this->dy;
        quad_pts_phys(nqx, nqy, nqz).z =
            ll_corner.z + this->z_quad_pts_ref(nqz) * this->get_dz(k + ks, n);
        quad_wts_phys(nqx, nqy, nqz) =
            this->x_quad_wts_ref(nqx) * this->dx * this->y_quad_wts_ref(nqy) *
            this->dy * this->z_quad_wts_ref(nqz) * this->get_dz(k + ks, n);
      }
    }
  }
}

template <class T>
template <class F>
void Geometry<T>::set_00form_values(F initial_value_function, Field &field,
                                    int ndof, int nz) const {
  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  parallel_for(
      "Set 00 form values",
      SimpleBounds<4>(nz, this->topology.n_cells_y, this->topology.n_cells_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
        SArray<CoordsXYZ, 1, 1> quad_pts_phys;
        SArray<real, 1, 1> quad_wts_phys;
        get_00form_quad_pts_wts(k, j, i, n, quad_pts_phys, quad_wts_phys);
        field.data(ndof, k + ks, j + js, i + is, n) =
            initial_value_function(quad_pts_phys(0).x, quad_pts_phys(0).y,
                                   quad_pts_phys(0).z) *
            quad_wts_phys(0);
      });
}

template <class T>
template <class F>
void Geometry<T>::set_00form_values(F initial_value_function, Field &field,
                                    int ndof) const {
  set_00form_values(initial_value_function, field, ndof, this->topology.ni);
}

template <class T>
template <class F>
void Geometry<T>::set_10form_values(F initial_value_function, Field &field,
                                    int ndof, int nz) const {
  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  int yedge_offset = std::numeric_limits<int>::min();
  int xedge_offset = std::numeric_limits<int>::min();

  if (ndims > 1) {
    // twisted edges
    yedge_offset = 0;
    xedge_offset = field.ndofs;

    // compared to twisted edges in 2D, straight edges are stored x/y (V, -U)
    // instead of y/x
    if (straight) {
      yedge_offset = field.ndofs;
      xedge_offset = 0;
    }
  } else {
    xedge_offset = 0;
  }

  parallel_for(
      "Set 10 form values",
      SimpleBounds<4>(nz, this->topology.n_cells_y, this->topology.n_cells_x,
                      this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
        SArray<CoordsXYZ, 1, ic_quad_pts_x> x_edge_quad_pts_phys;
        SArray<real, 1, ic_quad_pts_x> x_edge_quad_wts_phys;
        SArray<VecXYZ, 1, ic_quad_pts_x> x_edge_line_vec;

        SArray<CoordsXYZ, 1, ic_quad_pts_y> y_edge_quad_pts_phys;
        SArray<real, 1, ic_quad_pts_y> y_edge_quad_wts_phys;
        SArray<VecXYZ, 1, ic_quad_pts_y> y_edge_line_vec;

        get_10form_quad_pts_wts(k, j, i, n, x_edge_quad_pts_phys,
                                x_edge_quad_wts_phys, y_edge_quad_pts_phys,
                                y_edge_quad_wts_phys);
        get_10edge_tangents(k, j, i, n, x_edge_line_vec, y_edge_line_vec);

        // x edge
        real x_tempval = 0.0_fp;
        for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
          auto initval = initial_value_function(x_edge_quad_pts_phys(nqx).x,
                                                x_edge_quad_pts_phys(nqx).y,
                                                x_edge_quad_pts_phys(nqx).z);

          x_tempval += (initval.u * x_edge_line_vec(nqx).u +
                        initval.v * x_edge_line_vec(nqx).v +
                        initval.w * x_edge_line_vec(nqx).w) *
                       x_edge_quad_wts_phys(nqx);
        }
        field.data(ndof + xedge_offset, k + ks, j + js, i + is, n) = x_tempval;

        // y edge
        if (ndims > 1) {
          real y_tempval = 0.0_fp;
          for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
            auto initval = initial_value_function(y_edge_quad_pts_phys(nqy).x,
                                                  y_edge_quad_pts_phys(nqy).y,
                                                  y_edge_quad_pts_phys(nqy).z);

            y_tempval += (initval.u * y_edge_line_vec(nqy).u +
                          initval.v * y_edge_line_vec(nqy).v +
                          initval.w * y_edge_line_vec(nqy).w) *
                         y_edge_quad_wts_phys(nqy);
          }
          field.data(ndof + yedge_offset, k + ks, j + js, i + is, n) =
              y_tempval;
        }
      });
}

template <class T>
template <class F>
void Geometry<T>::set_10form_values(F initial_value_function, Field &field,
                                    int ndof) const {
  set_10form_values(initial_value_function, field, ndof, this->topology.ni);
}

template <class T>
template <class F>
void Geometry<T>::set_01form_values(F initial_value_function, Field &field,
                                    int ndof) const {

  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  parallel_for(
      "Set 01 form values",
      SimpleBounds<4>(this->topology.nl, this->topology.n_cells_y,
                      this->topology.n_cells_x, this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
        SArray<CoordsXYZ, 1, ic_quad_pts_z> z_edge_quad_pts_phys;
        SArray<real, 1, ic_quad_pts_z> z_edge_quad_wts_phys;
        SArray<VecXYZ, 1, ic_quad_pts_z> z_edge_line_vec;

        get_01form_quad_pts_wts(k, j, i, n, z_edge_quad_pts_phys,
                                z_edge_quad_wts_phys);
        get_01edge_tangents(k, j, i, n, z_edge_line_vec);

        real tempval = 0.0_fp;
        for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
          auto initval = initial_value_function(z_edge_quad_pts_phys(nqz).x,
                                                z_edge_quad_pts_phys(nqz).y,
                                                z_edge_quad_pts_phys(nqz).z);
          tempval += (initval.u * z_edge_line_vec(nqz).u +
                      initval.v * z_edge_line_vec(nqz).v +
                      initval.w * z_edge_line_vec(nqz).w) *
                     z_edge_quad_wts_phys(nqz);
        }
        field.data(ndof, k + ks, j + js, i + is, n) = tempval;
      });
}

template <class T>
template <class F>
void Geometry<T>::set_nm11form_values(F initial_value_function, Field &field,
                                      int ndof) const {
  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  int yz_face_offset = std::numeric_limits<int>::min();
  int xz_face_offset = std::numeric_limits<int>::min();

  if (ndims > 1) {
    // twisted edges
    yz_face_offset = 0;
    xz_face_offset = field.ndofs;

    // compared to twisted edges in 2D, straight edges are stored x/y (V, -U)
    // instead of y/x
    if (straight) {
      yz_face_offset = field.ndofs;
      xz_face_offset = 0;
    }
  } else {
    yz_face_offset = 0;
  }

  parallel_for(
      "Set nm11 form values",
      SimpleBounds<4>(this->topology.nl, this->topology.n_cells_y,
                      this->topology.n_cells_x, this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
        SArray<CoordsXYZ, 2, ic_quad_pts_x, ic_quad_pts_z>
            xz_face_quad_pts_phys;
        SArray<real, 2, ic_quad_pts_x, ic_quad_pts_z> xz_face_quad_wts_phys;
        SArray<VecXYZ, 2, ic_quad_pts_x, ic_quad_pts_z> xz_face_normal_vec;

        SArray<CoordsXYZ, 2, ic_quad_pts_y, ic_quad_pts_z>
            yz_face_quad_pts_phys;
        SArray<real, 2, ic_quad_pts_y, ic_quad_pts_z> yz_face_quad_wts_phys;
        SArray<VecXYZ, 2, ic_quad_pts_y, ic_quad_pts_z> yz_face_normal_vec;

        get_11form_quad_pts_wts(k, j, i, n, xz_face_quad_pts_phys,
                                xz_face_quad_wts_phys, yz_face_quad_pts_phys,
                                yz_face_quad_wts_phys);
        get_11face_normals(k, j, i, n, xz_face_normal_vec, yz_face_normal_vec);

        if (ndims > 1) {
          // xz face
          real xz_tempval = 0.0_fp;
          for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
            for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
              auto initval =
                  initial_value_function(xz_face_quad_pts_phys(nqx, nqz).x,
                                         xz_face_quad_pts_phys(nqx, nqz).y,
                                         xz_face_quad_pts_phys(nqx, nqz).z);

              xz_tempval += (initval.u * xz_face_normal_vec(nqx, nqz).u +
                             initval.v * xz_face_normal_vec(nqx, nqz).v +
                             initval.w * xz_face_normal_vec(nqx, nqz).w) *
                            xz_face_quad_wts_phys(nqx, nqz);
            }
          }
          field.data(ndof + xz_face_offset, k + ks, j + js, i + is, n) =
              xz_tempval;
        }

        // yz face
        real yz_tempval = 0.0_fp;
        for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
          for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
            auto initval =
                initial_value_function(yz_face_quad_pts_phys(nqy, nqz).x,
                                       yz_face_quad_pts_phys(nqy, nqz).y,
                                       yz_face_quad_pts_phys(nqy, nqz).z);

            yz_tempval += (initval.u * yz_face_normal_vec(nqy, nqz).u +
                           initval.v * yz_face_normal_vec(nqy, nqz).v +
                           initval.w * yz_face_normal_vec(nqy, nqz).w) *
                          yz_face_quad_wts_phys(nqy, nqz);
          }
        }
        field.data(ndof + yz_face_offset, k + ks, j + js, i + is, n) =
            yz_tempval;
      });
}

template <class T>
template <class F>
void Geometry<T>::set_n0form_values(F initial_value_function, Field &field,
                                    int ndof) const {
  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  parallel_for(
      "Set n0 form values",
      SimpleBounds<4>(this->topology.ni, this->topology.n_cells_y,
                      this->topology.n_cells_x, this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
        SArray<CoordsXYZ, 2, ic_quad_pts_x, ic_quad_pts_y>
            xy_face_quad_pts_phys;
        SArray<real, 2, ic_quad_pts_x, ic_quad_pts_y> xy_face_quad_wts_phys;
        SArray<VecXYZ, 2, ic_quad_pts_x, ic_quad_pts_y> xy_face_normal_vec;

        get_20form_quad_pts_wts(k, j, i, n, xy_face_quad_pts_phys,
                                xy_face_quad_wts_phys);
        get_20face_normals(k, j, i, n, xy_face_normal_vec);

        real xy_tempval = 0.0_fp;
        for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
          for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
            auto initval =
                initial_value_function(xy_face_quad_pts_phys(nqx, nqy).x,
                                       xy_face_quad_pts_phys(nqx, nqy).y,
                                       xy_face_quad_pts_phys(nqx, nqy).z);

            xy_tempval += (initval.u * xy_face_normal_vec(nqx, nqy).u +
                           initval.v * xy_face_normal_vec(nqx, nqy).v +
                           initval.w * xy_face_normal_vec(nqx, nqy).w) *
                          xy_face_quad_wts_phys(nqx, nqy);
          }
        }

        field.data(ndof, k + ks, j + js, i + is, n) = xy_tempval;
      });
}

template <class T>
template <class F>
void Geometry<T>::set_n1form_values(F initial_value_function, Field &field,
                                    int ndof) const {

  int is = this->topology.is;
  int js = this->topology.js;
  int ks = this->topology.ks;

  parallel_for(
      "Set n1 form values",
      SimpleBounds<4>(this->topology.nl, this->topology.n_cells_y,
                      this->topology.n_cells_x, this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
        SArray<CoordsXYZ, 3, ic_quad_pts_x, ic_quad_pts_y, ic_quad_pts_z>
            quad_pts_phys;
        SArray<real, 3, ic_quad_pts_x, ic_quad_pts_y, ic_quad_pts_z>
            quad_wts_phys;
        get_21form_quad_pts_wts(k, j, i, n, quad_pts_phys, quad_wts_phys);

        real tempval = 0.0_fp;
        for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
          for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
            for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
              tempval +=
                  initial_value_function(quad_pts_phys(nqx, nqy, nqz).x,
                                         quad_pts_phys(nqx, nqy, nqz).y,
                                         quad_pts_phys(nqx, nqy, nqz).z) *
                  quad_wts_phys(nqx, nqy, nqz);
            }
          }
        }

        field.data(ndof, k + ks, j + js, i + is, n) = tempval;
      });
}

template <class T>
template <class F>
void Geometry<T>::set_profile_00form_values(F initial_value_function,
                                            Profile &prof, int ndof) const {
  int ks = this->topology.ks;

  parallel_for(
      "Set profile 00 form values",
      SimpleBounds<2>(this->topology.ni, this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int n) {
        int i = 0;
        int j = 0;
        SArray<CoordsXYZ, 1, 1> quad_pts_phys;
        SArray<real, 1, 1> quad_wts_phys;
        get_00form_quad_pts_wts(k, j, i, n, quad_pts_phys, quad_wts_phys);
        prof.data(ndof, k + ks, n) =
            initial_value_function(quad_pts_phys(0).z) * quad_wts_phys(0);
      });
}

template <class T>
template <class F>
void Geometry<T>::set_profile_n1form_values(F initial_value_function,
                                            Profile &prof, int ndof) const {

  int ks = this->topology.ks;

  parallel_for(
      "Set profile n1 form values",
      SimpleBounds<2>(this->topology.nl, this->topology.nens),
      YAKL_CLASS_LAMBDA(int k, int n) {
        SArray<CoordsXYZ, 3, ic_quad_pts_x, ic_quad_pts_y, ic_quad_pts_z>
            quad_pts_phys;
        SArray<real, 3, ic_quad_pts_x, ic_quad_pts_y, ic_quad_pts_z>
            quad_wts_phys;

        int i = 0;
        int j = 0;
        get_21form_quad_pts_wts(k, j, i, n, quad_pts_phys, quad_wts_phys);

        real tempval = 0.0_fp;
        for (int nqx = 0; nqx < ic_quad_pts_x; nqx++) {
          for (int nqy = 0; nqy < ic_quad_pts_y; nqy++) {
            for (int nqz = 0; nqz < ic_quad_pts_z; nqz++) {
              tempval +=
                  initial_value_function(quad_pts_phys(nqx, nqy, nqz).z) *
                  quad_wts_phys(nqx, nqy, nqz);
            }
          }
        }

        prof.data(ndof, k + ks, n) = tempval;
      });
}

template <class T>
template <class F>
void Geometry<T>::set_0form_values(F initial_value_function, Field &field,
                                   int ndof) const {
  this->set_00form_values(
      YAKL_LAMBDA(real x, real y, real z) {
        return initial_value_function(x, y);
      },
      field, ndof, this->topology.nl);
}
template <class T>
template <class F>
void Geometry<T>::set_2form_values(F initial_value_function, Field &field,
                                   int ndof) const {
  this->set_n1form_values(
      YAKL_LAMBDA(real x, real y, real z) {
        return initial_value_function(x, y);
      },
      field, ndof);
}

template <class T>
template <class F>
void Geometry<T>::set_1form_values(F initial_value_function, Field &field,
                                   int ndof,
                                   LINE_INTEGRAL_TYPE line_type) const {

  auto f_3d = YAKL_LAMBDA(real x, real y, real z) {
    auto v_2d = initial_value_function(x, y);
    VecXYZ v_3d;
    v_3d.u = v_2d.u;
    v_3d.v = v_2d.v;
    v_3d.w = 0;
    return v_3d;
  };

  if (line_type == LINE_INTEGRAL_TYPE::NORMAL) {
    this->set_nm11form_values(f_3d, field, ndof);
  }
  if (line_type == LINE_INTEGRAL_TYPE::TANGENT) {
    this->set_10form_values(f_3d, field, ndof, this->topology.nl);
  }
}
