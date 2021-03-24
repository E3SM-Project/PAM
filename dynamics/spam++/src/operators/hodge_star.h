
#ifndef _HODGE_STARS_H_
#define _HODGE_STARS_H_

#include "common.h"
#include "geometry.h"


void YAKL_INLINE H(SArray<real,ndims> &var, SArray<real,ndims,1> const &velocity, SArray<real,ndims,1> const &Hgeom) {

  for (int d=0; d<ndims; d++) {
        var(d) = Hgeom(d,0) * velocity(d,0);
}
}

void YAKL_INLINE H(SArray<real,ndims> &var, SArray<real,ndims,3> const &velocity, SArray<real,ndims,3> const &Hgeom) {

  for (int d=0; d<ndims; d++) {
        var(d) = -1./24.* Hgeom(d,0) * velocity(d,0) + 26./24.* Hgeom(d,1) * velocity(d,1) - 1./24.* Hgeom(d,2) * velocity(d,2);
}
}

void YAKL_INLINE H(SArray<real,ndims> &var, SArray<real,ndims,5> const &velocity, SArray<real,ndims,5> const &Hgeom) {

  for (int d=0; d<ndims; d++) {
  var(d) = 9./1920.*Hgeom(d,0) * velocity(d,0) - 116./1920.*Hgeom(d,1) * velocity(d,1) + 2134./1920.* Hgeom(d,2) * velocity(d,2) - 116./1920.* Hgeom(d,3) * velocity(d,3) + 9./1920.* Hgeom(d,4) * velocity(d,4);
}
}


// IS THE INDEXING HERE CORRECT FOR NEW TOPOLOGY/GEOMETRY FORMATS?
template<uint ndofs, uint ord, uint off=ord/2 -1> void YAKL_INLINE compute_H(SArray<real,ndims> &u, const realArr vvar, Geometry<ndims,1,1,1> &pgeom, Geometry<ndims,1,1,1> &dgeom, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndims,ord-1> v;
  SArray<real,ndims,ord-1> Hgeom;
  for (int p=0; p<ord-1; p++) {
  for (int d=0; d<ndims; d++) {
    if (d==0) {
    v(d,p) = vvar(d, k+ks, j+js, i+is+p-off);
    Hgeom(d,p) = pgeom.get_area_1form(d, k+ks, j+js, i+is+p-off) / dgeom.get_area_1form(d, k+ks, j+js, i+is+p-off);
  }
  if (d==1) {
  v(d,p) = vvar(d, k+ks, j+js+p-off, i+is);
  Hgeom(d,p) = pgeom.get_area_1form(d, k+ks, j+js+p-off, i+is) / dgeom.get_area_1form(d, k+ks, j+js, i+is+p-off);
  }
  if (d==2) {
  v(d,p) = vvar(d, k+ks+p-off, j+js, i+is);
  Hgeom(d,p) = pgeom.get_area_1form(d, k+ks+p-off, j+js, i+is) / dgeom.get_area_1form(d, k+ks+p-off, j+js, i+is);
  }
  }}
  H(u, v, Hgeom);

}

template<uint ndofs, uint ord, ADD_MODE addmode=ADD_MODE::REPLACE, uint off=ord/2 -1> void YAKL_INLINE compute_H(realArr uvar, const realArr vvar, Geometry<ndims,1,1,1> &pgeom, Geometry<ndims,1,1,1> &dgeom, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndims> u;
  compute_H<ndofs, ord, off> (u, vvar, pgeom, dgeom, is, js, ks, i, j, k);
  if (addmode == ADD_MODE::REPLACE) {for (int d=0; d<ndims; d++) { uvar(d, k+ks, j+js, i+is) = u(d);}}
  if (addmode == ADD_MODE::ADD) {for (int d=0; d<ndims; d++) { uvar(d, k+ks, j+js, i+is) += u(d);}}
}



template<uint ndofs> void YAKL_INLINE I(SArray<real,ndofs> &var, SArray<real,ndofs,ndims,1> const &dens, SArray<real,ndims,1> const &Igeom) {

  for (int l=0; l<ndofs; l++) {
        var(l) = Igeom(0,0) * dens(l,0,0);
        }
}

template<uint ndofs> void YAKL_INLINE I(SArray<real,ndofs> &var, SArray<real,ndofs,ndims,3> const &dens, SArray<real,ndims,3> const &Igeom) {
  for (int l=0; l<ndofs; l++) {
    var(l) = Igeom(0,1) * dens(l,0,1);
    for (int d=0; d<ndims; d++) {
        var(l) += -1./48.* Igeom(d,0) * dens(l,d,0) + 2./48.* Igeom(d,1) * dens(l,d,1) - 1./48.* Igeom(d,2) * dens(l,d,2);
      }
    }
}

template<uint ndofs> void YAKL_INLINE I(SArray<real,ndofs> &var, SArray<real,ndofs,ndims,5> const &dens, SArray<real,ndims,5> const &Igeom) {
  for (int l=0; l<ndofs; l++) {
    var(l) = Igeom(0,2) * dens(l,0,2);
    for (int d=0; d<ndims; d++) {
          var(l) += 1./576.*Igeom(d,0) * dens(l,d,0) - 16./576.*Igeom(d,1) * dens(l,d,1) + 30./576.* Igeom(d,2) * dens(l,d,2) - 16./576.* Igeom(d,3) * dens(l,d,3) + 1./576.* Igeom(d,4) * dens(l,d,4);
}
}
}



// IS THE INDEXING HERE CORRECT FOR NEW TOPOLOGY/GEOMETRY FORMATS?

template<uint ndofs, uint ord, uint off=ord/2 -1> void YAKL_INLINE compute_I(SArray<real,ndofs> &x0, const realArr var, Geometry<ndims,1,1,1> &pgeom, Geometry<ndims,1,1,1> &dgeom, int is, int js, int ks, int i, int j, int k)
{
SArray<real,ndofs,ndims,ord-1> x;
SArray<real,ndims,ord-1> Igeom;
for (int p=0; p<ord-1; p++) {
for (int l=0; l<ndofs; l++) {
for (int d=0; d<ndims; d++) {
  if (d==0)
  {
x(l,d,p) = var(l, k+ks, j+js, i+is+p-off);
Igeom(d,p) = pgeom.get_area_0form(k+ks, j+js, i+is+p-off) / dgeom.get_area_2form(k+ks, j+js, i+is+p-off);
}
if (d==1)
{
x(l,d,p) = var(l, k+ks, j+js+p-off, i+is);
Igeom(d,p) = pgeom.get_area_0form(k+ks, j+js+p-off, i+is) / dgeom.get_area_2form(k+ks, j+js+p-off, i+is);
}
if (d==2)
{
x(l,d,p) = var(l, k+ks+p-off, j+js, i+is);
Igeom(d,p) = pgeom.get_area_0form(k+ks+p-off, j+js, i+is) / dgeom.get_area_2form(k+ks+p-off, j+js, i+is);
}
}}}
I<ndofs>(x0, x, Igeom);
}

template<uint ndofs, uint ord, ADD_MODE addmode=ADD_MODE::REPLACE, uint off=ord/2 -1> void YAKL_INLINE compute_I(realArr var0, const realArr var, Geometry<ndims,1,1,1> &pgeom, Geometry<ndims,1,1,1> &dgeom, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndofs> x0;
  compute_I<ndofs, ord, off> (x0, var, pgeom, dgeom, is, js, ks, i, j, k);
  if (addmode == ADD_MODE::REPLACE) {for (int l=0; l<ndofs; l++) {var0(l, k+ks, j+js, i+is) = x0(l);}}
  if (addmode == ADD_MODE::ADD) {for (int l=0; l<ndofs; l++) {var0(l, k+ks, j+js, i+is) += x0(l);}}

}


template<uint ndofs> void YAKL_INLINE J(SArray<real,ndofs> &var, SArray<real,ndofs,ndims,1> const &dens, SArray<real,ndims,1> const &Jgeom) {

  for (int l=0; l<ndofs; l++) {
        var(l) = Jgeom(0,0) * dens(l,0,0);
        }
}

template<uint ndofs> void YAKL_INLINE J(SArray<real,ndofs> &var, SArray<real,ndofs,ndims,3> const &dens, SArray<real,ndims,3> const &Jgeom) {
  for (int l=0; l<ndofs; l++) {
    var(l) = Jgeom(0,1) * dens(l,0,1);
    for (int d=0; d<ndims; d++) {
        var(l) += -1./48.* Jgeom(d,0) * dens(l,d,0) + 2./48.* Jgeom(d,1) * dens(l,d,1) - 1./48.* Jgeom(d,2) * dens(l,d,2);
      }
    }
}

template<uint ndofs> void YAKL_INLINE J(SArray<real,ndofs> &var, SArray<real,ndofs,ndims,5> const &dens, SArray<real,ndims,5> const &Jgeom) {
  for (int l=0; l<ndofs; l++) {
    var(l) = Jgeom(0,2) * dens(l,0,2);
    for (int d=0; d<ndims; d++) {
          var(l) += 1./576.*Jgeom(d,0) * dens(l,d,0) - 16./576.*Jgeom(d,1) * dens(l,d,1) + 30./576.* Jgeom(d,2) * dens(l,d,2) - 16./576.* Jgeom(d,3) * dens(l,d,3) + 1./576.* Jgeom(d,4) * dens(l,d,4);
}
}
}

// IS THE INDEXING HERE CORRECT FOR NEW TOPOLOGY/GEOMETRY FORMATS?

template<uint ndofs, uint ord, uint off=ord/2 -1> void YAKL_INLINE compute_J(SArray<real,ndofs> &x0, const realArr var, Geometry<ndims,1,1,1> &pgeom, Geometry<ndims,1,1,1> &dgeom, int is, int js, int ks, int i, int j, int k)
{
SArray<real,ndofs,ndims,ord-1> x;
SArray<real,ndims,ord-1> Jgeom;
for (int p=0; p<ord-1; p++) {
for (int l=0; l<ndofs; l++) {
for (int d=0; d<ndims; d++) {
  if (d==0)
  {
x(l,d,p) = var(l, k+ks, j+js, i+is+p-off);
Jgeom(d,p) = dgeom.get_area_0form(k+ks, j+js, i+is+p-off) / pgeom.get_area_2form(k+ks, j+js, i+is+p-off);
}
if (d==1)
{
x(l,d,p) = var(l, k+ks, j+js+p-off, i+is);
Jgeom(d,p) = dgeom.get_area_0form(k+ks, j+js+p-off, i+is) / pgeom.get_area_2form(k+ks, j+js+p-off, i+is);
}
if (d==2)
{
x(l,d,p) = var(l, k+ks+p-off, j+js, i+is);
Jgeom(d,p) = dgeom.get_area_0form(k+ks+p-off, j+js, i+is) / pgeom.get_area_2form(k+ks+p-off, j+js, i+is);
}
}}}
J<ndofs>(x0, x, Jgeom);
}

template<uint ndofs, uint ord, ADD_MODE addmode=ADD_MODE::REPLACE, uint off=ord/2 -1> void YAKL_INLINE compute_J(realArr var0, const realArr var, Geometry<ndims,1,1,1> &pgeom, Geometry<ndims,1,1,1> &dgeom, int is, int js, int ks, int i, int j, int k)
{
  SArray<real,ndofs> x0;
  compute_J<ndofs,ord,off> (x0, var, pgeom, dgeom, is, js, ks, i, j, k);
  if (addmode == ADD_MODE::REPLACE) {for (int l=0; l<ndofs; l++) {var0(l, k+ks, j+js, i+is) = x0(l);}}
  if (addmode == ADD_MODE::ADD) {for (int l=0; l<ndofs; l++) {var0(l, k+ks, j+js, i+is) += x0(l);}}
}
#endif
