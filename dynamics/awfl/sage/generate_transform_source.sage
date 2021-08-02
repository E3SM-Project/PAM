#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python

load("transformation_matrices.sage")
load("c_utils.sage")

N1 = 2
N2 = 9

print('#pragma once')

print('#include "awfl_const.h"')
print('#include <math.h>')
print('using yakl::SArray;\n')

print('namespace TransformMatrices {\n')

for N in range(N1,N2+1) :
    hs = (N-1)/2

    print('  template <class FP> YAKL_INLINE void get_gll_points(SArray<FP,1,%s> &rslt) {'%(N))
    pts,wts = lobatto_weights_nodes(N,129,False,1.e-35,100)
    pts = pts / 2
    wts = wts / 2
    print(add_spaces(4,c_vector('rslt',N,force_fp(pts,129),'none',200)))
    print('  }\n');

    print('  template <class FP> YAKL_INLINE void get_gll_weights(SArray<FP,1,%s> &rslt) {'%(N))
    pts,wts = lobatto_weights_nodes(N,129,False,1.e-35,100)
    pts = pts / 2
    wts = wts / 2
    print(add_spaces(4,c_vector('rslt',N,force_fp(wts,129),'none',200)))
    print('  }\n');

    print('  template <class FP> YAKL_INLINE void cells_equal_to_coefs(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
    p2c,c2p = points_cellsequal_to_coefs(N)
    print(add_spaces(4,c_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
    print('  }\n');
    
    print('  template <class FP> YAKL_INLINE void coefs_to_cells_equal(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
    p2c,c2p = points_cellsequal_to_coefs(N)
    print(add_spaces(4,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
    print('  }\n');

    print('  template <class FP> YAKL_INLINE void cells_equal_to_coefs(FP dx, SArray<FP,2,%s,%s> &rslt) {'%(N,N))
    var('dx')
    p2c,c2p = points_cellsequal_to_coefs_dx(N,dx)
    print(add_spaces(4,c_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
    print('  }\n');
    
    print('  template <class FP> YAKL_INLINE void coefs_to_cells_equal(FP dx, SArray<FP,2,%s,%s> &rslt) {'%(N,N))
    var('dx')
    p2c,c2p = points_cellsequal_to_coefs_dx(N,dx)
    print(add_spaces(4,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
    print('  }\n');

    print('  template <class FP> YAKL_INLINE void gll_to_coefs(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
    p2c,c2p = points_gll_to_coefs(N)
    print(add_spaces(4,c_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
    print('  }\n');

    print('  template <class FP> YAKL_INLINE void coefs_to_gll(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
    p2c,c2p = points_gll_to_coefs(N)
    print(add_spaces(4,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
    print('  }\n');

    print('  template <class FP> YAKL_INLINE void coefs_to_deriv(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
    c2d = coefs_to_deriv(N)
    print(add_spaces(4,c_matrix('rslt',N,N,force_fp(c2d,129),'none',200)))
    print('  }\n');

    print('  template <class FP> YAKL_INLINE FP coefs_to_tv(SArray<FP,1,%s> &a) {'%(N))
    print('    FP rslt;')
    rslt = coefs_to_TV(N)
    print(add_spaces(4,c_scalar('rslt',force_fp(rslt,129),'a',200)))
    print('    return rslt;')
    print('  }\n');

    if (N > 2) :
        print('  template <class FP> YAKL_INLINE void csfv_to_coefs(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
        p2c,c2p = csfv_to_coefs(N)
        print(add_spaces(2,c_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
        print(' }\n');
        
        print('  template <class FP> YAKL_INLINE void coefs_to_csfv(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
        p2c,c2p = csfv_to_coefs(N)
        print(add_spaces(2,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
        print(' }\n');

        print('  template <class FP> YAKL_INLINE void coefs_to_csfv_update(SArray<FP,2,%s,%s> &rslt) {'%(N,N+1))
        c2p = coefs_to_csfv_update(N)
        print(add_spaces(2,c_matrix('rslt',N+1,N,force_fp(c2p,129),'none',200)))
        print(' }\n');
        

    if (N%2 == 1) :
        print('  template <class FP> YAKL_INLINE void sten_to_coefs(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
        p2c,c2p = stencil_to_coefs(N)
        print(add_spaces(4,c_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
        print('  }\n');

        print('  template <class FP> YAKL_INLINE void coefs_to_sten(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
        p2c,c2p = stencil_to_coefs(N)
        print(add_spaces(4,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
        print('  }\n');

        s2g = sten_to_gll_lower(N)
        for M in range(1,N+1) :
            print('  template <class FP> YAKL_INLINE void sten_to_gll_lower(SArray<FP,2,%s,%s> &rslt) {'%(N,M))
            print(add_spaces(4,c_matrix_aoa('rslt',M,N,s2g[M-1],'none',200)))
            print('  }\n');

        c2g = coefs_to_gll_lower(N)
        for M in range(1,N+1) :
            print('  template <class FP> YAKL_INLINE void coefs_to_gll_lower(SArray<FP,2,%s,%s> &rslt) {'%(N,M))
            print(add_spaces(4,c_matrix_aoa('rslt',M,N,c2g[M-1],'none',200)))
            print('  }\n');

        print('  template <class FP> YAKL_INLINE void weno_lower_sten_to_coefs(SArray<FP,3,%s,%s,%s> &rslt) {'%(hs+1,hs+1,hs+1))
        weno = weno_lower_sten_to_coefs(N)
        print(add_spaces(4,c_3d('rslt',hs+1,hs+1,hs+1,weno,'none',200)))
        print('  }\n');

        print('  template <class FP> YAKL_INLINE void mcv_to_coefs(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
        p2c,c2p = mcv_to_coefs(N)
        print(add_spaces(4,c_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
        print('  }\n');

        print('  template <class FP> YAKL_INLINE void coefs_to_mcv(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
        p2c,c2p = mcv_to_coefs(N)
        print(add_spaces(4,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
        print('  }\n');

        print('  template <class FP> YAKL_INLINE void mcv_to_coefs(FP dx, SArray<FP,2,%s,%s> &rslt) {'%(N,N))
        var('dx')
        p2c,c2p = mcv_to_coefs_dx(N,dx)
        print(add_spaces(4,c_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
        print('  }\n');

        print('  template <class FP> YAKL_INLINE void coefs_to_mcv(FP dx, SArray<FP,2,%s,%s> &rslt) {'%(N,N))
        var('dx')
        p2c,c2p = mcv_to_coefs_dx(N,dx)
        print(add_spaces(4,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
        print('  }\n');
        
        print('  template <class FP> YAKL_INLINE void coefs_to_mcv_update(SArray<FP,2,%s,%s> &rslt) {'%(N,N+1))
        c2p = coefs_to_mcv_update(N)
        print(add_spaces(4,c_matrix('rslt',N+1,N,force_fp(c2p,129),'none',200)))
        print('  }\n');

        print('  template <class FP> YAKL_INLINE void coefs_to_mcv_update_R(SArray<FP,2,%s,%s> &rslt) {'%(N,(N+1)/2))
        c2p = coefs_to_mcv_update_R(N)
        print(add_spaces(4,c_matrix('rslt',(N+1)/2,N,force_fp(c2p,129),'none',200)))
        print('  }\n');

        print('  template <class FP> YAKL_INLINE void coefs_to_mcv_update_L(SArray<FP,2,%s,%s> &rslt) {'%(N,(N+1)/2))
        c2p = coefs_to_mcv_update_L(N)
        print(add_spaces(4,c_matrix('rslt',(N+1)/2,N,force_fp(c2p,129),'none',200)))
        print('  }\n');



N = 1

print('  template <class FP> YAKL_INLINE void get_gll_points(SArray<FP,1,%s> &rslt) {'%(N))
print('    rslt(0) = 0;')
print('  }\n');

print('  template <class FP> YAKL_INLINE void get_gll_weights(SArray<FP,1,%s> &rslt) {'%(N))
print('    rslt(0) = 1;')
print('  }\n');

print('  template <class FP> YAKL_INLINE void gll_to_coefs(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
print('    rslt(0,0) = 1;')
print('  }\n');

print('  template <class FP> YAKL_INLINE void coefs_to_gll(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
print('    rslt(0,0) = 1;')
print('  }\n');

print('  template <class FP> YAKL_INLINE void coefs_to_deriv(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
print('    rslt(0,0) = 0;')
print('  }\n');

print('  template <class FP> YAKL_INLINE FP coefs_to_tv(SArray<FP,1,%s> &a) {'%(N))
print('    return 0;')
print('  }\n');

if (N%2 == 1) :
    print('  template <class FP> YAKL_INLINE void sten_to_coefs(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
    print('    rslt(0,0) = 1;')
    print('  }\n');

    print('  template <class FP> YAKL_INLINE void coefs_to_sten(SArray<FP,2,%s,%s> &rslt) {'%(N,N))
    print('    rslt(0,0) = 1;')
    print('  }\n');

    s2g = sten_to_gll_lower(N)
    for M in range(1,N+1) :
        print('  template <class FP> YAKL_INLINE void sten_to_gll_lower(SArray<FP,2,%s,%s> &rslt) {'%(N,M))
        print('    rslt(0,0) = 1;')
        print('  }\n');

    c2g = coefs_to_gll_lower(N)
    for M in range(1,N+1) :
        print('  template <class FP> YAKL_INLINE void coefs_to_gll_lower(SArray<FP,2,%s,%s> &rslt) {'%(N,M))
        print('    rslt(0,0) = 1;')
        print('  }\n');

    print('  template <class FP> YAKL_INLINE void weno_sten_to_coefs(SArray<FP,3,%s,%s,%s> &rslt) {'%(N,N,N))
    print('    rslt(0,0,0) = 1;')
    print('  }\n');

print('}\n')

