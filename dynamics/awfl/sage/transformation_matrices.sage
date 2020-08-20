#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python

load("poly_utils.sage")
load("quadrature_utils.sage")


#Matrices that convert N GLL points on the domain [-dx/2,dx/2] into Nth-order-accurate polynomial coefficients and vice versa
def points_gll_to_coefs(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    gll = points_gll(N)
    pnts = vector([ p.subs(x=gll[i]) for i in range(N) ])
    coefs_to_pnts = jacobian(pnts,coefs)
    pnts_to_coefs = coefs_to_pnts^-1
    return pnts_to_coefs,coefs_to_pnts

#Matrices that convert N inner-element equally spaced cells on the domain [-dx/2,dx/2] into Nth-order-accurate polynomial coefficients and vice versa
def points_cellsequal_to_coefs(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    R = RealField( 129 )
    pnts = vector([ R(i)/R(N) - 1/2 for i in range(N+1) ])
    cells = vector([ integrate_poly(N,p,x,pnts[i],pnts[i+1]) / (pnts[i+1] - pnts[i]) for i in range(N) ])
    coefs_to_cells = jacobian(cells,coefs)
    cells_to_coefs = coefs_to_cells^-1
    return cells_to_coefs,coefs_to_cells

#Matrices that convert N inner-element equally spaced cells on the domain [-dx/2,dx/2] into Nth-order-accurate polynomial coefficients and vice versa
def points_cellsequal_to_coefs_dx(N,dx) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    R = RealField( 129 )
    pnts = vector([ R(i)/R(N) - 1/2 for i in range(N+1) ])
    cells = vector([ integrate_poly(N,p,x,pnts[i]*dx,pnts[i+1]*dx) / (pnts[i+1]*dx - pnts[i]*dx) for i in range(N) ])
    coefs_to_cells = jacobian(cells,coefs)
    cells_to_coefs = coefs_to_cells^-1
    return cells_to_coefs,coefs_to_cells

#Matrices that convert N stencil averages centered about zero with dx grid spacing into Nth-order-accurate polynomial coefficients and vice versa
def stencil_to_coefs(N) :
    var('x')
    hs = (N-1)/2
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    constr = vector([ integrate(p,x,(2*i-1)/2,(2*i+1)/2) for i in range(-hs,hs+1) ])
    coefs_to_constr = jacobian(constr,coefs)
    constr_to_coefs = coefs_to_constr^-1
    return constr_to_coefs,coefs_to_constr

# This will be an N x N matrix
# Matrix to transform DOFs between coefficients and CSFV DOFs inside a single element
# Cell averages on gll spaced control volumes
def csfv_to_coefs(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    # N-1 point will bound the N-2 cell averages inside the element
    gll = points_gll(N-1)
    # size N+1 because we'll have N-1 point values (fluxes) and two cell-edge derivatives
    constr = vector([ 0*x for i in range(N) ])
    constr[0  ] = p.subs(x=-1/2)                              #Left  point value
    constr[1:N-1] = vector([ integrate_poly(N,p,x,gll[i],gll[i+1]) / (gll[i+1]-gll[i]) for i in range(N-2) ])  #Cell averages
    constr[N-1] = p.subs(x= 1/2)                              #Right point value
    #Convert to a martix
    coefs_to_constr = jacobian(constr,coefs)
    constr_to_coefs = coefs_to_constr^-1
    return constr_to_coefs,coefs_to_constr

# This will be an (N+1) x N matrix
# Matrix to transform polynomial coefficients into DOFs used to update Constrained Spectral Finite-Volume
def coefs_to_csfv_update(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    # N-1 point will bound the N-2 cell averages inside the element
    gll = points_gll(N-1)
    # size N+1 because we'll have N-1 point values (fluxes) and two cell-edge derivatives
    constr = vector([ 0*x for i in range(N+1) ])
    constr[0] = p.diff(x).subs(x=-1/2)                              #Left  Derivative
    constr[1:N] = vector([ p.subs(x=gll[i]) for i in range(N-1) ])  #Flux  Values
    constr[N] = p.diff(x).subs(x= 1/2)                              #Right Derivative
    #Convert to a martix
    coefs_to_constr = jacobian(constr,coefs)
    return coefs_to_constr

# This will be an N x (N+1) matrix
# N always has to be an odd order of accuracy
# Matrix to transform polynomial coefficients into DOFs used to update Multi-moment Constrained finite-Volume (MCV)
def coefs_to_mcv_update(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    # N-Edge (Number of DOFs per edge)
    NE = (N + 1) / 2
    # size N+1 because we'll have N-1 point values (fluxes) and two cell-edge derivatives
    constr = vector([ 0*x for i in range(N+1) ])
    constr[0:NE] = vector([ p.diff(x,i).subs(x=-1/2) for i in range(NE) ])
    constr[NE:N+1] = vector([ p.diff(x,i).subs(x= 1/2) for i in range(NE) ])
    #Convert to a martix
    coefs_to_constr = jacobian(constr,coefs)
    return coefs_to_constr

# This will be an N x (N+1)/2 matrix
# N always has to be an odd order of accuracy
# Matrix to transform polynomial coefficients into DOFs used to update Weno MCV right flux
def coefs_to_mcv_update_R(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    # N-Edge (Number of DOFs on the edge)
    NE = (N + 1) / 2
    constr = vector([ 0*x for i in range((N+1)/2) ])

    R = RealField( 129 )
    pnts = vector([ R(i)/R(N) - 1/2 for i in range(N+1) ])
    constr[0:NE] = vector([ p.diff(x,i).subs(x=pnts[(N-1)/2]) for i in range(NE) ])
    coefs_to_constr = jacobian(constr,coefs)
    return coefs_to_constr

# This will be an N x (N+1)/2 matrix
# N always has to be an odd order of accuracy
# Matrix to transform polynomial coefficients into DOFs used to update Weno MCV left flux
def coefs_to_mcv_update_L(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    # N-Edge (Number of DOFs on the edge)
    NE = (N + 1) / 2
    constr = vector([ 0*x for i in range((N+1)/2) ])

    R = RealField( 129 )
    pnts = vector([ R(i)/R(N) - 1/2 for i in range(N+1) ])
    constr[0:NE] = vector([ p.diff(x,i).subs(x=pnts[(N+1)/2]) for i in range(NE) ])
    coefs_to_constr = jacobian(constr,coefs)
    return coefs_to_constr


# This will be an N x N matrix
# N always has to be an odd order of accuracy
# Matrix to transform DOFs between coefficients and Multi-moment Constrained finite-Volume (MCV) DOFs inside a single element
def mcv_to_coefs(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    # N-Edge (Number of DOFs per edge)
    NE = (N - 1) / 2
    # size N+1 because we'll have N-1 point values (fluxes) and two cell-edge derivatives
    constr = vector([ 0*x for i in range(N) ])
    constr[0:NE]   = vector([ p.diff(x,i).subs(x=-1/2) for i in range(NE) ])
    constr[NE]     = integrate_poly(N,p,x,-1/2,1/2)   #Cell average
    constr[NE+1:N] = vector([ p.diff(x,i).subs(x= 1/2) for i in range(NE) ])
    #Convert to a martix
    coefs_to_constr = jacobian(constr,coefs)
    constr_to_coefs = coefs_to_constr^-1
    return constr_to_coefs,coefs_to_constr


# This will be an N x N matrix
# N always has to be an odd order of accuracy
# Matrix to transform DOFs between coefficients and Multi-moment Constrained finite-Volume (MCV) DOFs inside a single element
def mcv_to_coefs_dx(N,dx) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    # N-Edge (Number of DOFs per edge)
    NE = (N - 1) / 2
    # size N+1 because we'll have N-1 point values (fluxes) and two cell-edge derivatives
    constr = vector([ 0*x for i in range(N) ])
    constr[0:NE]   = vector([ p.diff(x,i).subs(x=-dx/2) for i in range(NE) ])
    constr[NE]     = integrate_poly(N,p,x,-dx/2,dx/2)/dx   #Cell average
    constr[NE+1:N] = vector([ p.diff(x,i).subs(x= dx/2) for i in range(NE) ])
    #Convert to a martix
    coefs_to_constr = jacobian(constr,coefs)
    constr_to_coefs = coefs_to_constr^-1
    return constr_to_coefs,coefs_to_constr

#Matrix that converts polynomial coefficients into differentiated polynomial coefficients
def coefs_to_deriv(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    dp = diff(p,x)
    deriv_coefs = compute_coefs(N,dp,x)
    coefs_to_deriv = jacobian(deriv_coefs,coefs)
    return coefs_to_deriv


#Compute an Nth-order-accurage polynomial from N stencil averages. Then project that to fewer GLL point values
def sten_to_gll_lower(N) :
    var('x')
    #Compute a Nth-order-accurate polynomial from a stencil of N values
    gs = (N-1)/2
    coefs = coefs_1d(N,0,'a')
    uvals = coefs_1d(N,0,'u')
    p = poly_1d(N,coefs,x)
    constr = vector([ integrate(p,x,(2*i-1)/2,(2*i+1)/2) for i in range(-gs,gs+1) ])
    p = poly_1d( N , jacobian(constr,coefs)^-1 * uvals , x )
    sten_to_gll  = [ [ [ 0*x for i in range(N) ] for j in range(N) ] for k in range(N) ]
    for j in range(1,N+1) :
        #Compute j GLL points from the polynomial
        if (j == 1) :
            locs = vector([ 0*x ])
        else :
            locs = points_gll(j)
        gll = vector([ p.subs(x=locs[i]).expand() for i in range(j) ])
        s2g = jacobian(gll,uvals)
        for i2 in range(j) :
            for i1 in range(N) :
                sten_to_gll[j-1][i1][i2] = s2g[i2][i1].n(129)
    return sten_to_gll


#Project Nth-order-accurate polynomial coefficients to fewer GLL point values
def coefs_to_gll_lower(N) :
    var('x')
    #Compute a Nth-order-accurate polynomial from a stencil of N values
    gs = (N-1)/2
    coefs = coefs_1d(N,0,'a')
    p = poly_1d( N , coefs , x )
    coefs_to_gll  = [ [ [ 0*x for i in range(N) ] for j in range(N) ] for k in range(N) ]
    for j in range(1,N+1) :
        #Compute j GLL points from the polynomial
        if (j == 1) :
            locs = vector([ 0*x ])
        else :
            locs = points_gll(j)
        gll = vector([ p.subs(x=locs[i]).expand() for i in range(j) ])
        c2g = jacobian(gll,coefs)
        for i2 in range(j) :
            for i1 in range(N) :
                coefs_to_gll[j-1][i1][i2] = c2g[i2][i1]
    return coefs_to_gll


def weno_sten_to_coefs(N) :
    var('x')
    hs = (N-1)/2
    wenopolys = [[[ 0*x for i in range(N) ] for j in range(N) ] for k in range(hs+2) ]
    locs = vector([ (2*i-1)/2 for i in range(-hs,hs+2) ])
    #Lower-ordered polynomials
    for j in range(hs+1) :
        coefs = coefs_1d(hs+1,0,'a')
        p = poly_1d(hs+1,coefs,x)
        constr = vector([ integrate(p,x,locs[i],locs[i+1]) / (locs[i+1]-locs[i]) for i in range(j,j+hs+1) ])
        tmpl = force_fp( jacobian(constr,coefs)^-1 , 129 )
        for i2 in range(hs+1) :
            for i1 in range(hs+1) :
                wenopolys[j][i1][i2] = tmpl[i2][i1]
    #Higher-order polynomial
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    constr = vector([ integrate(p,x,locs[i],locs[i+1]) / (locs[i+1]-locs[i]) for i in range(N) ])
    tmph = force_fp( jacobian(constr,coefs)^-1 , 129 )
    for i2 in range(N) :
        for i1 in range(N) :
            wenopolys[hs+1][i1][i2] = tmph[i2][i1]
    return wenopolys


def coefs_to_TV(N) :
    var('x')
    hs = (N-1)/2
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    TV = sum( vector([ integrate_poly(N, p.diff(x,i)^2 ,x,-1/2,1/2) for i in range(1,N)]) )
    return TV


#Matrix that converts polynomial coefficients into differentiated polynomial coefficients
def coefs_to_prim(N) :
    var('x')
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    pp = integrate(p,x)
    prim_coefs = compute_coefs(N,pp,x)
    coefs_to_prim = jacobian(prim_coefs,coefs)
    return coefs_to_prim
