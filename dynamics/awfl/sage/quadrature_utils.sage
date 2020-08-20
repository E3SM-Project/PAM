#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python


############################################################################
## Computes Gauss-Lobatto quadrature weights and nodes to integrate over
## the domain [-1,1]. Adapted from matlab algorithm in
## http://www.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights
############################################################################
## nnodes  : Total number of integration points
## prec    : Floating Point Precision of the result
## verbose : Print things as you iterations proceed? (True or False)
## tol     : Tolerance for L-infinity of difference in old and new guesses
## mxcount : Maximum number of allowed iterations
############################################################################
def lobatto_weights_nodes(nnodes,prec,verbose,tol,mxcount) :
    N1 = nnodes
    N = N1-1
    x = vector( [ cos(pi*(i)/(N)+pi) for i in range(N1) ] )
    x = x.n(prec)
    P = matrix( N1 , N1 , [ 0.n(prec) for i in range(N1^2)] )
    xold = x * 2
    linf = max(x-xold)/max(x)
    count = 0
    while ( linf > tol and count < mxcount ) :
        xold = x
        P[:,0] = 1
        P[:,1] = x
        for k in range(1,N) :
            tmp = matrix(N1,1, [ x[i] * P[i,k] for i in range(N1) ] )
            P[:,k+1] = ( (2*(k+1)-1)*tmp - ((k+1)-1)*P[:,k-1] )/(k+1)
        tmp  = matrix(N1,1, [ x[i] * P[i,N] - P[i,N-1] for i in range(N1) ] )
        tmp2 = matrix(N1,1, [ N1*P[i,N] for i in range(N1) ] )
        tmp3 = vector( [ tmp[i,0] / tmp2[i,0] for i in range(N1) ] )
        x = xold - tmp3
        linf = max(x-xold)/max(x)
        count = count + 1
        if (verbose) :
            print(count,linf)
    w = vector( [ 2 / (N*N1*P[i,N]^2) for i in range(N1) ] )
    return x,w


def tanh_cluster(N,mu) :
    return vector( [ (tanh((2*mu*i)/(N-1)-mu)/tanh(mu)+1)/2-1/2 for i in range(N) ] )


#Returns N GLL point locations in [-0.5,0.5] (quad precision)
def points_gll(N) :
    x,w = lobatto_weights_nodes(N,129,False,1e-35,20)
    return vector([ x[i]/2 for i in range(N) ])
