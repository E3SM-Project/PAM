#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python


#Generates a set of coefficients as a sage "vector"
def coefs_1d(N,N0,lab) :
    return vector([ var(lab+'%s'%i) for i in range(N0,N0+N) ])


#Generates a 1-D polynomial as a sage expreession based on a set of coefficients
def poly_1d(N,coefs,x) :
    return sum( vector([ coefs[i]*x^i for i in range(N) ]) )


#Computes the coefficients of an Nth-order polynomial expression and returns a sage vector
#Mostly just a helper function for integrate_poly
def compute_coefs(N,p,x) :
    x = var(str(x))
    return vector([ diff(p,x,i).subs({x:0})/factorial(i) for i in range(N) ])


#Computes the integral of an Nth-order polynomial over the range [x1,x2].
#I had to do this myself because Sage's "integral" doesn't allow arbitrary accuracy
def integrate_poly(N,p,x,x1,x2) :
    var('a')
    a = compute_coefs(N,p,x)
    return sum( vector([ a[i]*(x2**(i+1) - x1**(i+1)) / (i+1) for i in range(N) ]) )
