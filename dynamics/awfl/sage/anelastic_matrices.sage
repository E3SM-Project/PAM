def coefs_1d(N,N0,lab) :
    return vector([ var(lab+'%s'%i) for i in range(N0,N0+N) ])

def poly_1d(N,coefs,x) :
    return sum( vector([ coefs[i]*x^i for i in range(N) ]) )

#########################################
## ORDER 1
#########################################
var('dx,cs')
print("*******************************")
print("** ORDER 1")
print("*******************************")
ru = coefs_1d(3,0,'ru')
p  = coefs_1d(3,0,'p' )
# x-direction
# Left interface upwind values
ru_m = ru[0]
ru_p = ru[1]
p_m  =  p[0]
p_p  =  p[1]
w1 = ru_p/2 - p_p/(2*cs)
w2 = ru_m/2 + p_m/(2*cs)
ru_L_upw =   w2 + w1
p_L_upw  = ( w2 - w1 ) * cs
#Right interface upwind values
ru_m = ru[1]
ru_p = ru[2]
p_m  =  p[1]
p_p  =  p[2]
w1 = ru_p/2 - p_p/(2*cs)
w2 = ru_m/2 + p_m/(2*cs)
ru_R_upw =   w2 + w1
p_R_upw  = ( w2 - w1 ) * cs
print("Mass Flux L: coefs p':",jacobian(ru_L_upw,p ))
print("Mass Flux L: coefs ru:",jacobian(ru_L_upw,ru))
print("Mass Flux R: coefs p':",jacobian(ru_R_upw,p ))
print("Mass Flux R: coefs ru:",jacobian(ru_R_upw,ru))
print("Pressure  L: coefs p':",jacobian( p_L_upw,p ))
print("Pressure  L: coefs ru:",jacobian( p_L_upw,ru))
print("Pressure  R: coefs p':",jacobian( p_R_upw,p ))
print("Pressure  R: coefs ru:",jacobian( p_R_upw,ru))
print("")



