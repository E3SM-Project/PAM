def coefs_1d(N,N0,lab) :
    return vector([ var(lab+'%s'%i) for i in range(N0,N0+N) ])

def poly_1d(N,coefs,x) :
    return sum( vector([ coefs[i]*x^i for i in range(N) ]) )

#########################################
## ORDER 1
#########################################
var('dx,dz')
print("*******************************")
print("** ORDER 1")
print("*******************************")
rw = coefs_1d(3,0,'rw')
ru = coefs_1d(3,0,'ru')
p  = coefs_1d(3,0,'p' )
# x-direction
# Left interface upwind values
ru_m = ru[0]
ru_p = ru[1]
p_m  =  p[0]
p_p  =  p[1]
p_x_L_upw  = ru_m/2 + p_m/2 - ru_p/2  + p_p/2
ru_x_L_upw = ru_m/2 + p_m/2 + ru_p/2  - p_p/2
#Right interface upwind values
ru_m = ru[1]
ru_p = ru[2]
p_m  =  p[1]
p_p  =  p[2]
p_x_R_upw  = ru_m/2 + p_m/2 - ru_p/2  + p_p/2
ru_x_R_upw = ru_m/2 + p_m/2 + ru_p/2  - p_p/2
p_expr = (ru_x_R_upw - ru_x_L_upw) / dx
ru_expr = ru[1] + (p_x_R_upw - p_x_L_upw) / dx
print("X-direction")
print("Mass Flux L: coefs p':",jacobian(ru_x_L_upw,p ))
print("Mass Flux L: coefs ru:",jacobian(ru_x_L_upw,ru))
print("Mass Flux R: coefs p':",jacobian(ru_x_R_upw,p ))
print("Mass Flux R: coefs ru:",jacobian(ru_x_R_upw,ru))
print("Pressure  L: coefs p':",jacobian( p_x_L_upw,p ))
print("Pressure  L: coefs ru:",jacobian( p_x_L_upw,ru))
print("Pressure  R: coefs p':",jacobian( p_x_R_upw,p ))
print("Pressure  R: coefs ru:",jacobian( p_x_R_upw,ru))
print("Mom Div    : coefs p':",jacobian( p_expr   ,p ))
print("Mom Div    : coefs ru:",jacobian( p_expr   ,ru))
print("Press Div  : coefs p':",jacobian(ru_expr   ,p ))
print("Press Div  : coefs ru:",jacobian(ru_expr   ,ru))
print("")



