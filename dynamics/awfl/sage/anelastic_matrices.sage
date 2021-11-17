def coefs_1d(N,N0,lab) :
    return vector([ var(lab+'%s'%i) for i in range(N0,N0+N) ])

def poly_1d(N,coefs,x) :
    return sum( vector([ coefs[i]*x^i for i in range(N) ]) )

def print_vector(M):
    print(M.rows())

#########################################
## ORDER 1
#########################################
var('dx,cs')
print("*******************************")
print("** ORDER 1")
print("*******************************")
ru_old = coefs_1d(5,0,'ru_old')
ru = coefs_1d(5,0,'ru')
p  = coefs_1d(5,0,'p' )
# Left interface upwind values
ru_m = ru[1]
ru_p = ru[2]
p_m  =  p[1]
p_p  =  p[2]
w1 = ru_p/2 - p_p/(2*cs)
w2 = ru_m/2 + p_m/(2*cs)
ru_L_upw = w2 + w1
#Right interface upwind values
ru_m = ru[2]
ru_p = ru[3]
p_m  =  p[2]
p_p  =  p[3]
w1 = ru_p/2 - p_p/(2*cs)
w2 = ru_m/2 + p_m/(2*cs)
ru_R_upw = w2 + w1


ru_L_upw = ru_L_upw.subs( ru1=ru_old1-(p[2]-p[0])/(2*dx) ,
                          ru2=ru_old2-(p[3]-p[1])/(2*dx) ,
                          ru3=ru_old3-(p[4]-p[2])/(2*dx) ).expand()

ru_R_upw = ru_R_upw.subs( ru1=ru_old1-(p[2]-p[0])/(2*dx) ,
                          ru2=ru_old2-(p[3]-p[1])/(2*dx) ,
                          ru3=ru_old3-(p[4]-p[2])/(2*dx) ).expand()

mom_div = ( (ru_R_upw - ru_L_upw) / dx ).expand()


pcoefs_L_upw  = jacobian(ru_L_upw,p)
pcoefs_R_upw  = jacobian(ru_R_upw,p)
pcoefs_momdiv = jacobian(mom_div ,p)

print("Mass Flux L: coefs p'    : ",pcoefs_L_upw .rows())
print("Mass Flux R: coefs p'    : ",pcoefs_R_upw .rows())
print("Mom. Div.  : coefs p'    : ",pcoefs_momdiv.rows())
print("RHS        : coefs ru_old: ",jacobian( -mom_div , ru_old ).rows() )


