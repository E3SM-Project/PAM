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
#Right interface upwind values
ru_m = ru[1]
ru_p = ru[2]
p_m  =  p[1]
p_p  =  p[2]
w1 = ru_p/2 - p_p/(2*cs)
w2 = ru_m/2 + p_m/(2*cs)
ru_R_upw =   w2 + w1
print("Mass Flux L: coefs p':",jacobian(ru_L_upw,p ).rows())
print("Mass Flux L: coefs ru:",jacobian(ru_L_upw,ru).rows())
print("Mass Flux R: coefs p':",jacobian(ru_R_upw,p ).rows())
print("Mass Flux R: coefs ru:",jacobian(ru_R_upw,ru).rows())
print("")



#########################################
## ORDER 3
#########################################
var('dx,cs')
print("*******************************")
print("** ORDER 3")
print("*******************************")
ru  = coefs_1d(5,0,'ru' )
p   = coefs_1d(5,0,'p'  )
wLm = coefs_1d(3,0,'wLm')
wLp = coefs_1d(3,0,'wLp')
wRm = coefs_1d(3,0,'wRm')
wRp = coefs_1d(3,0,'wRp')
# x-direction
# Left interface upwind values
ru_m = wLm[0]*ru[0] + wLm[1]*ru[1] + wLm[2]*ru[2]
ru_p = wLp[0]*ru[1] + wLp[1]*ru[2] + wLp[2]*ru[3]
p_m  = wLm[0]* p[0] + wLm[1]* p[1] + wLm[2]* p[2]
p_p  = wLp[0]* p[1] + wLp[1]* p[2] + wLp[2]* p[3]
w1 = ru_p/2 - p_p/(2*cs)
w2 = ru_m/2 + p_m/(2*cs)
ru_L_upw =   w2 + w1
#Right interface upwind values
ru_m = wRm[0]*ru[1] + wRm[1]*ru[2] + wRm[2]*ru[3]
ru_p = wRp[0]*ru[2] + wRp[1]*ru[3] + wRp[2]*ru[4]
p_m  = wRm[0]* p[1] + wRm[1]* p[2] + wRm[2]* p[3]
p_p  = wRp[0]* p[2] + wRp[1]* p[3] + wRp[2]* p[4]
w1 = ru_p/2 - p_p/(2*cs)
w2 = ru_m/2 + p_m/(2*cs)
ru_R_upw =   w2 + w1
print("Mass Flux L: coefs p':",jacobian(ru_L_upw,p ).rows())
print("Mass Flux L: coefs ru:",jacobian(ru_L_upw,ru).rows())
print("Mass Flux R: coefs p':",jacobian(ru_R_upw,p ).rows())
print("Mass Flux R: coefs ru:",jacobian(ru_R_upw,ru).rows())
print("")

