#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python

import sympy

#Single scalar value with a single array of coefficients
def c_scalar(retlab,val,coeflab) :
    import re
    code = ""
    s = retlab+"="
    #Remove spaces and convert to C code
    s = s+str( sympy.ccode(val) ).replace(' ','')
    #Replace coeficients (e.g., c1) with parentheses (e.g., c(1))
    s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
    #Add new line to the end.
    code = code + s + ";\n"
    return code


#Single-dimensional array with a single array of coefficients
def c_vector(retlab,N,vec,coeflab) :
    import re
    code = ""
    for i in range(N) :
        s = retlab+"("+str(i)+")="
        s = s+str( sympy.ccode(vec[i]) ).replace(' ','')
        s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
        code = code + s + ";\n"
    return code


#Two-dimensional array with a single array of coefficients
def c_matrix(retlab,M,N,mat,coeflab) :
    import re
    code = ""
    for j in range(N) :
        for i in range(M) :
            s = retlab+"("+str(j)+","+str(i)+")="
            s = s+str( sympy.ccode(mat[i,j]) ).replace(' ','')
            s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
            code = code + s + ";\n"
    return code


#Three-dimensional array with a single array of coefficients
def c_3d(retlab,N1,N2,N3,mat,coeflab) :
    import re
    code = ""
    for k in range(N3) :
        for j in range(N2) :
            for i in range(N1) :
                s = retlab+"("+str(k)+","+str(j)+","+str(i)+")="
                s = s+str( sympy.ccode(mat[k][j][i]) ).replace(' ','')
                s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
                code = code + s + ";\n"
    return code


#Two-dimensional matrix (array of arrays) with a single array of coefficients
def c_matrix_aoa(retlab,N1,N2,mat,coeflab) :
    import re
    code = ""
    for j in range(N2) :
        for i in range(N1) :
            s = retlab+"("+str(j)+","+str(i)+")="
            s = s+str( sympy.ccode(mat[j][i]) ).replace(' ','')
            s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
            code = code + s + ";\n"
    return code


#Add N spaces to the beginning of a block of code (string)
def add_spaces(N,code) :
    import re
    for i in range(N) :
      code = re.sub("([^\n]*)\n"," \\1\n",code,0,re.DOTALL)
    return code


#Force the 'expr' to be in floating point with the given precision
def force_fp(expr,prec) :
    R = RealField( prec )
    expr = expr * R(2.) / R(2.)
    return expr
