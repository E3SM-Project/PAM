#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python

#Single scalar value with a single array of coefficients
def c_scalar(retlab,val,coeflab,wraplen) :
    import re
    code = ""
    s = retlab+"="
    #Remove spaces
    s = s+str(val).replace(' ','')
    s = re.sub("([a-zA-Z()0-9_]*)\^2","(\\1*\\1)",s,0,re.DOTALL)
    #Replace exponentials with power functions
    s = re.sub("([a-zA-Z()0-9_]*)\^([0-9]*)","pow((double)\\1,(double)\\2)",s,0,re.DOTALL)
    #Replace coeficients (e.g., c1) with parentheses (e.g., c(1))
    s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
    #Break up code to the given wraplen for line length
    #s = re.sub("&\n&$","",re.sub("(.{"+str(wraplen)+"})", "\\1&\n&", s, 0, re.DOTALL))
    #Add new line to the end.
    code = code + s + ";\n"
    return code


#Single-dimensional array with a single array of coefficients
def c_vector(retlab,N,vec,coeflab,wraplen) :
    import re
    code = ""
    for i in range(N) :
        s = retlab+"("+str(i)+")="
        s = s+str(vec[i]).replace(' ','')
        s = re.sub("([a-zA-Z()0-9_]*)\^2","(\\1*\\1)",s,0,re.DOTALL)
        s = re.sub("([a-zA-Z()0-9_]*)\^([0-9]*)","pow((double)\\1,(double)\\2)",s,0,re.DOTALL)
        s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
        #s = re.sub("&\n&$","",re.sub("(.{"+str(wraplen)+"})", "\\1&\n&", s, 0, re.DOTALL))
        code = code + s + ";\n"
    return code


#Two-dimensional array with a single array of coefficients
def c_matrix(retlab,M,N,mat,coeflab,wraplen) :
    import re
    code = ""
    for j in range(N) :
        for i in range(M) :
            s = retlab+"("+str(j)+","+str(i)+")="
            s = s+str(mat[i,j]).replace(' ','')
            s = re.sub("([a-zA-Z()0-9_]*)\^2","(\\1*\\1)",s,0,re.DOTALL)
            s = re.sub("([a-zA-Z()0-9_]*)\^([0-9]*)","pow((double)\\1,(double)\\2)",s,0,re.DOTALL)
            s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
            #s = re.sub("&\n&$","",re.sub("(.{"+str(wraplen)+"})", "\\1&\n&", s, 0, re.DOTALL))
            code = code + s + ";\n"
    return code


#Three-dimensional array with a single array of coefficients
def c_3d(retlab,N1,N2,N3,mat,coeflab,wraplen) :
    import re
    code = ""
    for k in range(N3) :
        for j in range(N2) :
            for i in range(N1) :
                s = retlab+"("+str(k)+","+str(j)+","+str(i)+")="
                s = s+str(mat[k][j][i]).replace(' ','')
                s = re.sub("([a-zA-Z()0-9_]*)\^2","(\\1*\\1)",s,0,re.DOTALL)
                s = re.sub("([a-zA-Z()0-9_]*)\^([0-9]*)","pow((double)\\1,(double)\\2)",s,0,re.DOTALL)
                s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
                #s = re.sub("&\n&$","",re.sub("(.{"+str(wraplen)+"})", "\\1&\n&", s, 0, re.DOTALL))
                code = code + s + ";\n"
    return code


#Two-dimensional matrix (array of arrays) with a single array of coefficients
def c_matrix_aoa(retlab,N1,N2,mat,coeflab,wraplen) :
    import re
    code = ""
    for j in range(N2) :
        for i in range(N1) :
            s = retlab+"("+str(j)+","+str(i)+")="
            s = s+str(mat[j][i]).replace(' ','')
            s = re.sub("([a-zA-Z()0-9_]*)\^2","(\\1*\\1)",s,0,re.DOTALL)
            s = re.sub("([a-zA-Z()0-9_]*)\^([0-9]*)","pow((double)\\1,(double)\\2)",s,0,re.DOTALL)
            s = re.sub(coeflab+"([0-9]+)",coeflab+"(\\1)",s, 0, re.DOTALL)
            #s = re.sub("&\n&$","",re.sub("(.{"+str(wraplen)+"})", "\\1&\n&", s, 0, re.DOTALL))
            code = code + s + ";\n"
    return code


#Add N spaces to the beginning of a block of code (string)
def add_spaces(N,code) :
    import re
    code = re.sub("([^\n]*)\n","    \\1\n",code,0,re.DOTALL)
    return code


#Force the 'expr' to be in floating point with the given precision
def force_fp(expr,prec) :
    R = RealField( prec )
    expr = expr * R(2.) / R(2.)
    return expr
