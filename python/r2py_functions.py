import rpy2.robjects as robjects
r = robjects.r
import numpy as np

r['source']('R/functions/R_functions.R')
r['source']("R/functions/pomp10.R")
r['source']("R/functions/dacca_cholera.R")

if func=='push':
    f = push_py2R
    P = 14
elif func=='lunar':
    f = lunar
    P = 12
elif func=='rover':
    f = rover_py2R
    P = 60
elif func[:len('ackley')]=='ackley':
    P = int(func[len('ackley'):])
    r('ackley_u <- runif('+str(P)+',-0.5,0.5)')

    #x = np.random.uniform(size=P)
    def f(x):
        x_R = robjects.FloatVector(x)
        return r['ackley'](x_R)[0]
elif func[:len('levy')]=='levy':
    P = int(func[len('levy'):])

    #x = np.random.uniform(size=P)
    def f(x):
        x_R = robjects.FloatVector(x)
        return r['levy'](x_R)[0]
elif func[:len('rosen')]=='rosen':
    P = int(func[len('rosen'):])

    #x = np.random.uniform(size=P)
    def f(x):
        x_R = robjects.FloatVector(x)
        return r['rosen'](x_R)[0]
elif func=='pomp10log':
    P = 10
    def f(x):
        x_R = robjects.FloatVector(x)
        return r['pomp10_R_log'](x_R)[0]
elif func=='dacca':
    P = 23
    def f(x):
        x_R = robjects.FloatVector(x)
        return r['dacca_R'](x_R)[0]
else:
    raise Exception("Unknown function "+str(func))

#def 