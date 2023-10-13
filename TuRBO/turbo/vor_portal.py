import rpy2.robjects as robjects
r = robjects.r
r['source']('R/vornoi_cands.R')
import matplotlib.pyplot as plt

import numpy as np

#def vor_cands(X, y, ncands, onlylhs = False, lb = None, ub = None):
#    X_arr = robjects.FloatVector(X.T.flatten())
#    X_R = robjects.r['matrix'](X_arr, nrow = X.shape[0])
#
#    assert (lb is None and ub is None) or (lb is not None and ub is not None)
#
#    y_R = robjects.FloatVector(y.T.flatten())
#
#    if onlylhs:
#        st = 'lhs'
#    else:
#        st = 'rect' if X.shape[0]%2==0 else 'lhs'
#
#    if lb is None:
#        ret = r['vorwalkcands'](X_R, ncands, y_R, st = st, norm = 'linf')
#    else:
#        lb_R = robjects.FloatVector(lb.T.flatten())
#        ub_R = robjects.FloatVector(ub.T.flatten())
#        ret = r['vorwalkcands'](X_R, ncands, y_R, st = st, norm = 'linf', lb = lb_R, ub = ub_R)
#    cands = np.array(ret[0])
#    return cands

def vor_cands(X, y, ncands, onlylhs = False, lb = None, ub = None, pert = None):
    X_arr = robjects.FloatVector(X.T.flatten())
    X_R = robjects.r['matrix'](X_arr, nrow = X.shape[0])

    assert (lb is None and ub is None) or (lb is not None and ub is not None)

    y_R = robjects.FloatVector(y.T.flatten())

    if onlylhs:
        st = 'lhs'
        assert pert is None
    elif pert is not None:
        st_R_arr = robjects.FloatVector(pert.T.flatten())
        st_R = robjects.r['matrix'](st_R_arr, nrow = pert.shape[0])
        st = st_R
    else:
        st = 'rect' if X.shape[0]%2==0 else 'lhs'


    if lb is None:
        ret = r['vorwalkcands'](X_R, ncands, y_R, st = st, norm = 'linf')
    else:
        lb_R = robjects.FloatVector(lb.T.flatten())
        ub_R = robjects.FloatVector(ub.T.flatten())
        ret = r['vorwalkcands'](X_R, ncands, y_R, st = st, norm = 'linf', lb = lb_R, ub = ub_R)
    cands = np.array(ret[0])
    return cands


if __name__=='main':
    N = 10
    P = 2
    X = np.random.uniform(size=[N,P])
    y = X[:,0]
    ncands = 50

    cands = vor_cands(X, y, ncands)

    fig = plt.figure()
    plt.scatter(X[:,0], X[:,1])
    plt.scatter(cands[:,0], cands[:,1])
    plt.savefig("temp.pdf")
    plt.close()