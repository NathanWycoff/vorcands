#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  python/turbo_example.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 07.07.2022

# TODO: same initialization

import numpy as np
#import rpy2
#import rpy2.robjects.numpy2ri
#rpy2.robjects.numpy2ri.activate()

import sys
import os

import time

from turbo import Turbo1,TurboM
import numpy as np
import torch
import math
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

#exec(open("python/functions.py").read())
exec(open("./python/push.py").read())
exec(open("./python/bnn_fun.py").read())
exec(open("./python/rover.py").read())

func = sys.argv[1]
seed = sys.argv[2]
exec(open("./python/r2py_functions.py").read())

from turbo.utils import from_unit_cube, latin_hypercube, to_unit_cube

sim_name = func + '_' + str(seed) + '.csv'
init_file = 'sim_inits/'+sim_name
init_lhs = np.array(pd.read_csv(init_file).iloc[:,1:])
n_init = init_lhs.shape[0]

result_file = './sim_out/'+func+'/'+sim_name
outdf = pd.read_csv(result_file).iloc[:,1:]
time_file = './sim_out/'+func+'_time/'+sim_name
timedf = pd.read_csv(time_file)

max_evals = outdf.shape[0]

for vornoi in [False, True]:
    tt = time.time()

    turbo1 = Turbo1(
        f=f,  # Handle to objective function
        lb=np.zeros(P),  # Numpy array specifying lower bounds
        ub=np.ones(P),  # Numpy array specifying upper bounds
        #n_init=1000,  # Number of initial bounds from an Latin hypercube design
        #n_init=int(np.max([12,3*P])),  # Number of initial bounds from an Latin hypercube design
        #init_lhs = init_lhs,
        n_init=n_init,  # Number of initial bounds from an Latin hypercube design
        max_evals = max_evals,  # Maximum number of evaluations
        batch_size=1,  # How large batch size TuRBO uses
        #batch_size=50,  # How large batch size TuRBO uses
        verbose=False,  # Print information from each batch
        use_ard=True,  # Set to true if you want to use ARD for the GP kernel
        max_cholesky_size=2000,  # When we switch from Cholesky to Lanczos
        n_training_steps=50,  # Number of steps of ADAM to learn the hypers
        min_cuda=1024,  # Run on the CPU for small datasets
        device="cpu",  # "cpu" or "cuda"
        dtype="float64",  # float64 or float32
        vor_tr=vornoi
        #vornoi=vornoi
    )
    #print("Warn: not TR vornoi.")

    turbo1.optimize()

    X = turbo1.X  # Evaluated points
    fX = turbo1.fX  # Observed values
    ind_best = np.argmin(fX)
    f_best, x_best = fX[ind_best], X[ind_best, :]

    ## If necessary, crop to correct size.
    X = X[:max_evals,:]
    fX = fX[:max_evals]

    ## bov computation.
    for i in range(max_evals):
        fX[i] = np.min(fX[:(i+1)])
    fX = np.array(fX).flatten()
    name = 'TuRBO' if not vornoi else 'vor.TuRBO'
    outdf = pd.concat([outdf, pd.Series(fX, name = name)], axis=1)

    # Ellapsed Real Time
    ert = time.time()-tt
    timedf.loc[-1] = [name,ert]
    timedf.index = np.arange(timedf.shape[0])


outdf.to_csv(result_file)
timedf.to_csv(time_file, index = False)

#np.min(fX)


#xxf = X[np.argmin(fX),:]
#print("Our best value:")
#print(f(xxf))
#print("Optimizing params")
#print(xxf)
#print("Benchmark value:")
#print(f(testf['xtest']))
#testf['bt'](xxf)
#
