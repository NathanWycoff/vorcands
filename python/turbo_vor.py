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

tt = time.time()

#func = sys.argv[1]
#seed = sys.argv[2]
func = 'push'
seed = 1
#func = 'push'
#seed = 1

if func=='push':
    f = push_py2R
    P = 14
elif func=='lunar':
    f = lunar
    P = 12
elif func=='rover':
    f = rover_py2R
    P = 60
else:
    raise Exception("Unknown function.")

from turbo.utils import from_unit_cube, latin_hypercube, to_unit_cube

#sim_name = func + '_' + str(seed) + '.csv'
#init_file = 'sim_inits/'+sim_name
#init_lhs = np.array(pd.read_csv(init_file).iloc[:,1:])
#n_init = init_lhs.shape[0]

xx = np.random.uniform(size=P)
print(f(xx))

n_init = 20
#max_evals = 200
max_evals = 800

turbo1 = Turbo1(
    f=f,  # Handle to objective function
    lb=np.zeros(P),  # Numpy array specifying lower bounds
    ub=np.ones(P),  # Numpy array specifying upper bounds
    #n_init=1000,  # Number of initial bounds from an Latin hypercube design
    #n_init=int(np.max([12,3*P])),  # Number of initial bounds from an Latin hypercube design
    n_init=n_init,  # Number of initial bounds from an Latin hypercube design
    #max_evals = max_evals,  # Maximum number of evaluations
    max_evals = n_init+1,  # Maximum number of evaluations
    batch_size=1,  # How large batch size TuRBO uses
    verbose=True,  # Print information from each batch
    use_ard=True,  # Set to true if you want to use ARD for the GP kernel
    max_cholesky_size=2000,  # When we switch from Cholesky to Lanczos
    n_training_steps=50,  # Number of steps of ADAM to learn the hypers
    min_cuda=1024,  # Run on the CPU for small datasets
    device="cpu",  # "cpu" or "cuda"
    vornoi_tr=True,
    dtype="float64",  # float64 or float32
)

turbo1.n_cand =10
turbo1.optimize()

X = turbo1.X  # Evaluated points
fX = turbo1.fX  # Observed values
ind_best = np.argmin(fX)
f_best, x_best = fX[ind_best], X[ind_best, :]

ert = time.time()-tt
print("elapsed real time:")
print(ert)

xxf = X[np.argmin(fX),:]
print("Our best value:")
print(f(xxf))
#print("Optimizing params")
#print(xxf)
#print("Benchmark value:")
#print(f(testf['xtest']))
#testf['bt'](xxf)
#

fig = plt.figure()
plt.plot([np.sum(x) for x in turbo1.isins])
plt.ylabel("Design Points in Trust Region")
ax = plt.gca()
ax1 = ax.twinx()
ax1.plot(turbo1.lengths, color = 'orange')
ax1.set_ylabel("TR Size Coef")
plt.title("TR status on " + func)
plt.savefig("tr_size.pdf")
plt.close()
