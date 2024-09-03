#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: Oliver Rice

Python wrapper for running multiple instances of lare2d, varying certain parameters

"""


import os
import shutil
import numpy as np
import sys
from init import compute_initial_condition
from numpy import random
import time
from scipy.io import netcdf_file

#os.system('killall mf2d')

if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    raise Exception("Provide run ID and try again.")

if len(sys.argv) > 2:
    ncores = int(sys.argv[2])
else:
    ncores = 1

if os.uname()[1][-14:] == 'ham8.dur.ac.uk':
    hamilton_flag = 1
else:
    hamilton_flag = 0

time_unit = 1   #time unit in days. LARE seems to prefer values of around unity so do have to be careful
voutfact = 0.0

nx = 64
ny = 64

x0 = -0.5; x1 = 0.5
y0 = -1.0/ny; y1 = 1.0

shearfact = 1.0
bfact = 1.0

nplots = 50
ndiags = 600
tmax = 150.0/time_unit

eta = 1e-6

nu0 = np.geomspace(0.05,0.2)[5]
eta0 = np.geomspace(7.5e-4,1.25e-3,10)[5]

decay_type = 1   #1 for exponential, 0 for tanh

def lbound_fn(x):
    #Outputs lower boundary radial magnetic field as a function of position x
    return -bfact*np.sin(0.5*np.pi*x/x1) + bfact*0*random.rand(len(x))
    #return random.rand(len(x))

def zero_fn(x):
    return 0.0

if decay_type > 0.5: #exponential decay
    a = 0.1; b = np.linspace(0.0,1.0,10)[0]
    #a = 0.0; b =0.5
    deltay = 1.0; ystar = 1.0

else:   #for the tanh cutoff
    a = 0.25; b = 1.0
    ystar = np.linspace(0.0,0.5,10)[run//50]
    deltay = ystar/10

buoyant_factor = 0.01*run
nu0_decay = 0.0

variables = np.zeros((30))

variables[0] = run
variables[1] = nx
variables[2] = tmax
variables[3] = nplots
variables[4] = 5.0
variables[5] = bfact
variables[6] = shearfact
variables[7] = eta
variables[8] = nu0
variables[9] = eta0
variables[10] = x0
variables[11] = x1
variables[12] = y0
variables[13] = y1

variables[14] = a
variables[15] = b
variables[16] = deltay
variables[17] = ystar

variables[18] = hamilton_flag
variables[19] = ndiags

variables[20] = ny
variables[21] = nx   #init id  (id of initial condition)

variables[22] = nu0_decay

variables[23] = int(decay_type)
variables[24] = buoyant_factor

if hamilton_flag < 0.5:
    string = '/extra/tmp/trcn27/mf2d/%03d' % run
else:
    string = '/nobackup/trcn27/mf2d0/%03d' % run

if not os.path.exists(string[:-4]):
    os.mkdir(string[:-4])


if os.path.exists(string):
    shutil.rmtree(string)

os.mkdir(string)

if True:
    os.system('make')

class Grid():
    def __init__(self):
        self.x0 = x0; self.x1 = x1; self.y0 = y0; self.y1 = y1
        self.nx = nx ; self.ny = ny
        self.xs = np.linspace(self.x0,self.x1,self.nx+1)
        self.ys = np.linspace(self.y0,self.y1,self.ny+1)


if os.path.isdir('%s%03d/' % (string,run)):
    shutil.rmtree('%s%03d/' % (string,run))

time.sleep(2.0)

np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)

#Create initial condition using new init.py (potential field with arbitrary lower boundary and domain dimensions)
compute_initial_condition(Grid(), lbound_fn, run)

print('Output directory "%s" created' % (string))
if hamilton_flag < 0.5:
    os.system('/usr/lib64/openmpi/bin/mpiexec -np %d ./bin/mf2d %d' % (ncores, run))
