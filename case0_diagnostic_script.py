# Script for filtering data and claculating diagnostics

import traj_tools as trj
import numpy as np
import matplotlib.pyplot as plt
import cPickle

# Start script
f = open('/home/scratch/users/stephan.rasp/traj_data/IndexMatrix')
indmat = cPickle.load(f)
filelist = trj.common.CreateFileArray('/home/scratch/users/stephan.rasp/traj_data/test_')
traceind = 7   # 'P'
criterion1 = 600   # hPa
c1asc = trj.filters.MinXMatrix(filelist, traceind, criterion1, 
                                IndMatrix = indmat, Flat = True)
criterion2 = 100
c2asc = trj.filters.MinXMatrix(filelist, traceind, criterion2, 
                                IndMatrix = indmat, Flat = True)