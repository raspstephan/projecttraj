# Script for filtering data and calculating diagnostics

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
(nconv, start100, end100) = trj.filters.FilterStartEnd(filelist, traceind, 400, 576, 24)

# Merge index matrices
newindmat = []
newstartmat = []
newendmat = []
for i in range(len(nconv)):
    newindmat.append([])
    newstartmat.append([])
    newendmat.append([])
    for j in range(len(nconv[i])):
        if nconv[i][j] in indmat[i]:
            newindmat[i].append(nconv[i][j])
            newstartmat[i].append(start100[i][j])
            newendmat[i].append(end100[i][j])


#c1asc = trj.filters.MinXMatrix(filelist, traceind, criterion1, 
                                #IndMatrix = indmat, IndMatrix2 = nconv,
                                #Flat = True)
#criterion2 = 100
#c2asc = trj.filters.MinXMatrix(filelist, traceind, criterion2, 
                                #IndMatrix = indmat, IndMatrix2 = nconv,
                                #Flat = True)
#dt = 5
#c1asc = np.array(c1asc) * dt
#c2asc = np.array(c2asc) * criterion1 / criterion2 * dt


# Plot XY plot for all WCB
#trj.plots.DrawXYSingle([], 360, 
                       #SaveBase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/XY_start360_classicWCB',
                       #title = '260814b')
# Plot XY for non-convective
trj.plots.DrawXYSingle([], 360, newindmat,
                       SaveBase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/XY_start360_nonconv',
                       title = '260814c')

# Plot XY only of 100hPa ascent

trj.plots.DrawXYSingle([], 360, newindmat,
                       SaveBase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/XY_start360_onlyasc',
                       title = '260814d', only = True, StartM = newstartmat, EndM = newendmat)

## Plot Scatter
#plt.close('all')
##plt.figure(figsize = (10, 10))
#plt.scatter(c1asc, c2asc)
#plt.plot([0, 576 * dt], [0, 576 * dt])
#plt.xlim(0, 576 * dt)
#plt.ylim(0, 576 * dt)
#plt.xlabel('Time for 600hPa ascent [mins]')
#plt.ylabel('Time x 6 for 100hPa ascent [mins]')
#plt.text(0.94, 1.02, '260814a', transform = plt.gca().transAxes, fontsize = 6)
#plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/Scatter600v100')
