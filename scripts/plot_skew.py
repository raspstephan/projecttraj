# Plot skew T log P

import sys
import traj_tools as trj
from cosmo_utils.sounding import cosmo_sounding
import cosmo_utils.pywgrib as pwg
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import cosmo_utils as cu

matplotlib.rcParams['font.size'] = 8

plotdir = '/usr/users/stephan.rasp/repositories/rasp_selz_craig_2015_paper/figures/'
#plotdir = '/usr/users/stephan.rasp/tmp/paper_test/'

#JUL
which = sys.argv[1]
print which

if which == 1:
    time = 2100

    c1 = trj.loadme('c1new.trj')
    ind, flist = c1._get_index('T', time)


    fn = flist[ind]
    T = pwg.getfobj(fn, 'T')
    lon = -0
    lat = 0
    si = np.where(T.rlons[0, :] == lon)[0][0]
    sj = np.where(T.rlats[:, 0] == lat)[0][0]

    S = cosmo_sounding(fn, hhf = c1.cfile)

    fig = S.plot_skewt(i = int(si), j = int(sj), title='')
    fig.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/skew/200415u_')

    #c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], time, setting = 2, cbar = False, dot = (lon, lat), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/skew/080415r_')

#JAN
if which == 2:
    time = 1140

    c2 = trj.loadme('c2.trj')
    ind, flist = c2._get_index('T', time)


    fn = flist[ind]
    T = pwg.getfobj(fn, 'T')
    lon = -4.5
    lat = -17.5
    si = np.where(T.rlons[0, :] == lon)[0][0]
    sj = np.where(T.rlats[:, 0] == lat)[0][0]

    S = cosmo_sounding(fn, hhf = c2.cfile)

    fig = S.plot_skewt(i = int(si), j = int(sj), title='')
    fig.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/skew/200415v_')

    #c2.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], time, setting = 5, cbar = False, dot = (lon, lat), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/skew/080415t_')

#OCT
if which == 'OCTc' or which == 'OCTnc':
    if which == 'OCTc':
        time = 3600 - 180
        lon = 0
        lat = -14
        fn1 = plotdir + 'OCTc_sounding'
    elif which == 'OCTnc':
        time = 2760 -180
        lon = -6.5
        lat = -13.5
        fn1 = plotdir + 'OCTnc_sounding'

    c0 = trj.loadme('c0_rot.trj')
    ind, flist = c0._get_index('T', time)
    
    fn = flist[ind]
    
    # Get PS
    HH = cu.derive.hl_to_ml(pwg.getfobj(c0.cfile, 'HH').data)
    P0 = trj.utils.cosmo_ref_p(HH) / 100.   # hPa
    #print P0.shape
    #ATTENTION: Quick fix, not correct!!!
    tmp = np.ones((P0.shape[0], P0.shape[1], P0.shape[2] + 1))
    tmp[:, :, 0:-1] = P0
    tmp[:, :, -1] = P0[:, :, -1]
    P0 = tmp
    PP = pwg.getfobj(fn, 'PP').data / 100. # hPa
    PS = PP + P0
    
    T = pwg.getfobj(fn, 'T')
    si = np.where(T.rlons[0, :] == lon)[0][0]
    sj = np.where(T.rlats[:, 0] == lat)[0][0]

    S = cosmo_sounding(fn, hhf = c0.cfile, psext = PS)

    fig = S.plot_skewt(i = int(si), j = int(sj), title=which)
    fig.savefig(fn1)

    #c0.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], time, setting = 2, 
                    #cbar = False, dot = (lon, lat), savebase = fn2)
