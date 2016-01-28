# Plot vertical slice with trajectory postitions, Theta_e

import sys
import traj_tools as trj
import cosmo_utils.pywgrib as pwg
from cosmo_utils.plot import fig_contourf_vertslice_2sp
from cosmo_utils.plot import fig_vertslice_1sp
import cosmo_utils.plot as cplt
import matplotlib.pyplot as plt
from cosmo_utils.helpers  import latlon_distance
import numpy as np
import cosmo_utils as cu
import matplotlib as mpl

mpl.rcParams['font.size'] = 8
mpl.rcParams['font.family'] = 'sans-serif'

which = sys.argv[1]
print which
# Directory for figures
plotdir = '/usr/users/stephan.rasp/repositories/rasp_selz_craig_2015_paper/figures/'
#plotdir = '/usr/users/stephan.rasp/tmp/paper_test/'
cross_flag = True

if which == 'JUL':
    time = 2100
    obj = trj.loadme('c1new.trj')
    lon1 = -4
    lon2 = 1.5
    lat = 0
    dotlon = -1000
    dotlat = 0
    plevs = np.arange(300, 360, 2)
    klim = 10
    fn1 = plotdir + 'JUL_cross'
    fn2 = plotdir + 'JUL_RAP_loc_'
    filtername = ['WCB_Cy1']
    cplot = ['black']
    norm = plt.Normalize(0, 0.1)
    setting = 6
    lllimit = (100, 200) 
    urlimit = (800, 746)

if which == 'JAN':
    time = 1080
    obj = trj.loadme('c2.trj')
    lon1 = -10
    lon2 = -4
    lat = -18
    dotlon = -1000
    dotlat = 0
    plevs = np.arange(300, 330, 2)
    klim = 10
    fn1 = plotdir + 'JAN_cross'
    fn2 = plotdir + 'JAN_RAP_loc_'
    filtername = ['WCB_Cy1']
    cplot = ['black']
    norm = plt.Normalize(0, 0.1)
    setting = 6
    lllimit = (0, 0) 
    urlimit = (1250, 885)

if which == 'OCTnc':
    cross_flag = False
    time = 2760
    obj = trj.loadme('c0_rot.trj')
    lon1 = -1000
    lon2 = -1001
    lat = -13.5
    dotlon = -6.5
    dotlat = -13.5
    plevs = np.arange(300, 330, 2)
    klim = 10
    fn1 = plotdir + 'OCTnc_cross'
    fn2 = plotdir + 'OCTnc_RAP_loc_'
    filtername = ['Wnonconv2']
    cplot = ['black']
    norm = plt.Normalize(0, 0.1)
    lllimit = (400, 250) 
    urlimit = (1400, 900)
    
if which == 'OCTc':
    cross_flag = False
    time = 3600 - 180
    obj = trj.loadme('c0_rot.trj')
    lon1 = -1000
    lon2 = -1001
    lat = -14
    dotlon = 0
    dotlat = -14
    plevs = np.arange(295, 345, 2)
    klim = 10
    fn1 = plotdir + 'OCTc_cross'
    fn2 = plotdir + 'OCTc_RAP_loc_'
    filtername = ['Wconv2']
    cplot = ['black']
    norm = plt.Normalize(0, 0.1)
    lllimit = (400, 250) 
    urlimit = (1400, 900)
  


if which == 21:
    time = 1080
    obj = trj.loadme('c2.trj')
    lon1 = -10
    lon2 = -4
    lat = -18
    dotlon = -1000
    plevs = np.arange(300, 330, 2)
    klim = 10
    fn1 = '/usr/users/stephan.rasp/tmp/paper_test/c0_w2_'+ str(time)
    fn2 = '/usr/users/stephan.rasp/tmp/paper_test/c0_z3_'
    filtername = ['WCB_Cy1_new']
    cplot = ['black']
    norm = plt.Normalize(0, 0.1)
    setting = 5

if which == 22:
    time = 1440
    obj = trj.loadme('c2.trj')
    lon1 = -8
    lon2 = 0
    lat = -17
    dotlon = -1000
    plevs = np.arange(280, 350, 2)
    klim = 10
    fn1 = '/usr/users/stephan.rasp/tmp/paper_test/c0_w2_'+ str(time)
    fn2 = '/usr/users/stephan.rasp/tmp/paper_test/c0_z3_'
    filtername = ['WCB_Cy1_new']
    cplot = ['black']
    norm = plt.Normalize(0, 0.1)
    setting = 5

if which == 23:
    time = 1800
    obj = trj.loadme('c2.trj')
    lon1 = -8
    lon2 = 0
    lat = -17
    dotlon = -1000
    plevs = np.arange(280, 350, 2)
    klim = 10
    fn1 = '/usr/users/stephan.rasp/tmp/paper_test/c0_w2_'+ str(time)
    fn2 = '/usr/users/stephan.rasp/tmp/paper_test/c0_z3_'
    filtername = ['WCB_Cy1_new']
    cplot = ['black']
    norm = plt.Normalize(0, 0.1)
    setting = 5

if which == 24:
    time = 900
    obj = trj.loadme('c2.trj')
    lon1 = -8
    lon2 = 0
    lat = -17
    dotlon = -1000
    plevs = np.arange(280, 350, 2)
    klim = 10
    fn1 = '/usr/users/stephan.rasp/tmp/paper_test/c0_w2_'
    fn2 = '/usr/users/stephan.rasp/tmp/paper_test/c0_z3_'
    filtername = ['WCB_Cy1_new']
    cplot = ['black']
    norm = plt.Normalize(0, 0.1)
    setting = 5

elif which == 1:
    time = 2100
    obj = trj.loadme('c1new.trj')
    lon1 = -4.5
    lon2 = 0
    lat = 0
    dotlon = 0
    plevs = np.arange(280, 380, 2)
    klim = 10
    fn1 = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/skew/200415g_'
    fn2 = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/skew/200415h_slice_'
    filtername = ['WCB_Cy1']
    cplot = ['black']
    norm = plt.Normalize(0, 0.3)
    
elif which == 0:
    time = 3960
    obj = trj.loadme('c0_rot.trj')
    lon1 = -4
    lon2 = 4
    lat = -15
    dotlon = 3
    plevs = np.arange(280, 380, 2)
    klim = 10
    fn1 = '/home/users/stephan.rasp/Dropbox/figures/Case0_20121013/skew/200415a_'
    fn2 = '/home/users/stephan.rasp/Dropbox/figures/Case0_20121013/skew/200415b_slice'
    filtername = ['WCB']
    cplot = ['black']
    norm = plt.Normalize(0, 0.3)

elif which == 4:
    time = 2880
    obj = trj.loadme('c0_rot.trj')
    lon1 = -10
    lon2 = -2
    lat = -13
    dotlon = -3
    plevs = np.arange(280, 380, 2)
    klim = 10
    fn1 = '/home/users/stephan.rasp/Dropbox/figures/Case0_20121013/skew/200415c_'
    fn2 = '/home/users/stephan.rasp/Dropbox/figures/Case0_20121013/skew/200415d_slice'
    filtername = ['WCB']
    cplot = ['black']
    norm = plt.Normalize(0, 0.1)

dt = 15
ddeg = 0.025


# Start
if cross_flag:

    ind, flist = obj._get_index('T', time)
    fn = flist[ind]
    cfn = obj.cfile

    # Get theta_e object
    T = pwg.getfobj(fn, 'T')
    HH = pwg.getfobj(cfn, 'HH')
    HHdata = cu.derive.hl_to_ml(HH.data)
    try:
        P = pwg.getfobj(fn, 'PS').data
    except:
        P0 = trj.utils.cosmo_ref_p(HHdata)
        #print P0.shape
        #ATTENTION: Quick fix, not correct!!!
        tmp = np.ones((P0.shape[0], P0.shape[1], P0.shape[2] + 1))
        tmp[:, :, 0:-1] = P0
        tmp[:, :, -1] = P0[:, :, -1]
        P0 = tmp
        PP = pwg.getfobj(fn, 'PP').data
        P = P0 + PP
        del P0, PP

    QV = pwg.getfobj(fn, 'QV')
    TH_E = T

    mix = QV.data / (1 - QV.data) * 1000.
    del QV
    emat = (P / 100.) * mix / (622. + mix)
    tlclmat = 2840. / (3.5 * np.log(T.data) - np.log(emat) - 4.805) + 55.
    del emat
    tmpthetae = (T.data * (1.e5 / P) ** (0.2854 * 1 - 0.28e-3 *mix) * 
                np.exp((3.376 / tlclmat - 0.00254) * mix * (1 +0.81e-3 * mix)))
    TH_E.data = tmpthetae
    del tmpthetae, tlclmat, mix


    # plot vertical slice
    fig = plt.figure(figsize = (95./25.4, 2.8))
    ax = fig.add_subplot(111)

    lon1ind = np.where(T.rlons[0, :] == lon1)[0][0]
    lon2ind = np.where(T.rlons[0, :] == lon2)[0][0]
    latind = np.where(T.rlats[:, 0] == lat)[0][0]

    contour, cinds, coos = cplt.ax_vertslice(ax, TH_E, ji0 = (latind, lon1ind), 
                                            ji1 = (latind, lon2ind), 
                                hh = HH, pllevels=plevs, 
                                k0k1 = (klim, 49), sp_title = which, cmap = None,
                                annotate_ends = False)

    cb = fig.colorbar(contour, pad = 0.05)
    #cb.set_label(r'$\theta_e$ [K]')
    plt.text(1.17, -0.07, r'$\theta_e$', transform = ax.transAxes)
    plt.text(1.165, -0.145, r'[K]', transform = ax.transAxes)
    #plt.title(which)
    plt.xlim(0, 600)
    plt.ylim(0, 10000)
    # Get trj positions

    for filt in filtername:
        zarray, lonarray, carray = obj.get_cross_pos(lon1, lon2, lat, time - dt, time + dt, 
                                            filt, cname = 'P400withinP600')
        #zarray = zarray[::5]
        #lonarray = lonarray[::5]

        distlist = []
        for i in range(lonarray.shape[0]):
            distlist.append(latlon_distance((lat, lon1), (lat, lonarray[i])))

        distarray = np.array(distlist) / 1000.
        carray = 400. / carray / 60.
        
        sca = ax.scatter(distarray, zarray, color = 'black', s = 7, zorder = 2, 
                        linewidth = 0.1,
                cmap = 'hot', norm = norm) 
        #cb2 = plt.colorbar(sca, pad = 0.02)
        #cb2.set_label(r'$\bar{\omega}$ [hPa s$^{-1}$]')
        #plt.text(1.03, -0.07, r'$\bar{\omega}$', transform = ax.transAxes)
        #plt.text(0.98, -0.145, r'[-hPa s$^{-1}$]', transform = ax.transAxes)
    #plt.tight_layout()
    plt.savefig(fn1, bbox_inches = 'tight', dpi = 400)


# Plot contour with line
obj.draw_trj_dot(["TOT_PREC_S", 'var145_S'], tplus = time, 
                filtername=filtername[0], onlyasc='P600', ctracer = 'P_dt', 
                thinning = 1, setting = 6, cbar = False, cross = [dotlon, dotlat], 
                line = [lon1, lat, lon2, lat], savebase = fn2, 
                lllimit = lllimit, urlimit = urlimit, idtext = which)





