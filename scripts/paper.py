"""
Script to plot figures in Paper
"""
import sys
import os
# Import trajectory tools
import traj_tools as trj
import matplotlib.pyplot as plt
import numpy as np
import cPickle

# Load saved trajectory objects
path = '/home/users/stephan.rasp/repositories/projecttraj/scripts/'
c0 = trj.loadme(path + 'c0_rot.trj')
c1 = trj.loadme(path + 'c1new.trj')
c2 = trj.loadme(path + 'c2.trj')

# Directory for figures
plotdir = '/usr/users/stephan.rasp/repositories/rasp_selz_craig_2015_paper/figures/'
#plotdir = '/usr/users/stephan.rasp/tmp/paper_test/'

# Determine which figures to plot
plotlist = []
choicelist = ['spaghetti', 'RAP_hist', 'ascent_hist', 'trj_path', 'cross_skew',
              'temporal', 'pv_hist', 'vort_hist', 'theta_hist', 'temporal_JANJUL', 
              'pv_pdf', 'pv_vs_t', 'theta_pdf', 'theta_vs_t']
for arg in sys.argv[1:]:
    if arg in choicelist:
        plotlist.append(arg)
if len(plotlist) == 0:
    print 'Please chose a plot to be plotted. Choices are: '
    print choicelist

# plot spaghetti plots for all 3 cases
if 'spaghetti' in plotlist:
    c1.draw_spaghetti('P', 'WCB_Cy1', limit = [1, 12, 13, 4, 5, 6, 7, 18, 19, 20], 
                      locind = 20, 
                      xlim = (20, 60), savebase = plotdir + 'JUL_', 
                      title = 'JUL')

    c0.draw_spaghetti('P', 'Wconv2', limit = [1, 12, 13, 4, 5, 6, 7, 18, 19, 20],
                      locind = 51, 
                      xlim = (45, 85), savebase = plotdir + 'OCTc_',
                      title = 'OCTc')
    c0.draw_spaghetti('P', 'Wnonconv2', limit = [1, 12, 13, 4, 5, 6, 7, 18, 19, 20], 
                      locind = 30, 
                      xlim = (45, 85), savebase = plotdir + 'OCTnc_',
                      title = 'OCTnc', legend = False)

    c2.draw_spaghetti('P', 'WCB_Cy1_new', limit = [1, 12, 13, 4, 5, 6, 7, 18, 19, 20],
                      locind = 1, 
                      xlim = (5, 45), savebase = plotdir + 'JAN_',
                      title = 'JAN', legend = False)


# Plot RAP 2D histograms
if 'RAP_hist' in plotlist:
    c1.draw_slice_hist('P_dt', 0.1, 'WCB_Cy1', onlyasc = 'P600', mult = -1, 
                       savebase = plotdir + 'JUL_', title = 'JUL')
    
    c0.draw_slice_hist('P_dt', 0.1, 'Wconv2', onlyasc = 'P600', mult = -1, 
                       savebase = plotdir + 'OCTc_', title = 'OCTc')
    c0.draw_slice_hist('P_dt', 0.1, 'Wnonconv2', onlyasc = 'P600', mult = -1, 
                       savebase = plotdir + 'OCTnc_', title = 'OCTnc')

    c2.draw_slice_hist('P_dt', 0.1, 'WCB_Cy1_new', onlyasc = 'P600', mult = -1, 
                       savebase = plotdir + 'JAN_', title = 'JAN')
    
# Plot ascent histograms
if 'ascent_hist' in plotlist:
    c2.draw_hist('P600', (0, 48), (0, 250),filtername = 'WCB_Cy1_new',
                 bins = 96, ylog = False, idtext = 'JAN', mintohrs = True, 
                 savebase = plotdir + 'JAN_')
    c1.draw_hist('P600', (0, 48), (0, 5000), filtername = 'WCB_Cy1',
                 bins = 96, ylog = False, idtext = 'JUL', mintohrs = True, 
                 savebase = plotdir + 'JUL_')

    c0.draw_hist_stacked(['P600', 'P600'], ['Wconv2', 'Wnonconv2'], 
                         ['OCTc', 'OCTnc'], xlim = (0, 48), ylim = (0, 1400), 
                         bins = 96, ylog = False, idtext = '', 
                         savebase = plotdir + 'OCT_')

# Plot OCT trajectory path
if 'trj_path' in plotlist:
    c0.draw_trj_all([], 'Wconv2', thinning = 20, savebase=plotdir + 'OCTc_', 
                    title = 'OCTc', intersect = 750, onlyasc = 'P600',
                    lllimit = (300, 050), urlimit = (1900, 1400))
    c0.draw_trj_all([], 'Wnonconv2', thinning = 20, savebase=plotdir + 'OCTnc_', 
                    title = 'OCTnc', intersect = 750, onlyasc = 'P600',
                    lllimit = (300, 050), urlimit = (1900, 1400))

if 'cross_skew' in plotlist:
    os.system('python plot_skew.py OCTc')
    os.system('python plot_vert_slice.py OCTc')
    
    #os.system('python plot_skew.py OCTnc')
    #os.system('python plot_vert_slice.py OCTnc')
    
    #os.system('python plot_vert_slice.py JAN')

if 'temporal' in plotlist:
    conv = (c0._mask_array('Wconv2', 'Pcross750') * 5. + 
            c0._mask_array('Wconv2', 'startt'))
    nonconv = (c0._mask_array('Wnonconv2', 'Pcross750') * 5. + 
               c0._mask_array('Wnonconv2', 'startt'))
    fig = plt.figure(figsize = (95./25.4, 3))
    num1, xpos1 = np.histogram(conv, bins = 120, range = [0, 7200] )
    num2, xpos2 = np.histogram(nonconv, bins = 120, range = [0, 7200] )
    x = xpos1[:-1] / 60.
    print x.shape, num1.shape, np.zeros(x.shape[0]).shape
    trj.plots.fill_between_steps(x, num1, np.zeros(x.shape[0]), ax = plt.gca(),
                                 linewidth = 1, alpha = 1, edgecolor = 'red',
                                 color = 'none', label = 'OCTc')
    trj.plots.fill_between_steps(x, num2, np.zeros(x.shape[0]), ax = plt.gca(),
                                 linewidth = 1, alpha = 1, edgecolor = 'green',
                                 color = 'none', label = 'OCTnc')
    plt.plot([-1000, -1000], [-1000, -1000], color = 'red', label = 'OCTc')
    plt.plot([-1000, -1000], [-1000, -1000], color = 'green', label = 'OCTnc')
    plt.gca().set_xlim(0, 120)
    plt.gca().set_ylim(0, 600)
    plt.legend()
    plt.xlabel('Forecast lead time [h]')
    plt.ylabel('Frequency')
    plt.title('Time of 750 hPa crossing')
    #plt.plot(x, num1, color = 'red', label = 'OCTc'
    plt.tight_layout()
    plt.savefig(plotdir + 'OCT_temporal', dpi = 300, bbox_inches = 'tight')
    plt.close()
    
if 'temporal_JANJUL' in plotlist:
    conv = (c2._mask_array('WCB_Cy1_new', 'Pcross500') * 5. + 
            c2._mask_array('WCB_Cy1_new', 'startt'))
    
    fig = plt.figure(figsize = (95./25.4, 3))
    num1, xpos1 = np.histogram(conv, bins = 120, range = [0, 7200] )
    x = xpos1[:-1] / 60.
    print x.shape, num1.shape, np.zeros(x.shape[0]).shape
    trj.plots.fill_between_steps(x, num1, np.zeros(x.shape[0]), ax = plt.gca(),
                                 linewidth = 1, alpha = 1, edgecolor = 'red',
                                 color = 'none', label = 'JAN')
    
    plt.plot([-1000, -1000], [-1000, -1000], color = 'red', label = 'JAN')
    plt.gca().set_xlim(0, 120)
    plt.gca().set_ylim(0, 600)
    plt.legend()
    plt.xlabel('Forecast lead time [h]')
    plt.ylabel('Frequency')
    plt.title('Time of 500 hPa crossing')
    #plt.plot(x, num1, color = 'red', label = 'OCTc'
    plt.tight_layout()
    plt.savefig(plotdir + 'JAN_temporal', dpi = 300, bbox_inches = 'tight')
    plt.close()
    
    conv = (c1._mask_array('WCB_Cy1', 'Pcross500') * 5. + 
            c1._mask_array('WCB_Cy1', 'startt'))
    
    fig = plt.figure(figsize = (95./25.4, 3))
    num1, xpos1 = np.histogram(conv, bins = 120, range = [0, 7200] )
    x = xpos1[:-1] / 60.
    print x.shape, num1.shape, np.zeros(x.shape[0]).shape
    trj.plots.fill_between_steps(x, num1, np.zeros(x.shape[0]), ax = plt.gca(),
                                 linewidth = 1, alpha = 1, edgecolor = 'red',
                                 color = 'none', label = 'JUL')
    
    plt.plot([-1000, -1000], [-1000, -1000], color = 'red', label = 'JUL')
    plt.gca().set_xlim(0, 120)
    plt.gca().set_ylim(0, 2500)
    plt.legend()
    plt.xlabel('Forecast lead time [h]')
    plt.ylabel('Frequency')
    plt.title('Time of 500 hPa crossing')
    #plt.plot(x, num1, color = 'red', label = 'OCTc'
    plt.tight_layout()
    plt.savefig(plotdir + 'JUL_temporal', dpi = 300, bbox_inches = 'tight')
    plt.close()
    
    conv = (c0._mask_array('Wconv2', 'Pcross500') * 5. + 
            c0._mask_array('Wconv2', 'startt'))
    nonconv = (c0._mask_array('Wnonconv2', 'Pcross500') * 5. + 
               c0._mask_array('Wnonconv2', 'startt'))
    fig = plt.figure(figsize = (95./25.4, 3))
    num1, xpos1 = np.histogram(conv, bins = 120, range = [0, 7200] )
    num2, xpos2 = np.histogram(nonconv, bins = 120, range = [0, 7200] )
    x = xpos1[:-1] / 60.
    print x.shape, num1.shape, np.zeros(x.shape[0]).shape
    trj.plots.fill_between_steps(x, num1, np.zeros(x.shape[0]), ax = plt.gca(),
                                 linewidth = 1, alpha = 1, edgecolor = 'red',
                                 color = 'none', label = 'OCTc')
    trj.plots.fill_between_steps(x, num2, np.zeros(x.shape[0]), ax = plt.gca(),
                                 linewidth = 1, alpha = 1, edgecolor = 'green',
                                 color = 'none', label = 'OCTnc')
    plt.plot([-1000, -1000], [-1000, -1000], color = 'red', label = 'OCTc')
    plt.plot([-1000, -1000], [-1000, -1000], color = 'green', label = 'OCTnc')
    plt.gca().set_xlim(0, 120)
    plt.gca().set_ylim(0, 600)
    plt.legend()
    plt.xlabel('Forecast lead time [h]')
    plt.ylabel('Frequency')
    plt.title('Time of 500 hPa crossing')
    #plt.plot(x, num1, color = 'red', label = 'OCTc'
    plt.tight_layout()
    plt.savefig(plotdir + 'OCT_temporal_500', dpi = 300, bbox_inches = 'tight')
    plt.close()

if 'pv_hist' in plotlist:
    #c0.draw_pv_hist('var4', 'Wconv2', 'P600', 'OCTc', 
                    #savebase = plotdir + 'OCTc_4_', ex = 4)
    #c0.draw_pv_hist('var4', 'Wnonconv2', 'P600', 'OCTnc', 
                    #savebase = plotdir + 'OCTnc_4_', ex = 4)
    #c1.draw_pv_hist('POT_VORTIC', 'WCB_Cy1', 'P600', 'JUL', 
                    #savebase = plotdir + 'JUL_4_', ex = 4)
    #c2.draw_pv_hist('POT_VORTIC', 'WCB_Cy1_new', 'P600', 'JAN', 
                    #savebase = plotdir + 'JAN_4_', ex = 4)
    
    #c0.draw_pv_hist('var4', 'Wconv2', 'P600', 'OCTc', 
                    #savebase = plotdir + 'OCTc_1_', ex = 1)
    #c0.draw_pv_hist('var4', 'Wnonconv2', 'P600', 'OCTnc', 
                    #savebase = plotdir + 'OCTnc_1_', ex = 1)
    #c1.draw_pv_hist('POT_VORTIC', 'WCB_Cy1', 'P600', 'JUL', 
                    #savebase = plotdir + 'JUL_1_', ex = 1)
    c2.draw_pv_hist('POT_VORTIC', 'WCB_Cy1_new', 'P600', 'JAN', 
                    savebase = plotdir + 'JAN_1_', ex = 1)
    
if 'pv_pdf' in plotlist:
    #octc = c0.draw_pv_pdf('var4', 'Wconv2', 'P600', 'OCTc', 
                    #savebase = plotdir + 'OCTc_4_', ex = 4, mx = 20, ret = True)
    
    #jan = c2.draw_pv_pdf('POT_VORTIC', 'WCB_Cy1_new', 'P600', 'JAN', 
                    #savebase = plotdir + 'JAN_4_', ex = 4, mx = 20, ret = True)
    
    #octnc = c0.draw_pv_pdf('var4', 'Wnonconv2', 'P600', 'OCTnc', 
                    #savebase = plotdir + 'OCTnc_4_', ex = 4, mx = 20, ret = True)
    #jul = c1.draw_pv_pdf('POT_VORTIC', 'WCB_Cy1', 'P600', 'JUL', 
                    #savebase = plotdir + 'JUL_4_', ex = 4, mx = 20, ret = True)
    
    #f = open('./pv_pdf.cpkl', 'w')
    #cPickle.dump(octnc, f)
    #cPickle.dump(jan, f)
    #cPickle.dump(octc, f)
    #cPickle.dump(jul, f)
    #f.close()
    
    # Load data
    f = open('./pv_pdf.cpkl', 'r')
    octnc = cPickle.load(f)
    jan = cPickle.load(f)
    octc = cPickle.load(f)
    jul = cPickle.load(f)
    f.close()
    
    #### return [num1, num2, x1, mean1, mean2]
    # Plot data
    
    fig, axarr = plt.subplots(2, 2, figsize = (95./25.4*2, 4.2), 
                              sharex = True, sharey = True)
    zero = np.zeros(jul[0].shape[0])
    mx = 25
    
    # Top left
    trj.plots.fill_between_steps(jul[2][:-1], jul[0], zero, ax = axarr[0][0], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'steelblue', label = '3h pre-ascent')
    trj.plots.fill_between_steps(jul[2][:-1], jul[1], zero, ax = axarr[0][0], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'tomato', label = '3h post-ascent')
    axarr[0][0].plot([jul[3], jul[3]], [0, mx], color = 'steelblue', label = 'pre-ascent mean')
    axarr[0][0].plot([jul[4], jul[4]], [0, mx], color = 'tomato', label = 'post-ascent mean')
    plt.text(0.68, 0.78, 'Pre-ascent \nmean: {:.2f}pvu'.format(jul[3]), 
             transform = axarr[0][0].transAxes, 
             fontsize = 9, color = 'steelblue', ha = 'left')
    plt.text(0.68, 0.59, 'Post-ascent \nmean: {:.2f}pvu'.format(jul[4]), 
             transform = axarr[0][0].transAxes, 
             fontsize = 9, color = 'tomato', ha = 'left')
    plt.text(0.038, 0.85, 'JUL', 
             transform = axarr[0][0].transAxes, 
             fontsize = 10)
    
    # Top right
    trj.plots.fill_between_steps(jan[2][:-1], jan[0], zero, ax = axarr[0][1], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'steelblue', label = '3h pre-ascent')
    trj.plots.fill_between_steps(jan[2][:-1], jan[1], zero, ax = axarr[0][1], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'tomato', label = '3h post-ascent')
    axarr[0][1].plot([jan[3], jan[3]], [0, mx], color = 'steelblue', label = 'pre-ascent mean')
    axarr[0][1].plot([jan[4], jan[4]], [0, mx], color = 'tomato', label = 'post-ascent mean')
    plt.text(0.68, 0.78, 'Pre-ascent \nmean: {:.2f}pvu'.format(jan[3]), 
             transform = axarr[0][1].transAxes, 
             fontsize = 9, color = 'steelblue', ha = 'left')
    plt.text(0.68, 0.59, 'Post-ascent \nmean: {:.2f}pvu'.format(jan[4]), 
             transform = axarr[0][1].transAxes, 
             fontsize = 9, color = 'tomato', ha = 'left')
    plt.text(0.038, 0.85, 'JAN', 
             transform = axarr[0][1].transAxes, 
             fontsize = 10)
    
    # Bottom left
    trj.plots.fill_between_steps(octc[2][:-1], octc[0], zero, ax = axarr[1][0], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'steelblue', label = '3h pre-ascent')
    trj.plots.fill_between_steps(octc[2][:-1], octc[1], zero, ax = axarr[1][0], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'tomato', label = '3h post-ascent')
    axarr[1][0].plot([octc[3], octc[3]], [0, mx], color = 'steelblue', label = 'pre-ascent mean')
    axarr[1][0].plot([octc[4], octc[4]], [0, mx], color = 'tomato', label = 'post-ascent mean')
    plt.text(0.68, 0.78, 'Pre-ascent \nmean: {:.2f}pvu'.format(octc[3]), 
             transform = axarr[1][0].transAxes, 
             fontsize = 9, color = 'steelblue', ha = 'left')
    plt.text(0.68, 0.59, 'Post-ascent \nmean: {:.2f}pvu'.format(octc[4]), 
             transform = axarr[1][0].transAxes, 
             fontsize = 9, color = 'tomato', ha = 'left')
    plt.text(0.038, 0.85, 'OCTc', 
             transform = axarr[1][0].transAxes, 
             fontsize = 10)
    
    # Bottom right
    trj.plots.fill_between_steps(octnc[2][:-1], octnc[0], zero, ax = axarr[1][1], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'steelblue', label = '3h pre-ascent')
    trj.plots.fill_between_steps(octnc[2][:-1], octnc[1], zero, ax = axarr[1][1], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'tomato', label = '3h post-ascent')
    axarr[1][1].plot([octnc[3], octnc[3]], [0, mx], color = 'steelblue', label = 'pre-ascent mean')
    axarr[1][1].plot([octnc[4], octnc[4]], [0, mx], color = 'tomato', label = 'post-ascent mean')
    plt.text(0.68, 0.78, 'Pre-ascent \nmean: {:.2f}pvu'.format(octnc[3]), 
             transform = axarr[1][1].transAxes, 
             fontsize = 9, color = 'steelblue', ha = 'left')
    plt.text(0.68, 0.59, 'Post-ascent \nmean: {:.2f}pvu'.format(octnc[4]), 
             transform = axarr[1][1].transAxes, 
             fontsize = 9, color = 'tomato', ha = 'left')
    plt.text(0.038, 0.85, 'OCTnc', 
             transform = axarr[1][1].transAxes, 
             fontsize = 10)
    
    
    
    # Adjustments
    axarr[0, 0].set_ylabel('Frequency [%]')
    axarr[1, 0].set_ylabel('Frequency [%]')
    
    axarr[1, 0].set_yticks(range(0, 20, 5))
    
    axarr[1, 0].set_xlabel('PV [pvu]')
    axarr[1, 1].set_xlabel('PV [pvu]')
    
    axarr[1, 0].set_xticks(range(-3, 4, 1))
    
    axarr[0, 0].set_ylim(0, 20)
    axarr[0, 0].set_xlim(-4, 4)
    
    plt.suptitle('(a) Potential vorticity', fontsize = 10)
    plt.subplots_adjust(left = 0.07, right = 0.99, top = 0.925, bottom = 0.12, 
                        wspace = 0.03, hspace = 0.05)
    fig.savefig(plotdir + 'PV_pdf', dpi = 300)

if 'theta_pdf' in plotlist:
    #octc = c0.draw_pv_pdf('THETA', 'Wconv2', 'P600', 'OCTc', 
                    #savebase = plotdir + 'OCTc_4_', ex = 4, mx = 25, var = 'THETA',
                    #ret = True)
    
    
    #octnc = c0.draw_pv_pdf('THETA', 'Wnonconv2', 'P600', 'OCTnc', 
                    #savebase = plotdir + 'OCTnc_4_', ex = 4, mx = 25, var = 'THETA',
                    #ret = True)
    #jul = c1.draw_pv_pdf('THETA', 'WCB_Cy1', 'P600', 'JUL', 
                    #savebase = plotdir + 'JUL_4_', ex = 4, mx = 25, var = 'THETA',
                    #ret = True)
    #jan = c2.draw_pv_pdf('THETA', 'WCB_Cy1_new', 'P600', 'JAN', 
                    #savebase = plotdir + 'JAN_4_', ex = 4, mx = 25, var = 'THETA',
                    #ret = True)
    #f = open('./theta_pdf.cpkl', 'w')
    #cPickle.dump(octnc, f)
    #cPickle.dump(jan, f)
    #cPickle.dump(octc, f)
    #cPickle.dump(jul, f)
    #f.close()
    
    
        # Load data
    f = open('./theta_pdf.cpkl', 'r')
    octnc = cPickle.load(f)
    jan = cPickle.load(f)
    octc = cPickle.load(f)
    jul = cPickle.load(f)
    f.close()

    
    #### return [num1, num2, x1, mean1, mean2]
    # Plot data
    
    fig, axarr = plt.subplots(2, 2, figsize = (95./25.4*2, 4.2), 
                              sharex = True, sharey = True)
    zero = np.zeros(jul[0].shape[0])
    mx = 25
    
    # Top left
    trj.plots.fill_between_steps(jul[2][:-1], jul[0], zero, ax = axarr[0][0], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'steelblue', label = '3h pre-ascent')
    trj.plots.fill_between_steps(jul[2][:-1], jul[1], zero, ax = axarr[0][0], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'tomato', label = '3h post-ascent')
    axarr[0][0].plot([jul[3], jul[3]], [0, mx], color = 'steelblue', label = 'pre-ascent mean')
    axarr[0][0].plot([jul[4], jul[4]], [0, mx], color = 'tomato', label = 'post-ascent mean')
    plt.text(0.73, 0.78, 'Pre-ascent\nmean: {:.0f}K'.format(jul[3]), 
             transform = axarr[0][0].transAxes, 
             fontsize = 9, color = 'steelblue', ha = 'left')
    plt.text(0.73, 0.59, 'Post-ascent\nmean: {:.0f}K'.format(jul[4]), 
             transform = axarr[0][0].transAxes, 
             fontsize = 9, color = 'tomato', ha = 'left')
    plt.text(0.038, 0.85, 'JUL', 
             transform = axarr[0][0].transAxes, 
             fontsize = 10)
    
    # Top right
    trj.plots.fill_between_steps(jan[2][:-1], jan[0], zero, ax = axarr[0][1], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'steelblue', label = '3h pre-ascent')
    trj.plots.fill_between_steps(jan[2][:-1], jan[1], zero, ax = axarr[0][1], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'tomato', label = '3h post-ascent')
    axarr[0][1].plot([jan[3], jan[3]], [0, mx], color = 'steelblue', label = 'pre-ascent mean')
    axarr[0][1].plot([jan[4], jan[4]], [0, mx], color = 'tomato', label = 'post-ascent mean')
    plt.text(0.73, 0.78, 'Pre-ascent\nmean: {:.0f}K'.format(jan[3]), 
             transform = axarr[0][1].transAxes, 
             fontsize = 9, color = 'steelblue', ha = 'left')
    plt.text(0.73, 0.59, 'Post-ascent\nmean: {:.0f}K'.format(jan[4]), 
             transform = axarr[0][1].transAxes, 
             fontsize = 9, color = 'tomato', ha = 'left')
    plt.text(0.038, 0.85, 'JAN', 
             transform = axarr[0][1].transAxes, 
             fontsize = 10)
    
    # Bottom left
    trj.plots.fill_between_steps(octc[2][:-1], octc[0], zero, ax = axarr[1][0], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'steelblue', label = '3h pre-ascent')
    trj.plots.fill_between_steps(octc[2][:-1], octc[1], zero, ax = axarr[1][0], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'tomato', label = '3h post-ascent')
    axarr[1][0].plot([octc[3], octc[3]], [0, mx], color = 'steelblue', label = 'pre-ascent mean')
    axarr[1][0].plot([octc[4], octc[4]], [0, mx], color = 'tomato', label = 'post-ascent mean')
    plt.text(0.73, 0.78, 'Pre-ascent\nmean: {:.0f}K'.format(octc[3]), 
             transform = axarr[1][0].transAxes, 
             fontsize = 9, color = 'steelblue', ha = 'left')
    plt.text(0.73, 0.59, 'Post-ascent\nmean: {:.0f}K'.format(octc[4]), 
             transform = axarr[1][0].transAxes, 
             fontsize = 9, color = 'tomato', ha = 'left')
    plt.text(0.038, 0.85, 'OCTc', 
             transform = axarr[1][0].transAxes, 
             fontsize = 10)
    
    # Bottom right
    trj.plots.fill_between_steps(octnc[2][:-1], octnc[0], zero, ax = axarr[1][1], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'steelblue', label = '3h pre-ascent')
    trj.plots.fill_between_steps(octnc[2][:-1], octnc[1], zero, ax = axarr[1][1], 
                                 linewidth = 1., alpha = 0.7, 
                                 color = 'tomato', label = '3h post-ascent')
    axarr[1][1].plot([octnc[3], octnc[3]], [0, mx], color = 'steelblue', label = 'pre-ascent mean')
    axarr[1][1].plot([octnc[4], octnc[4]], [0, mx], color = 'tomato', label = 'post-ascent mean')
    plt.text(0.73, 0.78, 'Pre-ascent\nmean: {:.0f}K'.format(octnc[3]), 
             transform = axarr[1][1].transAxes, 
             fontsize = 9, color = 'steelblue', ha = 'left')
    plt.text(0.73, 0.59, 'Post-ascent\nmean: {:.0f}K'.format(octnc[4]), 
             transform = axarr[1][1].transAxes, 
             fontsize = 9, color = 'tomato', ha = 'left')
    plt.text(0.038, 0.85, 'OCTnc', 
             transform = axarr[1][1].transAxes, 
             fontsize = 10)
    
    
    
    # Adjustments
    axarr[0, 0].set_ylabel('Frequency [%]')
    axarr[1, 0].set_ylabel('Frequency [%]')
    
    axarr[1, 0].set_yticks([0, 5, 10, 15, 20])
    
    axarr[1, 0].set_xlabel('Theta [K]')
    axarr[1, 1].set_xlabel('Theta [K]')
    
    axarr[1, 1].set_xticks(range(290, 350, 10))
    
    axarr[0, 0].set_ylim(0, 25)
    axarr[0, 0].set_xlim(280, 350)
    
    plt.suptitle('(b) Potential temperature', fontsize = 10)
    plt.subplots_adjust(left = 0.07, right = 0.99, top = 0.925, bottom = 0.12, 
                        wspace = 0.03, hspace = 0.05)
    fig.savefig(plotdir + 'THETA_pdf', dpi = 300)
    
    
    
if 'theta_hist' in plotlist:
    c1.draw_dtheta_hist([2280, 2340, 2400, 2460, 2520], lllimit = (100, 200), urlimit = (800, 746), 
                      savebase = plotdir + 'JUL_',
                      level = 500, title = 'JUL')
    c2.draw_dtheta_hist([1200, 1260, 1320, 1380, 1440], lllimit = (0, 0), urlimit = (1250, 885), 
                      savebase = plotdir + 'JAN_',
                      level = 500, title = 'JAN')
    #c0.draw_dtheta_hist([2640, 2700, 2760, 2820, 2880], lllimit = (400, 250), urlimit = (1400, 900), 
                      #savebase = plotdir + 'OCTnc_',
                      #level = 500, title = 'OCTnc')
    #c0.draw_dtheta_hist(3600, lllimit = (400, 250), urlimit = (1400, 900), 
                      #savebase = plotdir + 'OCTc_',
                      #level = 500, title = 'OCTc')
if 'vort_hist' in plotlist:
    c1.draw_vort_hist([2280, 2340, 2400, 2460, 2520], lllimit = (100, 200), urlimit = (800, 746), 
                        savebase = plotdir + 'JUL_',
                        level = 500, title = 'JUL')
    c2.draw_vort_hist([1200, 1260, 1320, 1380, 1440], lllimit = (0, 0), urlimit = (1250, 885), 
                        savebase = plotdir + 'JAN_',
                        level = 500, title = 'JAN')
    #c0.draw_vort_hist(2760, lllimit = (400, 250), urlimit = (1400, 900), 
                      #savebase = plotdir + 'OCTnc_',
                      #level = 500, title = 'OCTnc')
    #c0.draw_vort_hist(3600, lllimit = (400, 250), urlimit = (1400, 900), 
                      #savebase = plotdir + 'OCTc_',
                      #level = 500, title = 'OCTc')
    
if 'pv_vs_t' in plotlist:
    
    #octnc= c0.draw_centered_vs_t(['var4'], ['Wnonconv2'], 'Pcross600', plottype = 'Paper', ylim = (-6, 6), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/tmp/paper_test/PV_vs_t_nc_', legnames = ['OCTnc', 'JAN'], letter = '', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3, ret = True)
    #print octnc
    #jan = c2.draw_centered_vs_t(['POT_VORTIC'], ['WCB_Cy1_new'], 'Pcross600', plottype = 'Paper', ylim = (-6, 6), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/tmp/paper_test/PV_vs_t_nc_', legnames = ['OCTnc', 'JAN'], letter = '', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3, ret = True)
    #print jan
    #octc= c0.draw_centered_vs_t(['var4'], ['Wconv2'], 'Pcross600', plottype = 'Paper', ylim = (-6, 6), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/tmp/paper_test/PV_vs_t_c_', legnames = ['OCTc', 'JUL'], letter = '', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3 , ret = True)
    
    #jul = c1.draw_centered_vs_t(['POT_VORTIC'], ['WCB_Cy1'], 'Pcross600', plottype = 'Paper', ylim = (-6, 6), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/tmp/paper_test/PV_vs_t_c_', legnames = ['OCTc', 'JUL'], letter = '', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3 , ret = True)

    #f = open('./mean_vs_t.cpkl', 'w')
    #cPickle.dump(octnc, f)
    #cPickle.dump(jan, f)
    #cPickle.dump(octc, f)
    #cPickle.dump(jul, f)
    #f.close()
    
    
    # Load data
    f = open('./mean_vs_t.cpkl', 'r')
    octnc = cPickle.load(f)
    jan = cPickle.load(f)
    octc = cPickle.load(f)
    jul = cPickle.load(f)
    f.close()
    
    # Plot data 
    fig, ax = plt.subplots(1, 1, figsize = (95./25.4, 3))
    
    ax.plot(jul[1], jul[0], color = 'orange', label = 'JUL', linewidth = 1.5)
    ax.plot(jan[1], jan[0], color = 'green', label = 'JAN', linewidth = 1.5)
    ax.plot(octc[1], octc[0], color = 'magenta', label = 'OCTc', linewidth = 1.5)
    ax.plot(octnc[1], octnc[0], color = 'indigo', label = 'OCTnc', linewidth = 1.5)
    
    ax.set_xlim(-12, 12)
    ax.set_xlabel('Time relative to 600 hPa crossing [h]')
    ax.set_ylabel('PV [pvu]')
    ax.set_ylim(-1, 5)
    #ax.grid()
    plt.legend(loc = 2)
    
    
    plt.text(0.9, 0.9, '(a)', 
             transform = ax.transAxes, 
             fontsize = 10)
    
    ax.plot([-12, 12], [0,0], linewidth = 0.5, c = 'gray', linestyle = '--')
    ax.plot([0, 0], [-1,5], linewidth = 0.5, c = 'gray', linestyle = '--', 
            zorder = 0.001)
    
    ax.set_xticks(range(-12, 15, 3))
    
    plt.subplots_adjust(left = 0.16, right = 0.97, top = 0.95, bottom = 0.15)
    fig.savefig(plotdir + 'PV_mean', dpi = 300)

if 'theta_vs_t' in plotlist:
    
    #octnc= c0.draw_centered_vs_t(['THETA'], ['Wnonconv2'], 'Pcross600', plottype = 'Paper', ylim = (-6, 6), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/tmp/paper_test/PV_vs_t_nc_', legnames = ['OCTnc', 'JAN'], letter = '', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3, ret = True)
    
    #jan = c2.draw_centered_vs_t(['THETA'], ['WCB_Cy1_new'], 'Pcross600', plottype = 'Paper', ylim = (-6, 6), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/tmp/paper_test/PV_vs_t_nc_', legnames = ['OCTnc', 'JAN'], letter = '', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3, ret = True)
    
    #octc= c0.draw_centered_vs_t(['THETA'], ['Wconv2'], 'Pcross600', plottype = 'Paper', ylim = (-6, 6), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/tmp/paper_test/PV_vs_t_c_', legnames = ['OCTc', 'JUL'], letter = '', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3 , ret = True)
    
    #jul = c1.draw_centered_vs_t(['THETA'], ['WCB_Cy1'], 'Pcross600', plottype = 'Paper', ylim = (-6, 6), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/tmp/paper_test/PV_vs_t_c_', legnames = ['OCTc', 'JUL'], letter = '', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3 , ret = True)

    #f = open('./theta_mean_vs_t.cpkl', 'w')
    #cPickle.dump(octnc, f)
    #cPickle.dump(jan, f)
    #cPickle.dump(octc, f)
    #cPickle.dump(jul, f)
    #f.close()
    
    
    # Load data
    f = open('./theta_mean_vs_t.cpkl', 'r')
    octnc = cPickle.load(f)
    jan = cPickle.load(f)
    octc = cPickle.load(f)
    jul = cPickle.load(f)
    f.close()
    
    # Plot data 
    fig, ax = plt.subplots(1, 1, figsize = (95./25.4, 3))
    
    ax.plot(jul[1], jul[0], color = 'orange', label = 'JUL', linewidth = 1.5)
    ax.plot(jan[1], jan[0], color = 'green', label = 'JAN', linewidth = 1.5)
    ax.plot(octc[1], octc[0], color = 'magenta', label = 'OCTc', linewidth = 1.5)
    ax.plot(octnc[1], octnc[0], color = 'indigo', label = 'OCTnc', linewidth = 1.5)
    
    ax.set_xlim(-12, 12)
    ax.set_xlabel('Time relative to 600 hPa crossing [h]')
    ax.set_ylabel('Theta [K]')
    ax.set_ylim(280, 340)
    #ax.grid()
    #plt.legend(loc = 2)
    
    
    plt.text(0.9, 0.9, '(b)', 
             transform = ax.transAxes, 
             fontsize = 10)
    
    ax.plot([-12, 12], [0,0], linewidth = 0.5, c = 'gray', linestyle = '--')
    ax.plot([0, 0], [280, 340], linewidth = 0.5, c = 'gray', linestyle = '--', 
            zorder = 0.001)
    
    ax.set_xticks(range(-12, 15, 3))
    
    plt.subplots_adjust(left = 0.16, right = 0.97, top = 0.95, bottom = 0.15)
    fig.savefig(plotdir + 'THETA_mean', dpi = 300)

    
    
    
    
