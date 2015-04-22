import traj_tools as trj

c0 = trj.loadme('c0_rot.trj')
#c0.create_filter('WCB_Conv', [('P600', 0, 2880), ('P400', 0, 120)])
c1 = trj.loadme('c1.trj')


## 061114
#c0.draw_centered_vs_t('P', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('P', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')


#c0.draw_centered_vs_t('var4', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (-3, 6), sigma = 1, idtext = '061114c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('var4', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (-3, 6), sigma = 1, idtext = '061114d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

#c1.draw_centered_vs_t('P', 'WCB', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
#c1.draw_centered_vs_t('var4', 'WCB', 'P600', plottype='Smooth', ylim = (-3, 6), sigma = 1, idtext = '061114f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')


## 111114
#c0.draw_centered_vs_t('THETA', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (280, 330), sigma = 1, idtext = '111114a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('THETA', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (280, 330), sigma = 1, idtext = '111114b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('QV', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (0, 0.02), sigma = 1, idtext = '111114c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('QV', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (0, 0.02), sigma = 1, idtext = '111114d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('THETAE', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (280, 330), sigma = 1, idtext = '111114e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('THETAE', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (280, 330), sigma = 1, idtext = '111114f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

#c1.draw_centered_vs_t('THETA', 'WCB', 'P600', plottype='Smooth', ylim = (280, 330), sigma = 1, idtext = '111114g', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
#c1.draw_centered_vs_t('QV', 'WCB', 'P600', plottype='Smooth', ylim = (0, 0.02), sigma = 1, idtext = '111114h', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')



# 121114
#c0.draw_centered_vs_t('Q1', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (0, 0.0002), 
#                      sigma = 1, idtext = '121114a', 
#                      savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('Q2', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (0, 0.0001), 
#                      sigma = 1, idtext = '121114b', 
#                      savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('Q3', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (0, 0.0003), 
#                      sigma = 1, idtext = '121114c', 
#                      savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('Q1', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (0, 0.0002), 
#                      sigma = 1, idtext = '12111d', 
#                      savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('Q2', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (0, 0.0001), 
#                     sigma = 1, idtext = '121114e', 
#                      savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('Q3', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (0, 0.0003), 
#                      sigma = 1, idtext = '121114f', 
#                      savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c1.draw_centered_vs_t('QR', 'WCB', 'P600', plottype='Smooth', ylim = (0, 0.0002), 
                      #sigma = 1, idtext = '12111g', 
                      #savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
#c1.draw_centered_vs_t('QC', 'WCB', 'P600', plottype='Smooth', ylim = (0, 0.0001), 
                      #sigma = 1, idtext = '121114h', 
                      #savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')

#c0.new_max_diff('z')
#c0.draw_hist('z_max_diff', savebase = 
             #'/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             #idtext = '121114i', range = (0, 15), filtername = 'WCB_Conv')
#c0.draw_hist('z_max_diff', savebase = 
             #'/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             #idtext = '121114j', range = (0, 15), filtername = 'WCB_NonConv')

c1.calc_theta()
c1.draw_centered_vs_t('THETAE', 'WCB', 'P600', plottype='Smooth', ylim = (290, 360), sigma = 1, idtext = '111114i', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/') 

# 131114
#with plt.style.context('fivethirtyeight'):
    #c0.draw_hist('P600', mintohrs = True, idtext = '131114a', filtername = 'WCB')
    #plt.gca().set_xlim(0, 48)
    #plt.xlabel('Ascent time for 600hPa [hrs]')
    #plt.xticks(range(0, 48 + 6, 6))
    #savename = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist_P600_WCB_131114a'
    #plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

#plt.close('all')
#with plt.style.context('fivethirtyeight'):
    #c0.draw_hist('P600', mintohrs = True, idtext = '131114b', filtername = 'WCB_Conv')
    #plt.gca().set_xlim(0, 48)
    #plt.xlabel('Ascent time for 600hPa [hrs]')
    #plt.xticks(range(0, 48 + 6, 6))
    #savename = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist_P600_WCB_Conv_131114b'
    #plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

#plt.close('all')
#with plt.style.context('fivethirtyeight'):
    #c0.draw_hist('P600', mintohrs = True, idtext = '131114c', filtername = 'WCB_NonConv')
    #plt.gca().set_xlim(0, 48)
    #plt.xlabel('Ascent time for 600hPa [hrs]')
    #plt.xticks(range(0, 48 + 6, 6))
    #savename = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist_P600_WCB_NonConv_131114c'
    #plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

#plt.close('all')
#with plt.style.context('fivethirtyeight'):
    #c1.draw_hist('P600', mintohrs = True, idtext = '131114d', filtername = 'WCB')
    #plt.gca().set_xlim(0, 48)
    #plt.xlabel('Ascent time for 600hPa [hrs]')
    #plt.xticks(range(0, 48 + 6, 6))
    #savename = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/hist_P600_WCB_131114d'
    #plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

# 141114

for t in range(360, c1.maxmins + 60, 60):
    c1.draw_trj_evo([('var4', 15)], filtername='WCB_2.5PV_360', tafter = t, 
                    savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
    
    
# 121214

c0.draw_vs_p('THETA', 'WCB_NonConv', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '121214a', ylim = (280, 350))
c0.draw_vs_p('THETA', 'WCB_Conv', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '121214b', ylim = (280, 350))
c1.draw_vs_p('THETA', 'WCB', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', idtext = '121214c', ylim = (280, 350))
c0.draw_vs_p('THETAE', 'WCB_NonConv', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '121214d', ylim = (300, 350))
c0.draw_vs_p('THETAE', 'WCB_Conv', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '121214e', ylim = (300, 350))
c1.draw_vs_p('THETAE', 'WCB', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', idtext = '121214f', ylim = (300, 350))  
c0.draw_vs_p('var4', 'WCB_NonConv', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '121214g', ylim = (-5, 10))
c0.draw_vs_p('var4', 'WCB_Conv', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '121214h', ylim = (-5, 10))
c1.draw_vs_p('POT_VORTIC', 'WCB', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', idtext = '121214i', ylim = (-5, 10))  

#171214
c0.draw_vs_p('THETA', 'WCB_NonConv', 'P600', (100, 1000), idtext = '171214a', ylim = (280, 330), ylabel = 'THETA [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_vs_p('THETA', 'WCB_Conv', 'P600', (100, 1000), idtext = '171214b', ylim = (285, 340), ylabel = 'THETA [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c1.draw_vs_p('THETA', 'WCB', 'P600', (50, 1000), idtext = '171214e', ylim = (280, 350), ylabel = 'THETA [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c0.draw_vs_p('THETAE_new', 'WCB_NonConv', 'P600', (100, 1000), idtext = '171214c', ylim = (300, 330), ylabel = 'THETA_E [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_vs_p('THETAE_new', 'WCB_Conv', 'P600', (100, 1000), idtext = '171214d', ylim = (300, 340), ylabel = 'THETA_E [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c1.draw_vs_p('THETAE_new', 'WCB', 'P600', (50, 1000), idtext = '171214f', ylim = (300, 350), ylabel = 'THETA_E [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c1.draw_vs_p('POT_VORTIC', 'WCB', 'P600', (50, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', idtext = '171214g', ylim = (-5, 10), ylabel = 'PV [PVU]')  
c0.draw_vs_p('QV', 'WCB_NonConv', 'P600', (100, 1000), idtext = '171214h', ylabel = 'QV [kg/kg]', ylim = (0.00001, 0.1), log = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_vs_p('QV', 'WCB_Conv', 'P600', (100, 1000), idtext = '171214i', ylabel = 'QV [kg/kg]', ylim = (0.00001, 0.1), log = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c1.draw_vs_p('QV', 'WCB', 'P600', (50, 1000), idtext = '171214j', ylabel = 'QV [kg/kg]', ylim = (0.00001, 0.1), log = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')

# 200115
with plt.style.context('fivethirtyeight'):
    c2.draw_hist('P600', mintohrs = True, idtext = '200115a', filtername = 'WCB')
    plt.gca().set_xlim(0, 48)
    plt.xlabel('Ascent time for 600hPa [hrs]')
    plt.xticks(range(0, 48 + 6, 6))
    savename = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/hist_P600_WCB_200115a'
    plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
# P
c0.draw_centered_vs_t('P', 'WCB_NonConv', 'P600', plottype='Fill', ylim = (200, 1000), sigma = 1, idtext = '200115b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('P', 'WCB_Conv', 'P600', plottype='Fill', ylim = (200, 1000), sigma = 1, idtext = '200115c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c1.draw_centered_vs_t('P', 'WCB', 'P600', plottype='Fill', ylim = (200, 1000), sigma = 1, idtext = '200115d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c2.draw_centered_vs_t('P', 'WCB', 'P600', plottype='Fill', ylim = (200, 1000), sigma = 1, idtext = '200115e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
#PV
c0.draw_centered_vs_t('var4', 'WCB_NonConv', 'P600', plottype='Fill', ylim = (-3, 6), sigma = 1, idtext = '200115f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('var4', 'WCB_Conv', 'P600', plottype='Fill', ylim = (-3, 6), sigma = 1, idtext = '20011g', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c1.draw_centered_vs_t('var4', 'WCB', 'P600', plottype='Fill', ylim = (-3, 6), sigma = 1, idtext = '200115h', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c2.draw_centered_vs_t('POT_VORTIC', 'WCB', 'P600', plottype='Fill', ylim = (-3, 6), sigma = 1, idtext = '200115i', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

# centered_THETA_P600_WCB_200115j.png
# centered_THETAE_P600_WCB_200115k.png

#xy_path
c2.draw_trj_all([], 'WCB', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/', idtext = '200115l', thinning = 10)

#210115
c1.draw_vs_p('POT_VORTIC', 'WCB', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', idtext = '210115a', ylim = (-5, 10), ylabel = 'PVU') 
c2.draw_vs_p('POT_VORTIC', 'WCB', 'P600', (100, 1000), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/', idtext = '210115b', ylim = (-5, 10), ylabel = 'PVU')

c1.draw_centered_vs_t('POT_VORTIC', 'WCB', 'P600', plottype='Fill', ylim = (-3, 6), sigma = 1, idtext = '210115c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')

c2.draw_centered_vs_t('P', 'WCB', 'P600', plottype='Fill', ylim = (200, 1000), xlim = (-24, 24), sigma = 1, idtext = '210115d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/xlim48')
c1.draw_centered_vs_t('P', 'WCB', 'P600', plottype='Fill', ylim = (200, 1000), xlim = (-24, 24), sigma = 1, idtext = '210115e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/xlim48')

c1.draw_centered_vs_t('POT_VORTIC', 'WCB', 'P600', plottype='Fill', ylim = (-3, 6), xlim = (-24, 24), sigma = 1, idtext = '210115f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/xlim48')

#220115
c2.new_cross_level('P', 'P600', 600)
c1.draw_centered_vs_t('P', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (200, 1000), xlim = (-36, 36), select=True, idtext = '220115a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c2.draw_centered_vs_t('P', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (200, 1000), xlim = (-36, 36), select=True, idtext = '220115b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c1.draw_centered_vs_t('POT_VORTIC', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (-3, 6), xlim = (-36, 36), select=True, idtext = '220115f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c2.draw_centered_vs_t('POT_VORTIC', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (-3, 6), xlim = (-36, 36), select=True, idtext = '220115c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c0.draw_centered_vs_t('P', 'WCB_Conv', 'Pcross600', plottype = 'Fill', ylim = (200, 1000), xlim = (-36, 36), select=True, idtext = '220115d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('P', 'WCB_NonConv', 'Pcross600', plottype = 'Fill', ylim = (200, 1000), xlim = (-36, 36), select=True, idtext = '220115e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

c0.draw_centered_vs_t('var4', 'WCB_Conv', 'Pcross600', plottype = 'Fill', ylim = (-3, 6), xlim = (-36, 36), select=True, idtext = '220115g', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('var4', 'WCB_NonConv', 'Pcross600', plottype = 'Fill', ylim = (-3, 6), xlim = (-36, 36), select=True, idtext = '220115h', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

c2.new_loc_filter(-10, 20, 10, 25)
c2.create_filter('WCB_Cy1_new', [('P600', 0, 2880), ('-10201025', True)])
c2.draw_trj_all([], 'WCB_Cy1_new', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/', idtext = '220115h', thinning = 10)

with plt.style.context('fivethirtyeight'):
    c2.draw_hist('P600', mintohrs = True, idtext = '220115i', filtername = 'WCB_Cy1_new')
    plt.gca().set_xlim(0, 48)
    plt.xlabel('Ascent time for 600hPa [hrs]')
    plt.xticks(range(0, 48 + 6, 6))
    savename = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/hist_P600_WCB_Cy1_220115i'
    plt.savefig(savename, dpi = 400, bbox_inches = 'tight')
    
c2.draw_centered_vs_t('P', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (200, 1000), xlim = (-36, 36), select=True, idtext = '220115j', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_centered_vs_t('POT_VORTIC', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (-3, 6), xlim = (-36, 36), select=True, idtext = '220115k', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c2.draw_centered_vs_t('THETA', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=True, idtext = '220115l', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_centered_vs_t('THETAE', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=True, idtext = '220115m', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c0.draw_centered_vs_t('THETA', 'WCB_Conv', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=True, idtext = '220115n', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('THETA', 'WCB_NonConv', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=True, idtext = '220115o', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

c1.draw_centered_vs_t('THETA', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=True, idtext = '220115p', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c1.draw_centered_vs_t('THETAE', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=True, idtext = '220115q', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')

c0.draw_centered_vs_t('THETAE', 'WCB_Conv', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=True, idtext = '220115r', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('THETAE', 'WCB_NonConv', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=True, idtext = '220115s', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

#230115

plt.style.use('fivethirtyeight')

c0.draw_hist('P400', mintohrs = True, idtext = '230115a', filtername = 'WCB')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 400hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist_P400_WCB_230115a'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c0.draw_hist('P500', mintohrs = True, idtext = '230115b', filtername = 'WCB')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 500hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist_P500_WCB_230115b'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c0.draw_hist('P600', mintohrs = True, idtext = '230115d', filtername = 'WCB')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 600hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist_P600_WCB_230115d'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c0.draw_hist('P700', mintohrs = True, idtext = '230115c', filtername = 'WCB')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 700hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist_P700_WCB_230115c'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')


c1.draw_hist('P400', mintohrs = True, idtext = '230115e', filtername = 'WCB')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 400hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/hist_P400_WCB_230115e'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c1.draw_hist('P500', mintohrs = True, idtext = '230115f', filtername = 'WCB')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 500hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/hist_P500_WCB_230115f'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c1.draw_hist('P600', mintohrs = True, idtext = '230115g', filtername = 'WCB')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 600hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/hist_P600_WCB_230115g'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c1.draw_hist('P700', mintohrs = True, idtext = '230115h', filtername = 'WCB')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 700hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/hist_P700_WCB_230115h'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')


c2.draw_hist('P400', mintohrs = True, idtext = '230115i', filtername = 'WCB_Cy1_new')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 400hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/hist_P400_WCB_Cy1_new_230115i'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c2.draw_hist('P500', mintohrs = True, idtext = '230115j', filtername = 'WCB_Cy1_new')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 500hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/hist_P500_WCB_Cy1_new_230115j'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c2.draw_hist('P600', mintohrs = True, idtext = '230115k', filtername = 'WCB_Cy1_new')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 600hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/hist_P600_WCB_Cy1_new_230115k'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c2.draw_hist('P700', mintohrs = True, idtext = '230115l', filtername = 'WCB_Cy1_new')
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 700hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/hist_P700_WCB_Cy1_new_230115l'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c0.draw_trj_dot(["PMSL", "TOT_PREC_S"], interval = 20, filtername='WCB', thinning = 10, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/dot/')

c0.draw_trj_all([], 'WCB', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '230115m', thinning = 10)
c0.draw_trj_all([], 'WCB_Conv', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '230115n', thinning = 10)
c0.draw_trj_all([], 'WCB_NonConv', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '230115o', thinning = 10)


#260115
c2.draw_centered_vs_t('P', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (200, 1000), xlim = (-12, 12), select=True, idtext = '260115a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_centered_vs_t('POT_VORTIC', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (-3, 6), xlim = (-12, 12), select=True, idtext = '260115b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c1.draw_centered_vs_t('QV', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (0, 0.02), xlim = (-36, 36), select=True, idtext = '260115c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')

c0.draw_centered_vs_t('QV', 'WCB_Conv', 'Pcross600', plottype = 'Fill', ylim = (0, 0.02), xlim = (-36, 36), select=True, idtext = '260115d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('QV', 'WCB_NonConv', 'Pcross600', plottype = 'Fill', ylim = (0, 0.02), xlim = (-36, 36), select=True, idtext = '260115e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

c2.draw_centered_vs_t('QV', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (0, 0.02), xlim = (-36, 36), select=True, idtext = '260115f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c2.draw_scatter('P600', 'P500withinP600', factor2 = 1.2, carray = 'P500withinP600', filtername='WCB_Cy1_new', idtext = '260115g', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_scatter('P600', 'P400withinP600', factor2 = 1.5, carray = 'P400withinP600', filtername='WCB_Cy1_new', idtext = '260115h', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_scatter('P600', 'P300withinP600', factor2 = 2, carray = 'P300withinP600', filtername='WCB_Cy1_new', idtext = '260115i', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_scatter('P600', 'P200withinP600', factor2 = 3, carray = 'P200withinP600', filtername='WCB_Cy1_new', idtext = '260115j', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_scatter('P600', 'P100withinP600', factor2 = 6, carray = 'P100withinP600', filtername='WCB_Cy1_new', idtext = '260115k', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c1.draw_centered_vs_t('z', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (0, 15000), xlim = (-36, 36), select=True, idtext = '260115l', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')

c0.draw_centered_vs_t('z', 'WCB_Conv', 'Pcross600', plottype = 'Fill', ylim = (0, 15000), xlim = (-36, 36), select=True, idtext = '260115m', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('z', 'WCB_NonConv', 'Pcross600', plottype = 'Fill', ylim = (0, 15000), xlim = (-36, 36), select=True, idtext = '260115n', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

c2.draw_centered_vs_t('z', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (0, 15000), xlim = (-36, 36), select=True, idtext = '260115o', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

#270115

c2.draw_centered_vs_t('P', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (200, 1000), xlim = (-6, 36), select=True, idtext = '270115a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c0.draw_centered_vs_t('P_dt', 'WCB_NonConv', 'Pcross600', plottype = 'Fill', ylim = (-0.10, 0.1), xlim = (-6, 6), select=True, idtext = '270115b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('P_dt', 'WCB_Conv', 'Pcross600', plottype = 'Fill', ylim = (-0.50, 0.1), xlim = (-6, 6), select=True, idtext = '270115c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c1.draw_centered_vs_t('P_dt', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (-0.5, 0.1), xlim = (-6, 6), select=True, idtext = '270115d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c2.draw_centered_vs_t('P_dt', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (-0.1, 0.1), xlim = (-6, 6), select=True, idtext = '270115e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c0.draw_vs_p('P_dt', 'WCB_NonConv', 'P600', (100, 1000), ylim = (-0.1, 0.05), idtext = '270115f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_vs_p('P_dt', 'WCB_Conv', 'P600', (100, 1000), ylim = (-0.5, 0.1), idtext = '270115g', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c1.draw_vs_p('P_dt', 'WCB', 'P600', (100, 1000), ylim = (-0.5, 0.1), idtext = '270115h', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c2.draw_vs_p('P_dt', 'WCB_Cy1_new', 'P600', (100, 1000), ylim = (-0.1, 0.05), idtext = '270115i', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

#280115

plt.style.use('fivethirtyeight')
c1.draw_hist('P600', mintohrs = True, idtext = '280115a', filtername = 'WCB', ylog = True)
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 600hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/hist_P600_WCB_280115a'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c0.draw_hist('P600', mintohrs = True, idtext = '280115b', filtername = 'WCB', ylog = True)
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 600hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist_P600_WCB_280115b'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')

c0.draw_hist('P600', mintohrs = True, idtext = '280115c', filtername = 'WCB_Conv', ylog = True)
plt.gca().set_xlim(0, 48)
plt.xlabel('Ascent time for 600hPa [hrs]')
plt.xticks(range(0, 48 + 6, 6))
savename = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist_P600_WCB_Conv_280115c'
plt.savefig(savename, dpi = 400, bbox_inches = 'tight')


c2.draw_hist_3d(['P200withinP600', 'P300withinP600', 'P400withinP600', 'P500withinP600', 'P600'], ['WCB', 'WCB', 'WCB', 'WCB', 'WCB'], idtext = '280115e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/', ylim = (0, 3000))
c1.draw_hist_3d(['P200withinP600', 'P300withinP600', 'P400withinP600', 'P500withinP600', 'P600'], ['WCB', 'WCB', 'WCB', 'WCB', 'WCB'], idtext = '280115f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', ylim = (0, 10000))
c0.draw_hist_3d(['P200withinP600', 'P300withinP600', 'P400withinP600', 'P500withinP600', 'P600'], ['WCB', 'WCB', 'WCB', 'WCB', 'WCB'], idtext = '280115g', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', ylim = (0, 3000))

#290115

c2.draw_centered_vs_t('POT_VORTIC', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (-3, 6), xlim = (-6, 36), select=True, idtext = '290115a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_centered_vs_t('P', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (200, 1000), xlim = (-12, 36), select=True, idtext = '290115b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_centered_vs_t('z', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (0, 15000), xlim = (-12, 36), select=True, idtext = '290115c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c2.draw_centered_vs_t('THETA', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-12, 36), select=True, idtext = '290115d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_centered_vs_t('THETAE', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (280, 350), xlim = (-12, 36), select=True, idtext = '290115e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c2.draw_centered_vs_t('POT_VORTIC', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (-3, 6), xlim = (-12, 36), select=True, idtext = '290115f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')
c2.draw_centered_vs_t('QV', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (0, 0.02), xlim = (-12, 36), select=True, idtext = '290115g', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

#030215

for t in range(0, c0.maxmins, 360):
    c0.draw_trj_dot([['var4', 300], 'PMSL'], t, filtername = 'WCB_NonConv', thinning = 10, idtext = '030215a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/synop/set1_', setting = 1, path = True)
    c0.draw_contour(['var145_S', 'PMSL', 'TOT_PREC_S'], t, idtext = '030215b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/synop/set2_', setting = 2)

for t in range(0, c1.maxmins, 360):
    c1.draw_contour([['var4', 300], 'PMSL'], t, idtext = '030215c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/synop/set1_', setting = 1)
    c1.draw_contour(['var145_S', 'PMSL', 'TOT_PREC_S'], t, idtext = '030215d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/synop/set2_', setting = 2)
    
for t in range(0, c2.maxmins, 360):
    c2.draw_trj_dot([['var4', 300], 'PMSL'], t, filtername = 'WCB_Cy1_new', thinning = 10, idtext = '030215e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/synop/set1_', setting = 1, path = True)
    c2.draw_contour(['var145_S', 'PMSL', 'TOT_PREC_S'], t, idtext = '030215f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/synop/set2_', setting = 2)

#040215 - Deepening
for t in range(0, c0.maxmins, 60):
    c0.draw_contour(['PMSL'], t, idtext = 'Deepening', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/deepening/', setting = 3)
for t in range(0, c2.maxmins, 60):
    c2.draw_contour(['PMSL'], t, idtext = 'Deepening', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/deepening/', setting = 3)

#100215
c1.draw_hist_stacked(['P600_start', 'P600_start'], ['WCB_Conv', 'WCB_NonConv'], ['WCB_Conv', 'WCB_NonConv'], xlim = (0, 6120 / 60.), bins = 102 / 3 + 1, realdate = True, legpos = 2, idtext = '100215a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')

c1.draw_hist('P600', [0, 48], [0, 12000], filtername='WCB', mintohrs=True, bins = 100, ylog = False, exp = [3150, -0.072], savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', idtext='100215b')
c2.draw_hist('P600', [0, 48], [0, 200], filtername='WCB_Cy1_new', mintohrs=True, bins = 100, ylog = False, idtext = '100215c', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

#110215
c2.create_filter('Wconv', [('P600', 0, 2880), ('P300withinP600', 0, 75), ('-10201025', True)])
c2.create_filter('Wnonconv', [('P600', 0, 2880), ('P300withinP600', 75, 2880), ('-10201025', True)])
c0.draw_hist_stacked(['P600', 'P600'], ['Wconv', 'Wnonconv'], ['WCB_Conv', 'WCB_NonConv'], xlim = (0, 48), ylim = (0, 1400), bins = 100, ylog = False, exp = [200, -0.072], idtext = '110215a', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

c1.draw_hist_stacked(['P600_start', 'P600_start'], ['Wconv', 'Wnonconv'], ['WCB_Conv', 'WCB_NonConv'], xlim = (0, 6120 / 60.), bins = 102 / 3 + 1, realdate = True, legpos = 2, idtext = '110215b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', ylim = (0, 15000))
c2.draw_hist_stacked(['P600_start', 'P600_start'], ['Wconv', 'Wnonconv'], ['WCB_Conv', 'WCB_NonConv'], xlim = (0, c2.maxmins / 60.), bins = 120 / 3 + 1, realdate = True, legpos = 1, ylim = (0,2000), savebase='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/', idtext = '110215c')

#120215

c0.draw_trj_all([], 'Wconv', thinning = 10, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '120215a')
c0.draw_trj_all([], 'Wnonconv', thinning = 10, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '120215b')

#180215
c0.draw_vs_p('THETA', ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '180215a', ylim = (280, 350), ylabel = 'THETA [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', extobj = c1, legnames = ['OCT convective', 'JUL'], legpos = 2, ax2upper = 300)
c0.draw_vs_p('THETA', ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '180215b', ylim = (275, 335), ylabel = 'THETA [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', extobj = c2, legnames = ['OCT non-convective', 'JAN'], legpos = 2, ax2upper = 500)
c0.draw_vs_p('THETAE', ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '180215c', ylim = (300, 350), ylabel = 'equivalent potential temperature [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', extobj = c1, legnames = ['OCT convective', 'JUL'], legpos = 2, ax2upper = 300)
c0.draw_vs_p('THETAE', ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '180215d', ylim = (285, 335), ylabel = 'equivalent potential temperature [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', extobj = c2, legnames = ['OCT nonconvective', 'JAN'], legpos = 2, ax2upper = 500)
c0.draw_vs_p('QV', ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '180215e', ylim = (-0.005, 0.015), ylabel = 'specific humidity [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', extobj = c1, legnames = ['OCT convective', 'JUL'], legpos = 1, ax2upper = 300)
c0.draw_vs_p('QV', ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '180215f', ylim = (-0.005, 0.015), ylabel = 'specific humidity [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', extobj = c2, legnames = ['OCT non-convective', 'JAN'], legpos = 1, ax2upper = 500)
c0.draw_vs_p(['Q1', 'QC'], ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '180215g', ylim = (-0.0003, 0.001), ylabel = 'LWC [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', extobj = c2, legnames = ['OCT non-convective', 'JAN'], legpos = 1, ax2upper = 500)
c0.draw_vs_p(['Q3', 'QS'], ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '180215h', ylim = (-0.0003, 0.001), ylabel = 'SWC [kg/kg]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', extobj = c2, legnames = ['OCT non-convective', 'JAN'], legpos = 1, ax2upper = 500)
c0.draw_vs_p(['Q2', 'QI'], ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '180215i', ylim = (-0.00004, 0.0001), ylabel = 'IWC [kg/kg]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', extobj = c2, legnames = ['OCT non-convective', 'JAN'], legpos = 1, ax2upper = 500)
c0.draw_vs_p(['var4', 'POT_VORTIC'], ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '180215j', ylim = (-5, 5), ylabel = 'potential vorticity [PVU]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', extobj = c2, legnames = ['OCT non-convective', 'JAN'], legpos = 1, ax2upper = 500)
c0.draw_vs_p(['var4', 'POT_VORTIC'], ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '180215k', ylim = (-7, 7), ylabel = 'potential vorticity [PVU]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', extobj = c1, legnames = ['OCT convective', 'JUL'], legpos = 2, ax2upper = 300)


#230215

c1.draw_contour([ 'PMSL', ['var4', 300]], 1440, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/230215a_')
c1.draw_contour([ 'PMSL', ['var4', 300]], 2880, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/230215b_')
c1.draw_contour([ 'PMSL', ['var4', 300]], 4320, setting = 1, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/230215c_')

c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 1440, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/230215d_')
c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 2880, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/230215e_')
c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 4320, setting = 2, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/230215f_')

c0.draw_trj_dot([['var4', 300], 'PMSL'], 1440, filtername = 'WCB_NonConv', thinning = 10, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/230215g_')
c0.draw_trj_dot([['var4', 300], 'PMSL'], 2880, filtername = 'WCB_NonConv', thinning = 10, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/230215h_')
c0.draw_trj_dot([['var4', 300], 'PMSL'], 4320, filtername = 'WCB_NonConv', thinning = 10, setting = 1, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/230215i_')

c0.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 1440, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/230215j_')
c0.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 2880, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/230215k_')
c0.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 4320, setting = 2, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/230215l_')

c2.draw_trj_dot([['var4', 300], 'PMSL'], 720, filtername = 'WCB_Cy1_new', thinning = 10, setting = 1, cbar = False,
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/230215m_')
c2.draw_trj_dot([['var4', 300], 'PMSL'], 2160, filtername = 'WCB_Cy1_new', thinning = 10, setting = 1, cbar = False,
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/230215n_')
c2.draw_trj_dot([['var4', 300], 'PMSL'], 3600, filtername = 'WCB_Cy1_new', thinning = 10, setting = 1, cbar = True,
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/230215o_')

c2.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 720, setting = 5, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/230215p_')
c2.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 2160, setting = 5, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/230215q_')
c2.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 3600, setting = 5, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/230215r_')

#030315
c0.draw_vs_p(['Q1', 'QC'], ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '030315a', ylim = (-0.0003, 0.001), ylabel = 'LWC [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', extobj = c1, legnames = ['OCT non-convective', 'JAN'], legpos = 1, ax2upper = 500)
c0.draw_vs_p(['Q3', 'QS'], ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '030315b', ylim = (-0.0003, 0.001), ylabel = 'SWC [kg/kg]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', extobj = c1, legnames = ['OCT non-convective', 'JAN'], legpos = 1, ax2upper = 500)
c0.draw_vs_p(['Q2', 'QI'], ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '030315c', ylim = (-0.00004, 0.0001), ylabel = 'IWC [kg/kg]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', extobj = c1, legnames = ['OCT non-convective', 'JAN'], legpos = 1, ax2upper = 500)

c0.draw_centered_vs_t('P', ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (200, 1300), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', legnames = ['OCT convective', 'JUL'], idtext = '030315d')
c0.draw_centered_vs_t('P', ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (200, 1300), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', legnames = ['OCT non-convective', 'JAN'], idtext = '030315e')

c0.draw_centered_vs_t('THETA', ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (270, 340), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', legnames = ['OCT convective', 'JUL'], idtext = '030315f')
c0.draw_centered_vs_t('THETA', ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (270, 340), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', legnames = ['OCT non-convective', 'JAN'], idtext = '030315g')

c0.draw_centered_vs_t('THETAE', ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (270, 340), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', legnames = ['OCT convective', 'JUL'], idtext = '030315h')
c0.draw_centered_vs_t('THETAE', ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (270, 340), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', legnames = ['OCT non-convective', 'JAN'], idtext = '030315i')

c0.draw_centered_vs_t('QV', ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-0.005, 0.015), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', legnames = ['OCT convective', 'JUL'], idtext = '030315j')
c0.draw_centered_vs_t('QV', ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-0.005, 0.015), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', legnames = ['OCT non-convective', 'JAN'], idtext = '030315k')

c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-6, 6), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', legnames = ['OCT convective', 'JUL'], idtext = '030315l')
c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-6, 6), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', legnames = ['OCT non-convective', 'JAN'], idtext = '030315m')

c0.draw_centered_vs_t(['Q1', 'QC'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-0.0005, 0.0015), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', legnames = ['OCT convective', 'JUL'], idtext = '030315n')
c0.draw_centered_vs_t(['Q1', 'QC'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-0.0005, 0.0015), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', legnames = ['OCT non-convective', 'JAN'], idtext = '030315o')

c0.draw_centered_vs_t(['Q3', 'QS'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-0.0005, 0.0015), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', legnames = ['OCT convective', 'JUL'], idtext = '030315p')
c0.draw_centered_vs_t(['Q3', 'QS'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-0.0005, 0.0015), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', legnames = ['OCT non-convective', 'JAN'], idtext = '030315q')

c0.draw_centered_vs_t(['Q2', 'QI'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-0.00005, 0.0001), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_', legnames = ['OCT convective', 'JUL'], idtext = '030315r')
c0.draw_centered_vs_t(['Q2', 'QI'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-0.00005, 0.0001), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_', legnames = ['OCT non-convective', 'JAN'], idtext = '030315s')

#100315
c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-6, 6), xlim = (-6, 72), select=(-6, 36), legnames = ['OCT non-convective', 'JAN'], idtext = '100315a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/nonconv_')
c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-6, 6), xlim = (-12, 72), select=(-6, 36), idtext = '100315b', legnames = ['OCT convective', 'JUL'], savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/conv_')

#120315
returnlist = c0.draw_vs_p('THETA', ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (280, 350), ylabel = r'$\theta\,[\mathrm{K}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/120315a_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None)
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (10, 3))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
zero = np.zeros(returnlist[1].shape[0])
for n in range(len(returnlist)-1):
    trj.plots.fill_between_steps(returnlist[0], np.array(returnlist[n+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 500
inc = 100
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('p [hPa]')
ax.set_xlim(200, 1000)
ax.invert_xaxis()
plt.subplots_adjust(left = 0.075, right = 0.95, top = 0.95, bottom = 0.2)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/120315b_p_density.png')
    
returnlist = c0.draw_vs_p('THETA', ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (280, 350), ylabel = r'$\theta\,[\mathrm{K}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/120315c', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None)  
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (10, 3))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
zero = np.zeros(returnlist[1].shape[0])
for n in range(len(returnlist)-1):
    trj.plots.fill_between_steps(returnlist[0], np.array(returnlist[n+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 500
inc = 100
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('p [hPa]')
ax.set_xlim(200, 1000)
ax.invert_xaxis()
plt.subplots_adjust(left = 0.075, right = 0.95, top = 0.95, bottom = 0.2)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/120315d_p_density.png')    
    
c0.draw_vs_p('THETAE', ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (290, 350), ylabel = r'$\theta_e\,[\mathrm{K}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/120315e_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None, letter = '(d)')    
    
c0.draw_vs_p('THETAE', ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (290, 350), ylabel = r'$\theta_e\,[\mathrm{K}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/120315f', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '(c)')    
    
c0.draw_vs_p('QV', ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (0, 15), ylabel = r'$q\,[\mathrm{g\,kg^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/120315g_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None, letter = '(f)', mult = 1000)    
    
c0.draw_vs_p('QV', ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (0, 15), ylabel = r'$q\,[\mathrm{g\,kg^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/120315h', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '(e)', mult = 1000)      
    
c0.draw_vs_p(['var4', 'POT_VORTIC'], ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (-6, 6), ylabel = 'PV [pvu]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/120315i_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None, letter = '(b)')    
    
c0.draw_vs_p(['var4', 'POT_VORTIC'], ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (-6, 6), ylabel = 'PV [pvu]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/120315j', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '(a)')     

c0.draw_vs_p(['Q1', 'QC'], ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (0, 1000), ylabel = r'$\mathrm{LWC}\,[\mathrm{mg\,kg^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/120315k_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None, letter = '(d)', mult = 1000000)    
    
c0.draw_vs_p(['Q1', 'QC'], ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (0, 1000), ylabel = r'$\mathrm{LWC}\,[\mathrm{mg\,kg^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/120315l', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '(c)', mult = 1000000)
    
c0.draw_vs_p(['Q3', 'QS'], ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (0, 1000), ylabel = r'$\mathrm{SWC}\,[\mathrm{mg\,kg^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/120315m_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None, letter = '(f)', mult = 1000000)    
    
c0.draw_vs_p(['Q3', 'QS'], ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (0, 1000), ylabel = r'$\mathrm{SWC}\,[\mathrm{mg\,kg^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/120315n', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '(e)', mult = 1000000)   
    
c0.draw_vs_p(['Q2', 'QI'], ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (0, 100), ylabel = r'$\mathrm{IWC}\,[\mathrm{mg\,kg^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/120315o_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None, letter = '(h)', mult = 1000000)    
    
c0.draw_vs_p(['Q2', 'QI'], ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (0, 100), ylabel = r'$\mathrm{IWC}\,[\mathrm{mg\,kg^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/120315p', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '(g)', mult = 1000000)     

#130315 vs_t plots 

returnlist = c0.draw_centered_vs_t('P', ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (200, 1000), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130315a_', legnames = ['OCTnc', 'JAN'], letter = '(b)', ax2 = None, ylabel = 'p [hPa]')
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (10, 3))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
for n in range(len(returnlist) / 2):
    zero = np.zeros(returnlist[n * 2].shape[0])
    trj.plots.fill_between_steps(returnlist[n * 2], np.array(returnlist[n*2+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 100
inc = 20
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.xaxis.set_ticks(np.arange(-120, 120, 6))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('time [h] relative to center')
ax.set_xlim(-36, 36)
plt.subplots_adjust(left = 0.075, right = 0.95, top = 0.95, bottom = 0.2)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130315b_t_density.png')    

returnlist = c0.draw_centered_vs_t('P', ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (200, 1000), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130315c_', legnames = ['OCTc', 'JUL'], letter = '(a)', ax2 = None, ylabel = 'p [hPa]')
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (10, 3))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
for n in range(len(returnlist) / 2):
    zero = np.zeros(returnlist[n * 2].shape[0])
    trj.plots.fill_between_steps(returnlist[n * 2], np.array(returnlist[n*2+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 100
inc = 20
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.xaxis.set_ticks(np.arange(-120, 120, 6))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('time [h] relative to center')
ax.set_xlim(-12, 12)
plt.subplots_adjust(left = 0.075, right = 0.95, top = 0.95, bottom = 0.2)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130315d_t_density.png')  

c0.draw_centered_vs_t('THETA', ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130315e_', legnames = ['OCTnc', 'JAN'], letter = '(d)', ax2 = None, ylabel = r'$\theta\,[\mathrm{K}]$')

c0.draw_centered_vs_t('THETA', ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (280, 350), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130315f_', legnames = ['OCTc', 'JUL'], letter = '(c)', ax2 = None, ylabel = r'$\theta\,[\mathrm{K}]$')

c0.draw_centered_vs_t('THETAE', ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (290, 350), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130315g_', legnames = ['OCTnc', 'JAN'], letter = '(f)', ax2 = None, ylabel = r'$\theta_e\,[\mathrm{K}]$')

c0.draw_centered_vs_t('THETAE', ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (290, 350), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130315h_', legnames = ['OCTc', 'JUL'], letter = '(e)', ax2 = None, ylabel = r'$\theta_e\,[\mathrm{K}]$')

c0.draw_centered_vs_t('QV', ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (0, 15), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130315i_', legnames = ['OCTnc', 'JAN'], letter = '(h)', ax2 = None, ylabel = r'$q\,[\mathrm{g\,kg^{-1}}]$', mult = 1000)

c0.draw_centered_vs_t('QV', ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (0, 15), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130315j_', legnames = ['OCTc', 'JUL'], letter = '(g)', ax2 = None, ylabel = r'$q\,[\mathrm{g\,kg^{-1}}]$', mult = 1000)

c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-6, 6), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130315k_', legnames = ['OCTnc', 'JAN'], letter = '(b)', ax2 = None, ylabel = 'PV [pvu]', mult = 1)

c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-6, 6), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130315l_', legnames = ['OCTc', 'JUL'], letter = '(a)', ax2 = None, ylabel = 'PV [pvu]', mult = 1)

c0.draw_centered_vs_t(['Q1', 'QC'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (0, 1000), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130315m_', legnames = ['OCTnc', 'JAN'], letter = '(d)', ax2 = None, ylabel = r'$\mathrm{LWC}\,[\mathrm{mg\,kg^{-1}}]$', mult = 1000000)

c0.draw_centered_vs_t(['Q1', 'QC'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (0, 1000), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130315n_', legnames = ['OCTc', 'JUL'], letter = '(c)', ax2 = None, ylabel = r'$\mathrm{LWC}\,[\mathrm{mg\,kg^{-1}}]$', mult = 1000000)

c0.draw_centered_vs_t(['Q3', 'QS'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (0, 1500), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130315o_', legnames = ['OCTnc', 'JAN'], letter = '(f)', ax2 = None, ylabel = r'$\mathrm{SWC}\,[\mathrm{mg\,kg^{-1}}]$', mult = 1000000)

c0.draw_centered_vs_t(['Q3', 'QS'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (0, 1500), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130315p_', legnames = ['OCTc', 'JUL'], letter = '(e)', ax2 = None, ylabel = r'$\mathrm{SWC}\,[\mathrm{mg\,kg^{-1}}]$', mult = 1000000)

c0.draw_centered_vs_t(['Q2', 'QI'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (0, 100), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130315q_', legnames = ['OCTnc', 'JAN'], letter = '(h)', ax2 = None, ylabel = r'$\mathrm{IWC}\,[\mathrm{mg\,kg^{-1}}]$', mult = 1000000)

c0.draw_centered_vs_t(['Q2', 'QI'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (0, 100), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130315r_', legnames = ['OCTc', 'JUL'], letter = '(g)', ax2 = None, ylabel = r'$\mathrm{IWC}\,[\mathrm{mg\,kg^{-1}}]$', mult = 1000000)

#190315
c1.draw_hist('P600', [0, 48], [0, 12000], filtername='WCB', mintohrs=True, bins = 100, ylog = False, exp = [3150, -0.072], savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', idtext='')

c2.draw_hist('P600', [0, 48], [0, 200], filtername='WCB_Cy1_new', mintohrs=True, bins = 100, ylog = False, idtext = '', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c0.draw_hist_stacked(['P600', 'P600'], ['Wconv', 'Wnonconv'], ['OCTc', 'OCTnc'], xlim = (0, 48), ylim = (0, 1400), bins = 100, ylog = False, exp = [200, -0.072], idtext = '', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

returnlist = c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wnonconv', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-6, 6), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/190315a', legnames = ['OCTnc', 'JAN'], letter = '(b)', ax2 = None, ylabel = 'PV [PVU]')
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (10, 3))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
for n in range(len(returnlist) / 2):
    zero = np.zeros(returnlist[n * 2].shape[0])
    trj.plots.fill_between_steps(returnlist[n * 2], np.array(returnlist[n*2+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 100
inc = 20
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.xaxis.set_ticks(np.arange(-120, 120, 6))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('time [h] relative to center')
ax.set_xlim(-36, 36)
plt.subplots_adjust(left = 0.075, right = 0.95, top = 0.95, bottom = 0.2)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/190315b_t_density_long.png')    

returnlist = c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wconv', 'WCB'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-6, 6), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/190315c_', legnames = ['OCTc', 'JUL'], letter = '(a)', ax2 = None, ylabel = 'PV [PVU]')
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (10, 3))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
for n in range(len(returnlist) / 2):
    zero = np.zeros(returnlist[n * 2].shape[0])
    trj.plots.fill_between_steps(returnlist[n * 2], np.array(returnlist[n*2+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 100
inc = 20
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.xaxis.set_ticks(np.arange(-120, 120, 6))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('time [h] relative to center')
ax.set_xlim(-12, 12)
plt.subplots_adjust(left = 0.075, right = 0.95, top = 0.95, bottom = 0.2)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/190315d_t_density_long.png') 

c0.draw_trj_dot(["PMSL", "TOT_PREC_S"], tplus = 3600, filtername=['Wconv', 'Wnonconv'], thinning = 10, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/dot/')

#190315
#c0.draw_centered_vs_t('P_dt', 'WCB_NonConv', 'Pcross600', plottype = 'Fill', ylim = (-0.10, 0.1), xlim = (-6, 6), select=True, idtext = '190315e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('P_dt', 'WCB_Conv', 'Pcross600', plottype = 'Fill', ylim = (-0.50, 0.1), xlim = (-6, 6), select=True, idtext = '190315f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c1.draw_centered_vs_t('P_dt', 'WCB', 'Pcross600', plottype = 'Fill', ylim = (-0.5, 0.1), xlim = (-6, 6), select=True, idtext = '190315g', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
#c2.draw_centered_vs_t('P_dt', 'WCB_Cy1_new', 'Pcross600', plottype = 'Fill', ylim = (-0.1, 0.1), xlim = (-6, 6), select=True, idtext = '190315h', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/')

c0.draw_vs_p('P_dt', ['Wnonconv', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (-0.1, 0.05), ylabel = r'$\omega\,[\mathrm{hPa\,s^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/190315e_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None, letter = '(b)')    
    
c0.draw_vs_p('P_dt', ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (-0.1, 0.05), ylabel = r'$\omega\,[\mathrm{hPa\,s^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/190315f_', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '(a)')    

for i in range(10):
    for tracer in ['P', 'QV', 'QC', 'QS', 'QI', 'POT_VORTIC']:
        plt.close('all')
        c2.draw_spaghetti(tracer, 'WCB_pick4', limit = [i * 6, i * 6 + 6])
        plt.xlim(0, 40)
        plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/spaghetti/'+str(i)+tracer)

#200315
c2.draw_centered_vs_t('P', ['WCB_pick6', 'WCB_pick8'], 'Pcross600', plottype = 'Fill', ylim = (200, 1000), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/200315a_', legnames = ['WCB_pick6', 'WCB_pick8'], ylabel = 'p [hPa]')

c2.draw_centered_vs_t('POT_VORTIC', ['WCB_pick6', 'WCB_pick8'], 'Pcross600', plottype = 'Fill', ylim = (-6, 6), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/200315b_', legnames = ['WCB_pick6', 'WCB_pick8'], ylabel = 'PV [pvu]')

for i in range(10):
    for tracer in ['P', 'QV', 'QC', 'QS', 'QI', 'POT_VORTIC']:
        plt.close('all')
        c2.draw_spaghetti(tracer, 'WCB_pick6', limit = [i * 6, i * 6 + 6])
        plt.xlim(0, 40)
        plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/spaghetti/WCB_pick6'+str(i)+tracer)
for i in range(10):
    for tracer in ['P', 'QV', 'QC', 'QS', 'QI', 'POT_VORTIC']:
        plt.close('all')
        c2.draw_spaghetti(tracer, 'WCB_pick8', limit = [i * 6, i * 6 + 6])
        plt.xlim(0, 40)
        plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/spaghetti/WCB_pick8'+str(i)+tracer)

for i in range(10):
    for tracer in ['P', 'QV', 'QC', 'QS', 'QI', 'POT_VORTIC']:
        plt.close('all')
        c2.draw_spaghetti(tracer, 'WCB', limit = [i * 6, i * 6 + 6])
        plt.xlim(0, 40)
        plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/spaghetti/WCB'+str(i)+tracer)



#250315

c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 720, filtername = 'WCB_Cy1_new', thinning = 10, setting = 1, diffobj = c2eu, diffvar = 'var4', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/250315a_diff_', ctracer = 'THETA')
c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 2160, filtername = 'WCB_Cy1_new', thinning = 10, setting = 1, diffobj = c2eu, diffvar = 'var4', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/250315b_diff_', ctracer = 'THETA')
c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 3600, filtername = 'WCB_Cy1_new', thinning = 10, setting = 1, diffobj = c2eu, diffvar = 'var4', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/250315c_diff_', ctracer = 'THETA')

c2eu.draw_contour([['var4', 310, 'THETA'], 'PMSL'], 720+360, setting = 1, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/250315d_eu_')
c2eu.draw_contour([['var4', 310, 'THETA'], 'PMSL'], 2160+360, setting = 1, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/250315e_eu_')
c2eu.draw_contour([['var4', 310, 'THETA'], 'PMSL'], 3600+360, setting = 1, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/250315f_eu_')

c2.draw_contour([['var4', 310, 'THETA'], 'PMSL'], 720, setting = 1, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/250315d_de_')
c2.draw_contour([['var4', 310, 'THETA'], 'PMSL'], 2160, setting = 1, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/250315e_de_')
c2.draw_contour([['var4', 310, 'THETA'], 'PMSL'], 3600, setting = 1, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/250315f_de_')

c1.draw_trj_dot([['var4', 330, 'THETA'], 'PMSL'], 1440, filtername = 'WCB', thinning = 20, setting = 1, diffobj = c1eu, diffvar = 'var4')

#300315
c0.draw_vs_p('THETA_dt', ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (-5, 30), ylabel = 'DHR [K/h]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/310315a_', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, mult = 3600)

c0.draw_vs_p('THETA_dt', ['Wconv', 'WCB'], 'P600', (200, 1000), idtext = '', ylim = (-5, 30), ylabel = 'DHR [K/h]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/310315b_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None, mult = 3600)


c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 2160, filtername = 'WCB_Cy1_new', thinning = 2, setting = 1, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/300315a_', ctracer = 'THETA')

c2eu.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 2160+360, filtername = 'WCB_Cy1_cropped', thinning = 2, setting = 1, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/300315b_eu_', ctracer = 'THETA')

#020415
# New synop PV plots

c0.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 1440, filtername = 'WCB', thinning = 10, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/020415a_')
c0.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 2880, filtername = 'WCB', thinning = 10, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/020415b_')
c0.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 4320, filtername = 'WCB', thinning = 10, setting = 1, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/020415c_')

c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 720, filtername = 'WCB_Cy1_new', thinning = 10, setting = 1, cbar = False,
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/020415d_')
c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 2160, filtername = 'WCB_Cy1_new', thinning = 10, setting = 1, cbar = False,
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/020415e_')
c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 3600, filtername = 'WCB_Cy1_new', thinning = 10, setting = 1, cbar = True,
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/020415f_')


#070415
c1new.draw_hist('P600', [0, 48], [0, 5000], filtername='WCB_Cy1', mintohrs=True, bins = 96, exp = [1500, -0.105], ylog = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/070415a_')
c1new.draw_hist('P600', [0, 48], [0, 500], filtername='WCB_start2', mintohrs=True, bins = 96, exp = [140, -0.085], ylog = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/070415b_')

c1eu.draw_hist('P600', [0, 48], [0, 5000], filtername='WCB_Cy1', mintohrs=True, bins = 96, exp = [1500, -0.105], ylog = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/070415c_eu_')
c1eu.draw_hist('P600', [0, 48], [0, 500], filtername='WCB_start2', mintohrs=True, bins = 96, exp = [140, -0.085], ylog = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/070415d_eu_')

c0.draw_trj_all([], 'Wconv2', thinning = 10, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/070415e_', idtext = '')
c0.draw_trj_all([], 'Wnonconv2', thinning = 10, savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/070415f_', idtext = '')

# vs_t
c1 = trj.loadme('c1new.trj')

returnlist = c0.draw_centered_vs_t('P', ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (200, 1000), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/070415g_', legnames = ['OCTnc', 'JAN'], letter = '(b)', ax2 = None, ylabel = 'p [hPa]')
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (3.15, 1))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
for n in range(len(returnlist) / 2):
    zero = np.zeros(returnlist[n * 2].shape[0])
    trj.plots.fill_between_steps(returnlist[n * 2], np.array(returnlist[n*2+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 100
inc = 20
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.xaxis.set_ticks(np.arange(-120, 120, 6))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('time [h] relative to center')
ax.set_xlim(-36, 36)
plt.subplots_adjust(left = 0.17, right = 0.95, top = 0.95, bottom = 0.4)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/070415h_t_density.png')    

returnlist = c0.draw_centered_vs_t('P', ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (200, 1000), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/070415i_', legnames = ['OCTc', 'JUL'], letter = '(a)', ax2 = None, ylabel = 'p [hPa]')
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (3.15, 1))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
for n in range(len(returnlist) / 2):
    zero = np.zeros(returnlist[n * 2].shape[0])
    trj.plots.fill_between_steps(returnlist[n * 2], np.array(returnlist[n*2+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 100
inc = 20
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.xaxis.set_ticks(np.arange(-120, 120, 12))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('time [h] relative to center')
ax.set_xlim(-12, 12)
plt.subplots_adjust(left = 0.17, right = 0.95, top = 0.95, bottom = 0.4)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/070415j_t_density.png')  

c0.draw_centered_vs_t('THETA', ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (280, 350), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/070415k_', legnames = None, letter = '(d)', ax2 = None, ylabel = r'$\theta$ [K]')

c0.draw_centered_vs_t('THETA', ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (280, 350), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/070415l_', legnames = None, letter = '(c)', ax2 = None, ylabel = r'$\theta$ [K]')

c0.draw_centered_vs_t('THETAE', ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (290, 350), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/070415m_', legnames = None, letter = '(f)', ax2 = None, ylabel = r'$\theta_e$ [K]')

c0.draw_centered_vs_t('THETAE', ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (290, 350), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/070415n_', legnames = None, letter = '(e)', ax2 = None, ylabel = r'$\theta_e$ [K]')

c0.draw_centered_vs_t('QV', ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (0, 15), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/070415o_', legnames = None, letter = '(h)', ax2 = None, ylabel = r'$q$ [g kg$^{-1}$]', mult = 1000)

c0.draw_centered_vs_t('QV', ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (0, 15), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/070415p_', legnames = None, letter = '(g)', ax2 = None, ylabel = r'$q$ [g kg$^{-1}$]', mult = 1000)

c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-6, 6), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/070415q_', legnames = ['OCTnc', 'JAN'], letter = '(b)', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3)

c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-6, 6), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/070415r_', legnames = ['OCTc', 'JUL'], letter = '(a)', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 3)

c0.draw_centered_vs_t(['Q1', 'QC'], ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (0, 1000), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/070415s_', legnames = None, letter = '(d)', ax2 = None, ylabel = r'LWC [mg kg$^{-1}$]', mult = 1000000)

c0.draw_centered_vs_t(['Q1', 'QC'], ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (0, 1000), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/070415t_', legnames = None, letter = '(c)', ax2 = None, ylabel = r'LWC [mg kg$^{-1}$]', mult = 1000000)

c0.draw_centered_vs_t(['Q3', 'QS'], ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (0, 1500), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/070415u_', legnames = None, letter = '(f)', ax2 = None, ylabel = r'SWC [mg kg$^{-1}$]', mult = 1000000)

c0.draw_centered_vs_t(['Q3', 'QS'], ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (0, 1500), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/070415v_', legnames = None, letter = '(e)', ax2 = None, ylabel = r'SWC [mg kg$^{-1}$]', mult = 1000000)

c0.draw_centered_vs_t(['Q2', 'QI'], ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (0, 100), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/070415w_', legnames = None, letter = '(h)', ax2 = None, ylabel = r'IWC [mg kg$^{-1}$]', mult = 1000000)

c0.draw_centered_vs_t(['Q2', 'QI'], ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (0, 100), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/070415x_', legnames = None, letter = '(g)', ax2 = None, ylabel = r'IWC [mg kg$^{-1}$]', mult = 1000000)


# JUL synop
c1.draw_trj_dot([ 'PMSL', ['var4', 325, 'THETA']], 900, filtername = 'WCB_Cy1', thinning = 10, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415y_')
c1.draw_trj_dot([ 'PMSL', ['var4', 325, 'THETA']], 2340, filtername = 'WCB_Cy1', thinning = 10, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415z_')
c1.draw_trj_dot([ 'PMSL', ['var4', 325, 'THETA']], 3780, filtername = 'WCB_Cy1', thinning = 10, setting = 1, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415aa_')

c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 900, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415ab_')
c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 2340, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415ac_')
c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 3780, setting = 2, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415ad_')

#OCThist
c0.draw_hist_stacked(['P600', 'P600'], ['Wconv2', 'Wnonconv2'], ['WCB_Conv', 'WCB_NonConv'], xlim = (0, 48), ylim = (0, 1400), bins = 96, ylog = False, exp = [200, -0.072], idtext = '', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/070415ae_')


#080415
returnlist = c0.draw_vs_p('THETA', ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (280, 350), ylabel = r'$\theta$ [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/080415a_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None)
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (3.15, 1))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
zero = np.zeros(returnlist[1].shape[0])
for n in range(len(returnlist)-1):
    trj.plots.fill_between_steps(returnlist[0], np.array(returnlist[n+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 500
inc = 100
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('p [hPa]')
ax.set_xlim(200, 1000)
ax.invert_xaxis()
plt.subplots_adjust(left = 0.17, right = 0.95, top = 0.95, bottom = 0.4)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/080415b_p_density.png')
    
returnlist = c0.draw_vs_p('THETA', ['Wconv', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (280, 350), ylabel = r'$\theta$ [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/080415c', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None)  
color = []
color.append(['#FF9999', '#FF3333', '#CC0000', '#FF3333'])
#color.append(['#8080FF', '#3333FF', '#0000CC', '#0000CC'])
color.append(['#E6E6E6', '#B8B8B8', '#000000', '#666666'])
hatch = ['/', '\\']
fig = plt.figure(figsize = (3.15, 1))
ax = plt.gca()
ax.grid(color = 'dimgrey', linestyle = '-')
plt.tick_params(axis = 'both', which = 'both', bottom = 'off', top = 'off',
                    left = 'off', right = 'off')
zero = np.zeros(returnlist[1].shape[0])
for n in range(len(returnlist)-1):
    trj.plots.fill_between_steps(returnlist[0], np.array(returnlist[n+1]), zero, 
                             ax = ax, linewidth = 1, 
                alpha = 0.5, edgecolor = color[n][2], 
                color = 'none', hatch = hatch[n])
maxbin = 500
inc = 100
ax.set_yticks(np.arange(inc, maxbin + inc, inc))
ax.set_ylim((0, maxbin))
ax.set_ylabel('%')
ax.set_xlabel('p [hPa]')
ax.set_xlim(200, 1000)
ax.invert_xaxis()
plt.subplots_adjust(left = 0.17, right = 0.95, top = 0.95, bottom = 0.4)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/080415d_p_density.png')    
    
c0.draw_vs_p('THETAE', ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (290, 350), ylabel = r'$\theta_e$ [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/080415e_', extobj = c2, legnames = None, legpos = 2, ax2upper = None, letter = '(d)')    
    
c0.draw_vs_p('THETAE', ['Wconv2', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (290, 350), ylabel = r'$\theta_e$ [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/080415f', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '(c)')    
    
c0.draw_vs_p('QV', ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (0, 15), ylabel = r'$q$ [g kg$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/080415g_', extobj = c2, legnames = None, legpos = 2, ax2upper = None, letter = '(f)', mult = 1000)    
    
c0.draw_vs_p('QV', ['Wconv2', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (0, 15), ylabel = r'$q$ [g kg$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/080415h', extobj = c1, legnames = None, legpos = 2, ax2upper = None, letter = '(e)', mult = 1000)      
    
c0.draw_vs_p(['var4', 'POT_VORTIC'], ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (-6, 6), ylabel = 'PV [pvu]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/080415i_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 3, ax2upper = None, letter = '(b)')    
    
c0.draw_vs_p(['var4', 'POT_VORTIC'], ['Wconv2', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (-6, 6), ylabel = 'PV [pvu]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/080415j', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 3, ax2upper = None, letter = '(a)')     

c0.draw_vs_p(['Q1', 'QC'], ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (0, 1000), ylabel = r'LWC [mg kg$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/080415k_', extobj = c2, legnames = None, legpos = 2, ax2upper = None, letter = '(d)', mult = 1000000)    
    
c0.draw_vs_p(['Q1', 'QC'], ['Wconv2', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (0, 1000), ylabel = r'LWC [mg kg$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/080415l', extobj = c1, legnames = None, legpos = 2, ax2upper = None, letter = '(c)', mult = 1000000)
    
c0.draw_vs_p(['Q3', 'QS'], ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (0, 1000), ylabel = r'SWC [mg kg$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/080415m_', extobj = c2, legnames = None, legpos = 2, ax2upper = None, letter = '(f)', mult = 1000000)    
    
c0.draw_vs_p(['Q3', 'QS'], ['Wconv2', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (0, 1000), ylabel = r'SWC [mg kg$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/080415n', extobj = c1, legnames = None, legpos = 2, ax2upper = None, letter = '(e)', mult = 1000000)   
    
c0.draw_vs_p(['Q2', 'QI'], ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (0, 100), ylabel = r'IWC [mg kg$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/080415o_', extobj = c2, legnames = None, legpos = 2, ax2upper = None, letter = '(h)', mult = 1000000)    
    
c0.draw_vs_p(['Q2', 'QI'], ['Wconv2', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (0, 100), ylabel = r'IWC [mg kg$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/080415p', extobj = c1, legnames = None, legpos = 2, ax2upper = None, letter = '(g)', mult = 1000000)     

#090415
c1.draw_scatter('P_max', 'P600', carray='P100withinP600', filtername = 'WCB_Cy1', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/090415a_')

c2.draw_scatter('P_max', 'P600', carray='P100withinP600', filtername = 'WCB_Cy1_new', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/090415b_')

c0.draw_scatter('P_max', 'P600', carray='P600', filtername = 'Wconv2', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/090415c_')
c0.draw_scatter('P_max', 'P600', carray='P600', filtername = 'Wnonconv2', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/090415d_')



c2.draw_hist('P_max', (0, 1), (0, 1000), filtername='WCB_Cy1_new', bins = 100, xlabel = 'Maximum ascent time [hPa/s]', legend = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/090415e_')

c1.draw_hist('P_max', (0, 1), (0, 1000), filtername='WCB_Cy1', bins = 100, xlabel = 'Maximum ascent time [hPa/s]', legend = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/090415f_')

c0.draw_hist('P_max', (0, 1), (0, 1000), filtername='Wnonconv2', bins = 100, xlabel = 'Maximum ascent time [hPa/s]', legend = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/090415g_')
c0.draw_hist('P_max', (0, 1), (0, 500), filtername='Wconv2', bins = 100, xlabel = 'Maximum ascent time [hPa/s]', legend = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/090415h_')


c0.draw_vs_p('P_dt', ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (-0.05, 0.1), ylabel = r'$\omega\,[\mathrm{hPa\,s^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/090415i_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 2, ax2upper = None, letter = '', mult = -1)    
    
c0.draw_vs_p('P_dt', ['Wconv2', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (-0.1, 0.5), ylabel = r'$\omega\,[\mathrm{hPa\,s^{-1}}]$', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/0', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '', mult = -1) 



#100415
c1.draw_spaghetti('P', 'WCB_Cy1', limit = (0, 20), locind = 14, xlim = (25, 45), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/100415a_')

c0.draw_spaghetti('P', 'Wconv2', limit = (0, 20), locind = 130, xlim = (60, 80), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/100415b_')
c0.draw_spaghetti('P', 'Wnonconv2', limit = (0, 20), locind = 36, xlim = (35, 75), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/100415c_')

c2.draw_spaghetti('P', 'WCB_Cy1_new', limit = (0, 20), locind = 1, xlim = (5, 35), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/100415d_')


#130415
c1.draw_hist('P600', [0, 48], [0, 5000], filtername='WCB_Cy1', mintohrs=True, bins = 96, ylog = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130415a_')
c2.draw_hist('P600', [0, 48], [0, 250], filtername='WCB_Cy1_new', mintohrs=True, bins = 96, ylog = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/130415b_')
c0.draw_hist_stacked(['P600', 'P600'], ['Wconv2', 'Wnonconv2'], ['OCTc', 'OCTnc'], xlim = (0, 48), ylim = (0, 1400), bins = 96, ylog = False, idtext = '', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/130415c_')
c0.draw_hist('P600', [0, 48], [0, 300], filtername='Wnonconv2', mintohrs=True, bins = 96, ylog = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/130415d_')

c1.draw_hist('P600', [0, 5], [0, 1500], filtername='WCB_Cy1', mintohrs=True, bins = 60, ylog = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130415e_')
c0.draw_hist('P600', [0, 5], [0, 350], filtername='Wconv2', mintohrs=True, bins = 60, ylog = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/130415f_')

c1.draw_hist('P600', [0, 48], [1, 10000], filtername='WCB_Cy1', mintohrs=True, bins = 96, ylog = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/130415g_', exp = [1550, -0.105])
c0.draw_hist('P600', [0, 48], [1, 10000], filtername='Wconv2', mintohrs=True, bins = 96, ylog = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/130415h_', exp = [300, -0.09])

c0.draw_centered_vs_t('latitude', ['Wconv2', 'Wnonconv2'], 'Pcross600', plottype = 'Fill', ylim = (30, 60), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/130415l_', legnames = ['OCTc', 'OCTnc'], ax2 = None, ylabel = 'latitude [deg]')


### Synop

c1.draw_trj_dot([ 'PMSL', ['var4', 325, 'THETA']], 900, filtername = 'WCB_Cy1', thinning = 20, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415y_')
c1.draw_trj_dot([ 'PMSL', ['var4', 325, 'THETA']], 2340, filtername = 'WCB_Cy1', thinning = 20, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415z_')
c1.draw_trj_dot([ 'PMSL', ['var4', 325, 'THETA']], 3780, filtername = 'WCB_Cy1', thinning = 20, setting = 1, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415aa_')

#c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 900, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415ab_')
#c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 2340, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415ac_')
#c1.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 3780, setting = 2, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case1/070415ad_')

c0.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 1440, filtername = 'WCB', thinning = 20, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/020415a_')
c0.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 2880, filtername = 'WCB', thinning = 20, setting = 1, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/020415b_')
c0.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 4320, filtername = 'WCB', thinning = 20, setting = 1, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/020415c_')

#c0.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 1440, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/230215j_')
#c0.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 2880, setting = 2, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/230215k_')
#c0.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 4320, setting = 2, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case0/230215l_')

c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 720, filtername = 'WCB_Cy1_new', thinning = 20, setting = 1, cbar = False,
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/020415d_')
c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 2160, filtername = 'WCB_Cy1_new', thinning = 20, setting = 1, cbar = False,
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/020415e_')
c2.draw_trj_dot([['var4', 310, 'THETA'], 'PMSL'], 3600, filtername = 'WCB_Cy1_new', thinning = 20, setting = 1, cbar = True,
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/020415f_')

#c2.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 720, setting = 5, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/230215p_')
#c2.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 2160, setting = 5, cbar = False, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/230215q_')
#c2.draw_contour(['PMSL', 'var145_S', "TOT_PREC_S"], 3600, setting = 5, cbar = True, savebase = '/usr/users/stephan.rasp/Dropbox/figures/thesis/Case2/230215r_')

#170415
c0.draw_slice_hist('P_dt', 0.1, 'Wconv2', onlyasc = 'P600', mult = -1, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/170415a_')
c0.draw_slice_hist('P_dt', 0.1, 'Wnonconv2', onlyasc = 'P600', mult = -1, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/170415b_')
c1.draw_slice_hist('P_dt', 0.1, 'WCB_Cy1', onlyasc = 'P600', mult = -1, savebase ='/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/170415c_')
c2.draw_slice_hist('P_dt', 0.1, 'WCB_Cy1_new', onlyasc = 'P600', mult = -1, savebase ='/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/170415d_')

#200415

c0.draw_vs_p('T', ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (200, 320), ylabel = 'T [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/200415a_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 3, ax2upper = None, letter = '(b)')    
    
c0.draw_vs_p('T', ['Wconv2', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (200, 320), ylabel = 'T [K]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/200415b_', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 3, ax2upper = None, letter = '(a)') 


c0.draw_vs_p('THETA_dt', ['Wnonconv2', 'WCB_Cy1_new'], 'P600', (200, 1000), idtext = '', ylim = (-5, 15), ylabel = r'DHR [K h$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_p/200415c_', extobj = c2, legnames = ['OCTnc', 'JAN'], legpos = 4, ax2upper = None, letter = '(d)', mult = 3600)    
    
c0.draw_vs_p('THETA_dt', ['Wconv2', 'WCB_Cy1'], 'P600', (200, 1000), idtext = '', ylim = (-5, 75), ylabel = r'DHR [K h$^{-1}$]', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_p/200415d_', extobj = c1, legnames = ['OCTc', 'JUL'], legpos = 2, ax2upper = None, letter = '(c)', mult = 3600)     


c0.draw_centered_vs_t('latitude', ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (30, 75), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/200415f_', legnames = ['OCTnc', 'JAN'], letter = '(b)', ax2 = None, ylabel = 'latitude [deg]', mult = 1, legpos = 2)

c0.draw_centered_vs_t('latitude', ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (30, 75), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/200415e_', legnames = ['OCTc', 'JUL'], letter = '(a)', ax2 = None, ylabel = 'latitude [deg]', mult = 1., legpos = 2)

c0.draw_centered_vs_t('longitude', ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-30, 50), xlim = (-36, 36), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/200415h_', legnames = None, letter = '(d)', ax2 = None, ylabel = 'longitude [deg]', mult = 1)

c0.draw_centered_vs_t('longitude', ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-30, 50), xlim = (-12, 12), select=(-6, 6), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/200415g_', legnames = None, letter = '(c)', ax2 = None, ylabel = 'longitude [deg]', mult = 1)

#220415
c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wnonconv2', 'WCB_Cy1_new'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-6, 6), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/220415a_', legnames = ['OCTnc', 'JAN'], letter = '(b)', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 4)

c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['Wconv2', 'WCB_Cy1'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-6, 6), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/220415b_', legnames = ['OCTc', 'JUL'], letter = '(a)', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 4)

c2.draw_hist_2d('P', 'POT_VORTIC', filtername = 'WCB_Cy1_new', after = 'Pcross600', plus = 2160, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/220415c_')

c1.draw_hist_2d('P', 'POT_VORTIC', filtername = 'WCB_Cy1', after = 'Pcross600', plus = 2160, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/220415d_')

#c1.draw_centered_vs_t('POT_VORTIC', ['5pvu36', '5pvu36'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-2, 10), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/220415e_', legnames = ['JUL', 'JAN'], letter = '(a)', ax2 = None, ylabel = 'PV [pvu]', mult = 1)

#c1.draw_centered_vs_t('P', ['5pvu36', '5pvu36'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (100, 1000), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/220415f_', legnames = None, letter = '(b)', ax2 = None, ylabel = 'p [hPa]', mult = 1)

#c1.draw_centered_vs_t('THETA', ['5pvu36', '5pvu36'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (280, 350), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/220415f_', legnames = None, letter = '(c)', ax2 = None, ylabel = r'$\theta$ [K]', mult = 1)

c0.draw_hist_2d('P', 'var4', filtername = 'Wconv2', after = 'Pcross600', plus = 2160, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/220415g_')

c0.draw_hist_2d('P', 'var4', filtername = 'Wnonconv2', after = 'Pcross600', plus = 2160, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/220415h_')

c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['5pvu36nonconv2', '5pvu36'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (-2, 10), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/220415i_', legnames = ['OCTnc', 'JAN'], letter = '(b)', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 4)

c0.draw_centered_vs_t(['var4', 'POT_VORTIC'], ['5pvu36conv2', '5pvu36'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (-2, 10), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/220415j_', legnames = ['OCTc', 'JUL'], letter = '(a)', ax2 = None, ylabel = 'PV [pvu]', mult = 1, legpos = 4)

c0.draw_centered_vs_t('P', ['5pvu36nonconv2', '5pvu36'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (100, 1000), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/220415k_', legnames = None, letter = '(d)', ax2 = None, ylabel = 'p [hPa]', mult = 1, legpos = 4)

c0.draw_centered_vs_t('P', ['5pvu36conv2', '5pvu36'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (100, 1000), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/220415l_', legnames = None, letter = '(c)', ax2 = None, ylabel = 'p [hPa]', mult = 1, legpos = 4)

c0.draw_centered_vs_t('THETA', ['5pvu36nonconv2', '5pvu36'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (280, 350), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/220415m_', legnames = None, letter = '(f)', ax2 = None, ylabel = r'$\theta$ [K]', mult = 1, legpos = 4)

c0.draw_centered_vs_t('THETA', ['5pvu36conv2', '5pvu36'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (280, 350), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/220415n_', legnames = None, letter = '(e)', ax2 = None, ylabel = r'$\theta$ [K]', mult = 1, legpos = 4)

c0.draw_centered_vs_t('QV', ['5pvu36nonconv2', '5pvu36'], 'Pcross600', extobj = c2, plottype = 'Fill', ylim = (0, 15), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case2_20090129/vs_t/220415o_', legnames = None, letter = '(h)', ax2 = None, ylabel = r'$q$ [g kg$^{-1}$]', mult = 1000, legpos = 4)

c0.draw_centered_vs_t('QV', ['5pvu36conv2', '5pvu36'], 'Pcross600', extobj = c1, plottype = 'Fill', ylim = (0, 15), xlim = (-12, 72), select=(-6, 36), savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/vs_t/220415p_', legnames = None, letter = '(g)', ax2 = None, ylabel = r'$q$ [g kg$^{-1}$]', mult = 1000, legpos = 4)





