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








