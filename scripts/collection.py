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
