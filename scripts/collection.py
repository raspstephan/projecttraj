import traj_tools as trj

c0 = trj.loadme('c0_rot.trj')
c0.create_filter('WCB_Conv', [('P600', 0, 2880), ('P400', 0, 120)])
c1 = trj.loadme('c1.trj')


## 061114
#c0.draw_centered_vs_t('P', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('P', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')


#c0.draw_centered_vs_t('var4', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (-3, 3), sigma = 1, idtext = '061114c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
#c0.draw_centered_vs_t('var4', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (-3, 3), sigma = 1, idtext = '061114d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

#c1.draw_centered_vs_t('P', 'WCB', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
#c1.draw_centered_vs_t('var4', 'WCB', 'P600', plottype='Smooth', ylim = (-3, 3), sigma = 1, idtext = '061114f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')


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
c1.draw_centered_vs_t('THETAE', 'WCB', 'P600', plottype='Smooth', ylim = (280, 330), sigma = 1, idtext = '111114i', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/') 