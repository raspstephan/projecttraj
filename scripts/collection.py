import traj_tools as trj

c0 = trj.loadme('c0_rot.trj')
c0.create_filter('WCB_Conv', [('P600', 0, 2880), ('P400', 0, 120)])
c1 = trj.loadme('c1.trj')

c0.draw_centered_vs_t('P', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114a', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('P', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114b', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')


c0.draw_centered_vs_t('var4', 'WCB_NonConv', 'P600', plottype='Smooth', ylim = (-3, 3), sigma = 1, idtext = '061114c', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')
c0.draw_centered_vs_t('var4', 'WCB_Conv', 'P600', plottype='Smooth', ylim = (-3, 3), sigma = 1, idtext = '061114d', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/')

c1.draw_centered_vs_t('P', 'WCB', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c1.draw_centered_vs_t('var4', 'WCB', 'P600', plottype='Smooth', ylim = (-3, 3), sigma = 1, idtext = '061114f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')



