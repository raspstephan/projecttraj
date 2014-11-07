import traj_tools as trj
c1 = trj.loadme('c1.trj')

c1.draw_centered_vs_t('P', 'WCB', 'P600', plottype='Smooth', ylim = (100, 1000), sigma = 1, idtext = '061114e', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')
c1.draw_centered_vs_t('var4', 'WCB', 'P600', plottype='Smooth', ylim = (-3, 3), sigma = 1, idtext = '061114f', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/')

