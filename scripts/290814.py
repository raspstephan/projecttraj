import traj_tools as trj
c0 = trj.loadme('c0.save')
# c0.create_filter('WCB_nonconv1', [('P600', 0, 2880), ('P400', 120, 2880)])
c0.draw_hist('P600', 'WCB_nonconv1', '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '290814a')
c0.draw_trj_all([], 'WCB_nonconv1_1080', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '290814b')
c0.draw_trj_all([], 'WCB_nonconv1_1620', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', idtext = '290814c')
c0.draw_trj_all([], 'WCB_nonconv1_1080', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/onlyasc_', idtext = '290814d', onlyasc='P100', linewidth = 1.5)
c0.draw_trj_all([], 'WCB_nonconv1_1620', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/onlyasc_', idtext = '290814e', onlyasc='P100', linewidth = 1.5)
c0.draw_scatter('P600', 'P100', factor2 = 6, filtername='WCB_nonconv1', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/scatter_P600_P100_WCB_nonconv1.png', idtext = '290814f')
c0.draw_scatter('P600', 'P100', factor2 = 6, filtername='WCB_nonconv1_1080', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/scatter_P600_P100_WCB_nonconv1_1080.png', idtext = '290814g')
c0.draw_scatter('P600', 'P100', factor2 = 6, filtername='WCB_nonconv1_1620', savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/scatter_P600_P100_WCB_nonconv1_1620.png', idtext = '290814h')
