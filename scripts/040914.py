# Script 040914

import traj_tools as trj


c0 = trj.loadme('c0.trj')
c1 = trj.loadme('c1.trj')

c0.draw_contour(["PMSL", "TOT_PREC_S", 'var145_S'], 0, interval = 240, 
                savebase = '/usr/users/stephan.rasp/tmp/evo_test/', 
                idtext = '040914a')
c1.draw_contour(["PMSL", "TOT_PREC_S", 'var145_S'], 0, interval = 240, 
                savebase = '/usr/users/stephan.rasp/tmp/evo_test2/', 
                idtext = '040914b')



c0.draw_hist(['P600', 'P100', 1, 6], 'WCB_NonConv', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             log = True, idtext = '040914c')
c0.draw_hist(['P600', 'P100', 1, 6], 'WCB_Conv', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             log = True, idtext = '040914d')

c1.draw_hist(['P600', 'P100', 1, 6], 'WCB', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
             log = True, idtext = '040914e')