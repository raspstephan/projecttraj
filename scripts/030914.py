# Script for 030914
import traj_tools as trj


c0 = trj.loadme('c0.trj')
c1 = trj.loadme('c1.trj')

c0.new_asc(100)
c0.new_asc(400)
c0.new_asc(600)

c0.create_filter('WCB', [('P600', 0, 2880)])
c0.create_filter('WCB_NonConv', [('P600', 0, 2880), ('P400', 120, 2880)])
c0.create_filter('WCB_Conv', [('P600', 0, 2880), ('P400', 0, 120)])

# Case 0, P600vP100 Scatter
c0.draw_scatter('P600', 'P100', factor2 = 6, carray = 'P100', 
                filtername = 'WCB_NonConv', 
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
                idtext = '030914a')
c0.draw_scatter('P600', 'P100', factor2 = 6, carray = 'P100', 
                filtername = 'WCB_Conv', 
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
                idtext = '030914b')
c0.draw_scatter('P600', 'P100', factor2 = 6, carray = 'P100', 
                filtername = 'WCB', 
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
                idtext = '030914c')



c1.new_asc(100)
c1.new_asc(400)
c1.new_asc(600)

c1.create_filter('WCB', [('P600', 0, 2880)])
c1.create_filter('WCB_NonConv', [('P600', 0, 2880), ('P400', 120, 2880)])
c1.create_filter('WCB_Conv', [('P600', 0, 2880), ('P400', 0, 120)])

# Case 1, P600vP100 Scatter
c1.draw_scatter('P600', 'P100', factor2 = 6, carray = 'P100', 
                filtername = 'WCB_NonConv', 
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
                idtext = '030914d')
c1.draw_scatter('P600', 'P100', factor2 = 6, carray = 'P100', 
                filtername = 'WCB_Conv', 
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
                idtext = '030914e')
c1.draw_scatter('P600', 'P100', factor2 = 6, carray = 'P100', 
                filtername = 'WCB', 
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
                idtext = '030914f')