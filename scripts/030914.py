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

c0.draw_scatter('P600', 'P100', factor2 = 6, carray = 'P100', 
                filtername = 'WCB_NonConv', 
                savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
                idtext = '030914a')



#c1.new_asc(100)
#c1.new_asc(400)
#c1.new_asc(600)