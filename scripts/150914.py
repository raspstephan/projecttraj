# Script for 15.09.2014

c0 = trj.loadme('c0.trj')
c0.new_allasct(100, 30)
c0.new_allasct(100, 10)

for t in range(720, 7200, 180):
     c0.draw_asc_loc('P100in30', ["PMSL", "TOT_PREC_S"], 'WCB',t, 90, 
                     idtext = '150914a', 
                     savebase = 
                     '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/P100in30/')
     
c0.create_filter('WCB_NonConv_1620', [('P600', 0, 2880), ('P400', 120, 2880), ('startt', 1620)])

c0.draw_trj_all([], 'WCB_NonConv_1620', '/usr/users/stephan.rasp/tmp/w_test', carray = 'z', centdiff = True, sigma=5)
