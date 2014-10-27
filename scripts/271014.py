import traj_tools as trj
c0 = trj.loadme('c0_rot.trj')
for level in [700, 500, 300, 100]:
    array = c0.draw_intersect_hor(level, idtext='271014a', filtername='WCB', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/intersect/')
    trj.plots.draw_hist(array, idtext='271014b', xlabel = ('Vertical velocity at hPa: ' + str(level)), savename = ('/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/intersect/' + 'hist_' + str(level)) + '_271014b', range=(-15,15))
    
c1 = trj.loadme('c1.trj')
for level in [700, 500, 300, 100]:
    array = c1.draw_intersect_hor(level, idtext='271014c', filtername='WCB', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720//intersect/')
    trj.plots.draw_hist(array, idtext='271014d', xlabel = ('Vertical velocity at hPa: ' + str(level)), savename = ('/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720//intersect/' + 'hist_' + str(level)) + '_271014d', range=(-15,15))