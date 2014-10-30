import traj_tools as trj
c0 = trj.loadme('c0_rot.trj')
for level in [700, 500, 300]:
    array = c0.draw_intersect_hor(level, idtext='301014a', filtername='WCB', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/intersect/')
c0.create_filter('WCB_Conv', [('P600', 0, 2880), ('P400', 0, 120)])    
c1 = trj.loadme('c1.trj')
for level in [700, 500, 300]:
    array = c1.draw_intersect_hor(level, idtext='301014b', filtername='WCB', savebase='/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720//intersect/')

levellist = [300, 400, 500, 600, 700, 800, 900]
vellist1 = []
vellist2 = []
vellist3 = []
vellist4 = []
for level in levellist:
    vellist1.append((c0.draw_intersect_hor(level, filtername='WCB')).mean())
    vellist2.append((c1.draw_intersect_hor(level, filtername='WCB')).mean())
    vellist3.append((c0.draw_intersect_hor(level, filtername='WCB_NonConv')).mean())
    vellist4.append((c0.draw_intersect_hor(level, filtername='WCB_Conv')).mean())
    plt.close('all')

plt.plot(vellist1, levellist, 'r', label = 'Case0 P600')
plt.plot(vellist2, levellist, 'g', label = 'Case1 P600')
plt.plot(vellist3, levellist, 'b', label = 'Case0 P600 NonConvective')
plt.plot(vellist4, levellist, 'y', label = 'Case0 P600 Convective')
plt.gca().invert_yaxis()
plt.legend(loc=4)

plt.text(0.94, 1.02, '301014c', transform = plt.gca().transAxes, 
             fontsize = 6)
plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/vel_vs_p_301014c.png')