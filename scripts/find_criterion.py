import traj_tools as trj
import matplotlib.pyplot as plt

c1 = trj.loadme('c1.trj')
c2 = trj.loadme('c2.trj')

plist = ['100', '200', '300', '400', '500']
hlist = [0.5, 1, 2, 3, 4, 5, 10]

totc1 = c1.count_trjs('WCB')
totc2 = c2.count_trjs('WCB_Cy1_new')

alist = []
blist = []
clist = []

for p in plist:
    print 'Testing:', p
    for h in hlist:
        print 'Maxtime:', h
        c1.create_filter('Test'+ p + str(h), [('P600', 0, 2880), 
                                    ('P' + p +'withinP600', 0, h * 60)])
        c2.create_filter('Test'+ p + str(h), [('P600', 0, 2880), ('-10201025', True), 
                                    ('P' + p +'withinP600', 0, h * 60)])
        partc1 = c1.count_trjs('Test'+p+ str(h))
        partc2 = c2.count_trjs('Test'+p+ str(h))
        
        alist.append(partc1 / float(totc1))
        blist.append(1 - partc2 / float(totc2))
        clist.append(partc1 / float(totc1) - partc2 / float(totc2) + 1)

h4list = hlist * 4

plt.plot(alist, marker = '*')
plt.plot(blist, marker = '*')
plt.plot(clist, marker = '*')
plt.show()