import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['font.size'] = 8

jan = [997, 995, 992, 989, 987, 985, 983, 980, 979, 977, 974, 970, 966, 961, 
       956, 956, 954, 953, 953, 953, 953, 953, 953, 953]

oct = [1004, 1003, 1002, 1001, 1000, 999, 998, 997, 996, 995, 995, 994, 994, 
       994, 994, 994, 994, 994, 994, 993, 993, 993, 993, 993, 993, 994, 994, 
       994, 994, 994, 993, 993, 993, 992, 992, 992, 991, 991, 991, 990, 990, 
       990, 990, 990]

jul = [1011, 1011, 1010, 1010, 1008, 1007, 1007, 1006, 1006, 1006, 1005, 1004, 1004, 1004, 1003 ,1003, 1002, 
       1002, 1002, 1002, 1002, 1003, 1003, 1003, 1003, 1003, 1004,1004,1004,1004,1004,1004,1004,1004,1004,1004,
       1003,1003,1003,1003,1002, 1003,1002, 1002, 1001, 1001, 1001, 1000, 1000, 1000, 999]
fig = plt.figure(figsize = (3.15, 3.15))


plt.plot(range(len(jul)), jul, c = 'black', linestyle = ':', label = 'JUL',
         linewidth = 2)
plt.plot(range(len(oct)), oct, c = 'black', linestyle = '--', label = 'OCT',
         linewidth = 2)
plt.plot(range(len(jan)), jan, c = 'black', linestyle = '-', label = 'JAN', 
         linewidth = 2)
plt.legend(loc = 4, fontsize = 8)
plt.grid()

plt.xlabel('Time [h]')
plt.ylabel('Minimum surface pressure [hPa]')

plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/thesis/deepening_rates',
            bbox_inches = 'tight', dpi = 300)