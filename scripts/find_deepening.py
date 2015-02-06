import matplotlib.pyplot as plt
import numpy as np

jan = [997, 995, 992, 989, 987, 985, 983, 980, 979, 977, 974, 970, 966, 961, 
       956, 956, 954, 953, 953, 953, 953, 953, 953, 953]

oct = [1004, 1003, 1002, 1001, 1000, 999, 998, 997, 996, 995, 995, 994, 994, 
       994, 994, 994, 994, 994, 994, 993, 993, 993, 993, 993, 993, 994, 994, 
       994, 994, 994, 993, 993, 993, 992, 992, 992, 991, 991, 991, 990, 990, 
       990, 990, 990]
fig = plt.figure()

plt.plot(range(len(jan)), jan, c = 'black', linestyle = '-', label = 'JAN', 
         linewidth = 2)
plt.plot(range(len(oct)), oct, c = 'black', linestyle = '--', label = 'OCT',
         linewidth = 2)
plt.legend()
plt.grid()

plt.xlabel('Time [h]')
plt.ylabel('Minimum surface pressure [hPa]')

plt.savefig('/usr/users/stephan.rasp/Dropbox/figures/thesis/deepening_rates',
            bbox_inches = 'tight')