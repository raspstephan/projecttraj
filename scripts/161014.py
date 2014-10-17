import traj_tools as trj
c0 = trj.loadme('c0_rot.trj')

for i in range(720, 7200+360, 360):
    c0.draw_hist_2d(["TOT_PREC_S"], i, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist2d/all_161014a_', idtext='161014a')
c1 = trj.TrjObj('/home/scratch/users/stephan.rasp/Case1_20070720/d4deout/')
for i in range(720, 6120 + 360, 360):
    c1.draw_hist_2d(["TOT_PREC_S"], i, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720//hist2d/all_161014b_', idtext='161014b')
    
c0.new_delta('P', mode = 'mintomax_r')
c0.create_filter('plus300', [('deltaPmintomax_r', -1000, -300)])
for i in range(1440, 5760 + 1440, 1440):
    c0.draw_hist_2d(["TOT_PREC_S"], i, savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/hist2d/P300_161014c_', idtext='161014c', filtername='plus300')

