import traj_tools as trj
c0 = trj.loadme('c0_rot.trj')
c1 = trj.TrjObj('/home/scratch/users/stephan.rasp/Case1_20070720/d4deout/')

c0.new_delta('QV', mode = 'mintomax_r')
c1.new_delta('QV', mode = 'mintomax_r')

c0.draw_hist('deltaQVmintomax_r', savebase = 
             '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             idtext = '221014a', range = (-0.025, 0))
c1.draw_hist('deltaQVmintomax_r', savebase = 
             '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
             idtext = '221014b', range = (-0.025, 0))

c0.new_max_diff('z')
c1.new_max_diff('z')

c0.draw_hist('z_max_diff', savebase = 
             '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             idtext = '221014c', range = (0, 15))
c1.draw_hist('z_max_diff', savebase = 
             '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
             idtext = '221014d', range = (0, 15))

c0.new_delta('THETA', mode = 'mintomax')
c1.new_delta('THETA', mode = 'mintomax')
c0.create_filter('deltaTHETA5', [('deltaTHETAmintomax', 5, 1000)])
c1.create_filter('deltaTHETA5', [('deltaTHETAmintomax', 5, 1000)])

c0.draw_hist('z_max_diff', savebase = 
             '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             idtext= '221014e', filtername='deltaTHETA5', range = (0, 15))
c1.draw_hist('z_max_diff', savebase = 
             '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
             idtext= '221014f', filtername='deltaTHETA5', range = (0, 15))



