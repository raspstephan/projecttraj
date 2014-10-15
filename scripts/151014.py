

import traj_tools as trj
c0 = trj.loadme('c0_rot.trj')
c1 = trj.TrjObj('/home/scratch/users/stephan.rasp/Case1_20070720/d4deout/')

c0.new_delta('THETA', mode = 'climb')
c1.new_delta('THETA', mode = 'climb')


c0.draw_hist('deltaTHETAclimb', idtext='151014a', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             range = (0, 100))
c1.draw_hist('deltaTHETAclimb', idtext='151014b', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/',
             range = (0, 100))

c0.new_delta('P', mode = 'climb_r')
c1.new_delta('P', mode = 'climb_r')


c0.draw_hist('deltaPclimb_r', idtext='151014c', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             range = (-2500, 0))
c1.draw_hist('deltaPclimb_r', idtext='151014d', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
             range = (-2500, 0))

c0.new_delta('QV', mode = 'climb_r')
c1.new_delta('QV', mode = 'climb_r')

c0.draw_hist('deltaQVclimb_r', idtext='151014e', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             range = (-0.035, 0))
c1.draw_hist('deltaQVclimb_r', idtext='151014f', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
             range = (-0.035, 0))

c0.new_delta('P', mode = 'mintomax_r')
c1.new_delta('P', mode = 'mintomax_r')

c0.draw_hist('deltaPmintomax_r', idtext='151014g', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case0_20121013/', 
             range = (-1000, 0))
c1.draw_hist('deltaPmintomax_r', idtext='151014h', 
             savebase = '/usr/users/stephan.rasp/Dropbox/figures/Case1_20070720/', 
             range = (-1000, 0))