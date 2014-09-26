indir = '/home/scratch/users/stephan.rasp/traj_data/'
outdir = '/usr/users/stephan.rasp/tmp/traj_files/'

import traj_tools as trj

trj.utils.convert_pickle2netcdf(indir, outdir)

c0 = trj.loadme('c0_rot.trj')

#c0.draw_trj_evo(["PMSL", "TOT_PREC_S"], 'WCB', 360, 180, '260914a', 
                #'/usr/users/stephan.rasp/tmp/c0evo/')
                
c0.new_asc(400)
c0.create_filter('WCB_NonConv', [('P600', 0, 2880), ('P400', 120, 2880)])


c0.draw_trj_evo(["PMSL", "TOT_PREC_S"], 'WCB_NonConv', 360, 180, '260914b', 
                '/usr/users/stephan.rasp/tmp/c0evo2/')
