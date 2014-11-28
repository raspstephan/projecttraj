import traj_tools as trj

c0 = trj.loadme('c0_rot.trj')

c0.calc_cape(200 * 1000., getp = True)