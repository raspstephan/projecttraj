import traj_tools as trj
c0 = trj.loadme('c0_rot.trj')

c0.draw_contour(["PMSL", "TOT_PREC_S", 'var145_S'],1440, savebase = '/home/users/stephan.rasp/Dropbox/figures/thesis/Case0/', idtext = '081014a')
c0.draw_contour(["PMSL", "TOT_PREC_S", 'var145_S'],2880, savebase = '/home/users/stephan.rasp/Dropbox/figures/thesis/Case0/', idtext = '081014b')
c0.draw_contour(["PMSL", "TOT_PREC_S", 'var145_S'],4320, savebase = '/home/users/stephan.rasp/Dropbox/figures/thesis/Case0/', idtext = '081014c')

c1 = trj.loadme('c1.trj')
c1.draw_contour(["PMSL", "TOT_PREC_S", 'var145_S'],1080, savebase = '/home/users/stephan.rasp/Dropbox/figures/thesis/Case1/', idtext = '081014d')
c1.draw_contour(["PMSL", "TOT_PREC_S", 'var145_S'],2520, savebase = '/home/users/stephan.rasp/Dropbox/figures/thesis/Case1/', idtext = '081014e')
c1.draw_contour(["PMSL", "TOT_PREC_S", 'var145_S'],3960, savebase = '/home/users/stephan.rasp/Dropbox/figures/thesis/Case1/', idtext = '081014f')
