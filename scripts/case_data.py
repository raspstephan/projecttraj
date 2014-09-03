# This file contains info on how to initialize the case object files
# Name = trj.TrjObj(datadir, pollon, pollat)
import traj_tools as trj

# Case 0 
c0 = trj.TrjObj('/home/scratch/users/stephan.rasp/Case0_20121013/d4deout/', 
                -165, 35)

# Case 1
c1 = trj.TrjObj('/home/scratch/users/stephan.rasp/Case1_20070720/d4deout/', 
                -170, 40)


# Save objects
c0.saveme('c0.trj')
c1.saveme('c1.trj')
