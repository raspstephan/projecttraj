"""
Helper functions for traj tools
submodule. 
Do not call functions explicitely
"""

import cPickle
import glob
import struct
import numpy as np


####################
### Helper Functions
####################

def CreateFileArray(FileDir):
    """
    Creates list of files in a directory
    """
    files = FileDir + "*"
    allfiles = glob.glob(files)
    allfiles = sorted(allfiles)
    print "Read", len(allfiles), "files"

    return allfiles


def LoadCaseSpec(CaseName="test"):
    """
    Loads Specifics for certain case
    """
    # global pfiles, rfiles, cfile, TrjOffset, WCBIndM, FileList, rvar, pvar, DefFile, dt
    if CaseName == "test":
        print "Loading Specifications for TestCase"
        SourceDir = "/home/cosmo/tobias.selz/cosmo_data/caseWCB/d4deout/"
        COSMOName = "lfff"
        dt = 5.   # Minutes
        rsuff = "_5m"
        psuff = "p_5m"
        rfiles = glob.glob(SourceDir + "*" + rsuff)
        rfiles = sorted(rfiles)
        rvar = ['CLCT_S', 'PMSL', 'TOT_PREC_S', 'var145_S', 'var146_S']
        pfiles = glob.glob(SourceDir + "*" + psuff)
        pfiles = sorted(pfiles)
        pvar = ['OMEGA', 'U', 'V', 'T', 'FI']
        rfiles = [x for x in rfiles if x not in pfiles]
        cfile = SourceDir + "lfff00000000c_1h"
        DefFile = "/home/cosmo/tobias.selz/cosmo_data/caseWCB/d4deout/ltradef"
        TrjStart = LoadTrjDef(DefFile)[0, 4]
        TrjOffset = int(TrjStart / dt)
        f = open("/home/scratch/users/stephan.rasp/traj_data/IndexMatrix", "r")
        WCBIndM = cPickle.load(f)
        f.close()
        FileList = CreateFileArray("/home/scratch/users/stephan.rasp/traj_data/test_case_")
        return pfiles, rfiles, cfile, TrjOffset, WCBIndM, FileList, rvar, pvar, DefFile, dt
    else:
        return None
    
    
def LoadTrjDef(Deffile):
    " Returns Matrix [Trajectory nr starting with 0][nr, starting x, y, z, t]"
    f = open(Deffile, "r")
    buf = f.read(4)
    ntra = struct.unpack('>i', buf)[0]
    buf = f.read(4)
    ntrace = struct.unpack('>i', buf)[0]
    DefM = np.zeros((ntra, 5), float)
    lrec = 20
    for i in xrange(ntra):
        pos = 8 + i * lrec
        f.seek(pos, 0)
        buf = f.read(4)
        DefM[i, 0] = struct.unpack('>i', buf)[0]
        buf = f.read(4)
        DefM[i, 1] = struct.unpack('>f', buf)[0]
        buf = f.read(4)
        DefM[i, 2] = struct.unpack('>f', buf)[0]
        buf = f.read(4)
        DefM[i, 3] = struct.unpack('>f', buf)[0]
        buf = f.read(4)
        DefM[i, 4] = struct.unpack('>i', buf)[0]
    f.close()
    return DefM
