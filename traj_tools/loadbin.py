"""
bin submodule
Submodule of traj_tools
Opening and Saving binary trajectory files
"""
from common import *
import cPickle
import numpy as np
import struct
import glob
import time



#############################################
### Loading and Saving Data from Binary Files
#############################################


def LoadTrjSpecif(Deffile):
    " Returns array [number trajectories, number tracers]"
    f = open(Deffile, "r")
    buf = f.read(8)
    ntra = struct.unpack('>i', buf[:4])[0]
    ntrace = struct.unpack('>i', buf[4:])[0]
    f.close()
    return [ntra, ntrace]


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


def SaveBin2File(Files, MSize, SaveSize, SaveBase, StartDate, dt):
    """
    Saves all trajectories in Matrix Files with Specification file
    For StartDate and dt use datetime module!
    Attention: MSize must be a multiple of SaveSize!!!
    """

    # Check for correct input
    if MSize % SaveSize != 0:
        raise ValueError('MSize is not a multiple of SaveSize!')

    # Read ntra, ntraces
    [ntra, ntrace] = LoadTrjSpecif(Files[0])

    # Save files with cPickle
    irange = int(np.ceil(ntra / MSize)) + 1
    jrange = int(MSize / SaveSize)
    print irange, jrange
    for i in range(irange):    # Iterate over read matrices
        M = ReadBin2Matrix(Files, i*MSize, MSize)    # Read to Matrix

        for j in range(jrange):    # Iterate over save matrices
            print "Saving file:", str(i*jrange + j).zfill(4)
            savefile = open(SaveBase + str(i*jrange + j).zfill(4), "w")

            # First, save Specification file
            cPickle.dump([ntra, ntrace, i*jrange + j, i*MSize + j*SaveSize, i*MSize + (j+1)*SaveSize - 1,
                          len(Files), StartDate, dt], savefile, 2)
            # Second, save Matrix
            cPickle.dump(M[(j*SaveSize):((j+1)*SaveSize)], savefile, 2)
            savefile.close()
        del M    # delete Matrix to save RAM
    pass


def ReadBin2Matrix(Files, NrStart, NrMax):
    """
    Function for reading binary trajectory files to matrix.
    Files: list of binary file locations
    NrStart: lowest wanted trj_nr
    NrMax: number of trj's
    """

    [ntra, ntrace] = LoadTrjSpecif(Files[0])    # Read ntra, ntraces
    M = np.zeros((NrMax, ntrace + 3, len(Files)))    # Create Matrix [trj_nr, tracer, time]
    told = time.time()    # Initialize time counter

    for t in range(len(Files)):    # Iterate over time, third dim in M
        if t % 50 == 0:    # Time module
            clock = time.time() - told
            told = time.time()
            print "Step =", t, clock

        lrec = 4 + (3+ntrace) *4    # Length of one trj
        i = 1    # Index for going through binary data
        f = open(Files[t])    # Open file
        buf = f.read()    # Read entire string

        while(True):
            if(buf[((i-1)*lrec)+8 : ((i-1)*lrec)+12] == ''):
                # Break for end of string
                x = nan
                y = nan
                z = nan
                break

            # Read trj nr = ind
            ind = int(struct.unpack('>i', buf[((i-1)*lrec)+8 : ((i-1)*lrec)+12])[0]) - 1

            # If ind is in required range [NrStart,NrMax+NrStart[, write content to Matrix
            if ind < NrMax + NrStart and ind >= NrStart:
                for j in range(ntrace + 3):
                    M[ind-NrStart, j, t] = struct.unpack('>f', buf[((i-1)*lrec + (j+1)*4) + 8 :
                                                         (i-1)*lrec + (j+1)*4 + 12])[0] 
            i += 1    # Forward index
        f.close()
    return M

#End of submodule
