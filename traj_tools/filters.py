"""
filters module!!
Submodule of traj_tools
Opening and Filtering saved trajectory files
Naming conventions:
Files = all trajectories
Matrix = one trajectory file (i.e. 1000 trajectories)
Array = one trajectory
"""
from common import *
import cPickle
import numpy as np
import loadbin




def OpenSaveFile(FileName):
    """
    Unpacks Saved Trajectory files
    Returns tuple (Specification file, Matrix)
    """
    f = open(FileName, "r")
    Spec = cPickle.load(f)
    Matrix = cPickle.load(f)
    f.close()
    return (Spec, Matrix)


def FilterFiles(FileList, TraceInd, Criterion, LenMax, LenMin = 0):
    """
    Filter goes through all files, returns matrix
    """
    IndMatrix = []
    for i in range(len(FileList)):
        print "Completed", i, "out of", len(FileList)
        Matrix = OpenSaveFile(FileList[i])[1]
        IndMatrix.append(FilterMatrix(Matrix, TraceInd, Criterion, LenMax, LenMin))
        del Matrix    # Delete Matrix to save RAM
    return IndMatrix


def FilterMatrix(Matrix, TraceInd, Criterion, LenMax, LenMin = 0):
    """
    Filters one 3D matrix by trajectories fulfilling Criterion 
    in Len Interval given by LenMax and LenMin for given TraceInd
    """
    IndList = []
    for i in range(len(Matrix[:, 0, 0])):
        if np.average(Matrix[i, TraceInd, :]) == 0:    # Checks for all Zero arrays"
            print "Given array is all zeros, break!"
            break
        elif MinXSpan(Matrix[i, TraceInd, :], Criterion) <= LenMax and MinXSpan(Matrix[i, TraceInd, :], Criterion) >= LenMin:
            IndList.append(i)
    return IndList


def IndData(FileList, IndList, MaxFile="default"):
    """
    Gives Matrix of data for all indexed trajectories
    """
    f = open(FileList[0], "r")
    Spec = cPickle.load(f)
    f.close()
    NrTrj = 0
    if MaxFile == "default":
        MaxFile = min(len(FileList), len(IndList))
    for i in range(MaxFile):
        NrTrj += len(IndList[i])
    FilterM = np.zeros((NrTrj, Spec[1]+3, Spec[5])) 
    ind = 0
    for i in range(MaxFile):
        print "Opening file", i
        RawM = OpenSaveFile(FileList[i])[1]
        for j in IndList[i]:
            FilterM[ind, :, :] = RawM[j, :, :]
            ind += 1    
    return FilterM


def StartPos(Matrix2D):
    """
    Returns Starting index, x,y,z,p of trajectory
    """
    StartInd = np.where(Matrix2D[7] > 7)[0][0]
    return (StartInd, Matrix2D[0, StartInd], Matrix2D[1, StartInd], Matrix2D[2, StartInd], Matrix2D[7, StartInd])


def EndPos(Matrix2D):
    """
    Returns End Index, x, y, z, p of trajcetory
    """
    EndInd = np.where(Matrix2D[7] > 7)[0][-1]
    return (EndInd, Matrix2D[0, EndInd], Matrix2D[1, EndInd], Matrix2D[2, EndInd], Matrix2D[7, EndInd])

def StartIndList(CaseName = "test"):
    """
    Returns list with all occuring starting positions
    """
    pfiles, rfiles, cfile, TrjOffset, WCBIndM, FileList, rvar, pvar, DefFile, dt = LoadCaseSpec()
    StartIndList = list((np.unique(loadbin.LoadTrjDef(DefFile)[:,4]).astype(int)-np.unique(loadbin.LoadTrjDef(DefFile)[:,4]).astype(int)[0]) / 5)
    return StartIndList


def MinXMatrix(FileList, TraceInd, Criterion, StartInd = False, IndMatrix = False, Flat = False):
    """
    Returns 2D list of fastest ascent times for given criterion
    if Index Matrix is given, search only for indexed trajectories
    if Flat == True, returns one long array. To be used for histograms!
    Note: StartInd not compatible with Flat == False!!!
    """
    AscMatrix = []
    if IndMatrix == False:
        for i in range(len(FileList)):
            if i % 10 == 0:
                print "Going through file %i of %i" % (i,len(FileList))
            if StartInd:
                AscArray = [[], []]
            else:
                AscArray = []
            M = OpenSaveFile(FileList[i])[1]
            for j in range(M.shape[0]):
                if StartInd:
                    AscArray[1].append(StartPos(M[j])[0])
                    AscArray[0].append(MinXSpan(M[j, TraceInd, :], Criterion))
                else:
                    AscArray.append(MinXSpan(M[j, TraceInd, :], Criterion))
            AscMatrix.append(AscArray)
                                
    elif type(IndMatrix) == list:
        for i in range(min(len(FileList), len(IndMatrix))):
            if i % 10 == 0:
                print "Going through file %i of %i" % (i,min(len(FileList), len(IndMatrix)))
            if StartInd:
                AscArray = [[], []]
            else:
                AscArray = []
            M = OpenSaveFile(FileList[i])[1]
            for j in IndMatrix[i]:
                if StartInd:
                    AscArray[1].append(StartPos(M[j])[0])
                    AscArray[0].append(MinXSpan(M[j, TraceInd, :], Criterion))
                else:
                    AscArray.append(MinXSpan(M[j, TraceInd, :], Criterion))   
            AscMatrix.append(AscArray)
        
    else: 
        raise ValueError('IndMatrix type not compatible')
    
    if Flat == False:
        if StartInd != False:
            raise Exception ('Options not compatible!')
        return AscMatrix
    
    elif Flat == True:
        print "Flattening Matrix to Array"
        if StartInd:
            AscFlat = [[], []]
            for i in range(len(AscMatrix)):
                AscFlat[0] += AscMatrix[i][0]
                AscFlat[1] += AscMatrix[i][1]
            AscFlat[0] = np.array(AscFlat[0])
            AscFlat[1] = np.array(AscFlat[1])
            StartList = StartIndList()
            nstart = len(StartList)
            AscFlatInd = []
            for i in range(nstart):
                if StartList[i] in list(AscFlat[1]):
                    AscFlatInd.append(AscFlat[0][AscFlat[1] == StartList[i]])
            return AscFlatInd
    
        else:
            AscFlat = []
            for i in range(len(AscMatrix)):
                AscFlat += AscMatrix[i]
            return AscFlat
    
    
def MinXSpan(Array, Criterion):
    """
    Returns the min x span for required criterion
    Automatically converts min to max
    """
    Array = Array[Array != 0]   # Removes Zero values
    if np.amax(Array) - np.amin(Array) < Criterion:
        asc_span = np.nan
    else:
        #import pdb; pdb.set_trace()
        min_tot = Array.argmin()
        max_tot = Array.argmax()
        if min_tot > max_tot: 
            Array = Array[::-1]
            min_tot = Array.argmin()
            max_tot = Array.argmax()
        asc_span = max_tot - min_tot
        for i in range(max_tot - min_tot):
            where = np.where(Array > (Array[min_tot+i] + Criterion))[0]
            where = where[where > min_tot]
            if where.shape[0] == 0:
                break
            asc_ind = where[0]   
            if (asc_ind - (min_tot+i)) < asc_span:
                asc_span = asc_ind - (min_tot+i)
                a = asc_ind
                b = min_tot+i
    assert asc_span > 0 or np.isnan(asc_span), "asc_span is 0 or negative"
    return asc_span

    
def VertVelMatrix(FileList, TraceInd, IntSpan, StartInd = False, IndMatrix = False, Flat = False):
    """
    Returns 2D list of fastest ascent times for given criterion
    if Index Matrix is given, search only for indexed trajectories
    if Flat == True, returns one long array. To be used for histograms!
    Note: StartInd not compatible with Flat == False!!!
    """
    AscMatrix = []
    if IndMatrix == False:
        for i in range(len(FileList)):
            if i % 10 == 0:
                print "Going through file %i of %i" % (i,len(FileList))
            if StartInd:
                AscArray = [[], []]
            else:
                AscArray = []
            M = OpenSaveFile(FileList[i])[1]
            for j in range(M.shape[0]):
                if StartInd:
                    AscArray[1].append(StartPos(M[j])[0])
                    AscArray[0].append(VertVel(M[j, TraceInd, :], IntSpan))
                else:
                    AscArray.append(VertVel(M[j, TraceInd, :], IntSpan))
            AscMatrix.append(AscArray)
                                
    elif type(IndMatrix) == list:
        for i in range(min(len(FileList), len(IndMatrix))):
            if i % 10 == 0:
                print "Going through file %i of %i" % (i,min(len(FileList), len(IndMatrix)))
            if StartInd:
                AscArray = [[], []]
            else:
                AscArray = []
            M = OpenSaveFile(FileList[i])[1]
            for j in IndMatrix[i]:
                if StartInd:
                    AscArray[1].append(StartPos(M[j])[0])
                    AscArray[0].append(VertVel(M[j, TraceInd, :], IntSpan))
                else:
                    AscArray.append(VertVel(M[j, TraceInd, :], IntSpan))   
            AscMatrix.append(AscArray)
        
    else: 
        raise ValueError('IndMatrix type not compatible')
    
    if Flat == False:
        if StartInd != False:
            raise Exception ('Options not compatible!')
        return AscMatrix
    
    elif Flat == True:
        print "Flattening Matrix to Array"
        if StartInd:
            AscFlat = [[], []]
            for i in range(len(AscMatrix)):
                AscFlat[0] += AscMatrix[i][0]
                AscFlat[1] += AscMatrix[i][1]
            AscFlat[0] = np.array(AscFlat[0])
            AscFlat[1] = np.array(AscFlat[1])
            StartList = StartIndList()
            nstart = len(StartList)
            AscFlatInd = []
            for i in range(nstart):
                if StartList[i] in list(AscFlat[1]):
                    AscFlatInd.append(AscFlat[0][AscFlat[1] == StartList[i]])
            return AscFlatInd
    
        else:
            AscFlat = []
            for i in range(len(AscMatrix)):
                AscFlat += AscMatrix[i]
            return AscFlat


def VertVel(Array, IntSpan):
    """
    Returns the maximum vertical ascent in int span	
    """
    Array = Array[Array != 0]
    maxvel = 0
    for i in range(Array.shape[0]-IntSpan):
        tmp = (Array[i + IntSpan] - Array[i]) / IntSpan
        maxvel = max(maxvel, tmp)
    assert maxvel > 0, "Maximum velocity is 0 of negative"
    return maxvel
    

    
    
# End of submodule!
