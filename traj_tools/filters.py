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
import fortran.futils as futils




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
    flip = False
    if TraceInd == 7:
        flip = True
    IndList = []
    for i in range(len(Matrix[:, 0, 0])):
        if np.average(Matrix[i, TraceInd, :]) == 0:    # Checks for all Zero arrays"
            print "Given array is all zeros, break!"
            break
        elif MinXSpan(Matrix[i, TraceInd, :], Criterion, flip = flip) <= LenMax and MinXSpan(Matrix[i, TraceInd, :], Criterion, flip = flip) >= LenMin:
            IndList.append(i)
    return IndList

def FilterStartEnd(FileList, TraceInd, Criterion, LenMax, LenMin = 0):
    """
    Filter plus start and end for ascent
    """
    IndMatrix = []
    StartMatrix = []
    EndMatrix = []
    for i in range(len(FileList)):
        print "Completed", i, "out of", len(FileList)
        Matrix = OpenSaveFile(FileList[i])[1]
        IndMatrix.append(StartEndM(Matrix, TraceInd, Criterion, LenMax, LenMin = 0)[0])
        StartMatrix.append(StartEndM(Matrix, TraceInd, Criterion, LenMax, LenMin = 0)[1])
        EndMatrix.append(StartEndM(Matrix, TraceInd, Criterion, LenMax, LenMin = 0)[2])
        del Matrix
    return (IndMatrix, StartMatrix, EndMatrix)


def StartEndM(Matrix, TraceInd, Criterion, LenMax, LenMin = 0):
    """
    for function above
    """
    flip = False
    if TraceInd == 7:
        flip = True
    IndList = []
    StartList = []
    EndList = []
    for i in range(len(Matrix[:, 0, 0])):
        if np.average(Matrix[i, TraceInd, :]) == 0:    # Checks for all Zero arrays"
            print "Given array is all zeros, break!"
            break
        else:
            Span, Start, End = MinXSpan(Matrix[i, TraceInd, :], Criterion, mode = 2,
                                        flip = flip)
            if Span <= LenMax and Span >= LenMin:
                IndList.append(i)
                StartList.append(Start)
                EndList.append(End)
    return (IndList, StartList, EndList)

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


def MinXMatrix(FileList, TraceInd, Criterion, StartInd = False, IndMatrix = False, IndMatrix2 = False, 
               Flat = False):
    """
    Returns 2D list of fastest ascent times for given criterion
    if Index Matrix is given, search only for indexed trajectories
    if Flat == True, returns one long array. To be used for histograms!
    Note: StartInd not compatible with Flat == False!!!
    """
    AscMatrix = []
    flip = False
    if TraceInd == 7:
        flip = True
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
                    AscArray[0].append(MinXSpan(M[j, TraceInd, :], Criterion,
                                                flip = flip))
                else:
                    AscArray.append(MinXSpan(M[j, TraceInd, :], Criterion, 
                                             flip = flip))
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
                if IndMatrix2 != False:
                    if j in IndMatrix2[i]:
                        if StartInd:
                            AscArray[1].append(StartPos(M[j])[0])
                            AscArray[0].append(MinXSpan(M[j, TraceInd, :], Criterion,
                                                        flip = flip))
                        else:
                            AscArray.append(MinXSpan(M[j, TraceInd, :], Criterion,
                                                    flip= flip))
                else:
                    
                    if StartInd:
                        AscArray[1].append(StartPos(M[j])[0])
                        AscArray[0].append(MinXSpan(M[j, TraceInd, :], Criterion,
                                                    flip = flip))
                    else:
                        AscArray.append(MinXSpan(M[j, TraceInd, :], Criterion,
                                                flip= flip))   
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

def MinXSpan(array, yspan, mode = 1, flip = False):
    """
    Returns the minimum time steps needed to conver given criterion.
    Automatically filters out zero and nan values. If yspan is not fulfilled,
    returns np.nan.
    Parameters
    ----------
    array : np.array
    Input array
    yspan : float
    Ascent criterion in y-direction
    flip : bool
    If True array will be flipped (e.g. use for 'P')
    Returns
    -------
    xspan : int
    Time steps for criterion
    istart : int
    Start index of ascent
    istop : int
    End index of ascent
    """

    # Filter out nans and zeros
    array = array[np.isfinite(array)]
    array = array[array != 0]
    # Flip array if needed
    if flip:
        array = -array
    # Check if criterion is met
    #print array
    if array.shape[0] == 0:
        xspan = np.nan
        istart = np.nan
        istop = np.nan
    elif np.amax(array) - np.amin(array) < yspan:
        xspan = np.nan
        istart = np.nan
        istop = np.nan
    else:
        # Use Fortran implementation, use 0 as error values
        xspan, istart, istop = futils.futils.minxspan(array, yspan,
        len(array) + 1, 0, 0)
        # Check if output is correct. NOTE: THIS IS A POTENTIAL BUG!!!
        if (istart < 0) and (istop < 0):
            xspan = np.nan
            istart = np.nan
            istop = np.nan
    #assert ((xspan > 0) and (xspan <= len(array)) and (istart >= 0)
    #and (istop >= 0)), \
    #'One of the minxspan outputs is zero or negative.'
    if mode == 1:
        return xspan
    elif mode ==2:
        return (xspan, istart, istop)

# OLD IMPLEMENTATION 
#def MinXSpan(Array, Criterion, mode = 1):
    #"""
    #Returns the min x span for required criterion
    #Automatically converts min to max
    #"""
    #a = np.nan
    #b = np.nan
    #off = np.where(Array != 0)[0][0]
    #Array = Array[Array != 0]   # Removes Zero values
    #if np.amax(Array) - np.amin(Array) < Criterion:
        #asc_span = np.nan
        #a = np.nan
        #b = np.nan
    #else:
        ##import pdb; pdb.set_trace()
        #min_tot = Array.argmin()
        #max_tot = Array.argmax()
        #if min_tot > max_tot: 
            #Array = Array[::-1]
            #min_tot = Array.argmin()
            #max_tot = Array.argmax()
        #asc_span = max_tot - min_tot
        #for i in range(max_tot - min_tot):
            #where = np.where(Array > (Array[min_tot+i] + Criterion))[0]
            #where = where[where > min_tot]
            #if where.shape[0] == 0:
                
                #break
            #asc_ind = where[0]   
            #if (asc_ind - (min_tot+i)) < asc_span:
                #asc_span = asc_ind - (min_tot+i)
                #a = asc_ind
                #b = min_tot+i
    #assert asc_span > 0 or np.isnan(asc_span), "asc_span is 0 or negative"
    #if mode == 1:
        #return asc_span
    #elif mode == 2:
        ##print a, Array.shape[0] ,Array.shape[0] - a
        #return (asc_span, Array.shape[0] - a + off, Array.shape[0] - b + off)

    
def VertVelMatrix(FileList, TraceInd, IntSpan, StartInd = False, IndMatrix = False, Flat = False, mode = 1):
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
                    AscArray[0].append(VertVel(M[j, TraceInd, :], IntSpan, mode = mode))
                else:
                    AscArray.append(VertVel(M[j, TraceInd, :], IntSpan,mode =mode))
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
                    AscArray[0].append(VertVel(M[j, TraceInd, :], IntSpan, mode = mode))
                else:
                    AscArray.append(VertVel(M[j, TraceInd, :], IntSpan, mode = mode))   
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


def VertVel(Array, IntSpan, mode = 1):
    """
    Returns the maximum vertical ascent in int span	
    """
    Array = Array[Array != 0]
    maxvel = 0
    if mode ==2:
        for i in range(Array.shape[0]):
            maxvel = max(maxvel, Array[i])
    elif mode == 3:
        mn = np.amin(Array)
        mni = np.argmin(Array)
        mx = np.amax(Array)
        mxi = np.argmax(Array)
        slope = (mx - mn) / Array.shape[0]
        for i in range(Array.shape[0]):
            maxvel += abs(Array[i] - (i * slope))
    elif mode ==4:
        for i in range(Array.shape[0]-IntSpan):
            tmp = np.amax(Array[i:i+IntSpan]) - Array[i]
            #print np.amax(Array[i:i+IntSpan]), Array[i], tmp
            maxvel = max(maxvel, tmp)
                          
            
    elif mode == 1:
        for i in range(Array.shape[0]-IntSpan):
            tmp = (Array[i + IntSpan] - Array[i]) / IntSpan
            maxvel = max(maxvel, tmp)
    #assert maxvel > 0, "Maximum velocity is 0 of negative"
    return maxvel

    

    
    
# End of submodule!
