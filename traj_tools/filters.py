"""
filters module!!
Submodule of traj_tools
Opening and Filtering saved trajectory files
Naming conventions:
Files = all trajectories
Matrix = one trajectory file (i.e. 1000 trajectories)
Array = one trajectory
"""
import common
import cPickle
import numpy as np
import loadbin
import netCDF4 as nc


class TrajProp(object):
    """
    Class containing relevant information for 
    a set of trajectories.
    
    This class is the heart of this module. 
    It links the raw data and the evaluation and plotting
    functions.
    
    All diagnostics are also performed in a class object.
    
    Inside the class, arrays are numpy.arrays to make masking 
    easier. The output type is as specified under attributes.
    
    
    Parameters
    ----------
    CaseSpecs : Class object
      Class object containing constant specifications of case
      
    
    Attributes
    ----------
    CaseSpecs : Class object
      Class object containing constant specifications of case
    inddict : dictionary
      Dictionary connecting variable name to index in data list
    data : list
      List of numpy.arrays with information about each trajectory

    """
    
    def __init__(self,
                 CaseSpecs):
        
        self.CaseSpecs = CaseSpecs
        
        self._init_prop()
        
        
    def _init_prop(self):
        """
        Initialize dictionary and list of numpy arrays with fixed properties:
        
        * filename: location of file containing respective trajectory
        * trjind: index of trajectory in file
        * startt: starting time of trajectory after model start
        
        """
        
        self.inddict = dict(filename = 0, trjind = 1, startt = 2)
        self.data = [[], [], []]
        
        # Looping over all files, initializing lists instead of np.arrays
        # in order to dynamically add files
        
        nrtot = 0
        
        for i in range(len(CaseSpecs.filelist)):
            rootgrp = nc.Dataset(CaseSpecs.filelist[i], 'r')
            nrtrj = len(rootgrp.dimensions['id'])
            tmpt = int(CaseSpecs.filelist[i].split('-')[1].lstrip('t'))
            self.data[0].extend([CaseSpecs.filelist[i]] * nrtrj)
            self.data[1].extend(range(nr.trj))
            self.data[2].extend([tmpt] * nrtrj)
            nrtot += nrtrj
            
        # Convert lists to np.array
        
        for i in range(len(self.data)):
            self.data[i] = np.array(self.data[i])

        assert (self.data[0].shape[0] == self.data[1].shape[0] == 
                self.data[2].shape[0] == nrtot), \
                    "Error while initializing properties for class: Attribute arrays do not have same shape!"
                               
    def new_filter_function(name, function, *args):
        """
        Add a new filter to data list from custom function. 
        
        
        Parameters
        ----------
        name : string or tuple
          Name of the filter. To be used in inddict.
          If more than one return value in function, use tuple of names
          in correct order.
        function : function
          The Python function to be applied as a filter.
          Has to return a single value for one trajectory array.
          If function returns more than one single value, suffixes have to 
          be specified. The arguments for the function have to be given 
          in *args
          
        ADD EXAMPLE!!!
        """
        
        # Add entries to dictionary
        
        if name == str:
            self.inddict[name] = len(self.data)
        elif name == tuple or name == list:
            self.inddict.update(dict(zip(name), 
                                     range(len(self.data), 
                                           len(self.data) + len(name))))
        else:
            raise TypeError('Invalid input type for name')
        
        # Evaluate variables to be extracted from NetCDF files
        
        filevars = []
        for item in args:
            if item in nc.Dataset(CaseSpecs.filelist[i], 'r').variables:
                filevars.append(item)
        print('Following variables are used:', filevars)
        assert len(filevars) > 0, "No correct variables chosen"
        
        datam = []
        for i in range(len(CaseSpecs.filelist)):
            rootgrp = nc.Dataset(CaseSpecs.filelist[i], 'r')
            for n in range(len(filevars)):
                datam.append(rootgrp.variables[filevars[n]][:, :])
            
            for j in range(len(rootgrp.dimensions['id'])):
                # function call, return output
                if name == str:
                    pass
    
    
    def new_filter_array(name, array):
        """
        Add one or more arrays to data. 
        If more than one array is added, names and data arrays must
        be in tuples in the right order.
        
        Parameters
        ----------
        name : string or tuple
          name(s) of the array(s) to be added
        array : np.array of tuple
          data array(s) to be added         
        
        """
        
        if type(name) == str:
            self.inddict[name] = len(self.data)
            self.data.append(array)
            
        elif type(name) == tuple or type(name) == list:
            self.inddict.update(dict(zip(name), 
                                     range(len(self.data), 
                                           len(self.data) + len(name))))
            self.data.extend(array)
        
    
    
    def new_filter_asc(pspan, tracer = 'P'):
        """
        Adds a new filter to class instance. 
        """
        
        mint = np.zeros(self.data[0].shape[0])
        ascstart = np.zeros(self.data[0].shape[0])
        ascend = np.zeros(self.data[0].shape[0])
        
        
        # Add to dic and data * 3
       
    
    
    
    def apply_filter(self, FilterTUPLE):
        """
        Returns the filtered filename and trajid list.
        Minimum and maximum criteria can be applied for 
        ascent time and vertical velocity.
        
        Parameters
        ----------
        minasct : float
          Lower end of ascent time criterion
        maxasct : float
          Upper end of ascent time criterion
        minvertvel : float
          Lower end of vertical velocity criterion
        maxvertvel : float
          Upper end of vertical velocity criterion
          
        Returns
        -------
        tuple
          tuple containing one list and one list of lists:
          1. list of unique file locations.
          2. list of lists of trajectory IDs for each element in list 1. 
        
        
        """
        
        mask = np.array([True] * self.len)   # Initialize mask 
        
        if minasct == maxasct == minvertvel == maxvertvel == None:
            print("No filter criteria chosen," 
                  "return lists for all trajectories")
        
        if not minasct == None:
            mask &= self.asct >= minasct
        if not maxasct == None:
            mask &= self.asct <= maxasct
        if not minvertvel == None:
            mask &= self.vertvel >= minvertvel
        if not maxvertvel == None:
            mask &= self.vertvel <= maxvertvel
        print mask
        
        uniqueloc = np.unique(self.filename[mask])
        idlist = []
        for i in range(len(uniqueloc)):
            locmask = self.filename == uniqueloc[i]
            idlist.append(int(self.trajid[locmask & mask]))
        
        
        return (uniqueloc, idlist)
            






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
    
    
    
# End of submodule!
