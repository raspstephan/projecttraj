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
        
        self.inddict = dict(startt = 0)
        self.filename = []
        self.trajid = []
        self.data = [[]]
        self.filtdict = dict()
        self.filtlist = []
        
        # Looping over all files, initializing lists instead of np.arrays
        # in order to dynamically add files
        
        nrtot = 0   # Counter for assert check below
        
        for i in range(len(self.CaseSpecs.filelist)):
            rootgrp = nc.Dataset(self.CaseSpecs.filelist[i], 'r')
            nrtrj = len(rootgrp.dimensions['id'])
            tmpt = self.CaseSpecs.filelist[i].split('/')[-1]
            tmpt = int(tmpt.split('_')[1].lstrip('t'))
            self.filename.extend([self.CaseSpecs.filelist[i]] * nrtrj)
            self.trajid.extend(range(nrtrj))
            self.data[0].extend([tmpt] * nrtrj)
            nrtot += nrtrj
            
        # Convert lists to np.array
        
        self.filename = np.array(self.filename)
        self.trajid = np.array(self.trajid)
        self.data[0] = np.array(self.data[0])

        assert (self.data[0].shape[0] == self.filename.shape[0] == 
                self.trajid.shape[0] == nrtot), \
                "Error while initializing properties for class: Attribute arrays do not have same shape!"
                               
    
    
    def new_prop_array(self, name, array):
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
            assert (array.shape[0] == self.data[0].shape[0]), \
                    'Array shaped do not match.'
            assert (name not in self.inddict), 'Name already exists.' 
            self.inddict[name] = len(self.data)
            self.data.append(array)
            print(name, 'has been added.')
            
            
        elif type(name) == tuple or type(name) == list:
            for i in range(len(name)):
                assert (name[i] not in self.inddict), \
                        'One of the names already exists.'
                assert (array[i].shape[0] == self.data[0].shape[0]), \
                        'One or all of the array shapes do not match.'
            self.inddict.update(dict(zip(name), 
                                     range(len(self.data), 
                                           len(self.data) + len(name))))
            self.data.extend(array)
            print(name, 'have been added.')
        
    
    
    def new_prop_asc(self, yspan, tracer = 'P'):
        """
        Adds a new ascent filter to data.
        
        Parameters
        ----------
        yspan : float
          Ascent criterion in y-direction
        tracer : str (default = 'P')
          COSMO name of y-axis variable
        
        """
        
        # Update dictionary
        code = tracer + str(yspan)
        self.inddict[code] = len(self.data)
        self.inddict[code + '_start'] = len(self.data) + 1
        self.inddict[code + '_stop'] = len(self.data) + 2
        
        ascdata = minasct(self.CaseSpecs.filelist, yspan, tracer)
        assert (ascdata[0].shape[0] == self.data[0].shape[0]), \
                'Array shapes do not match. Look for error in source code.'
        self.data.extend(ascdata)
    
       
    def create_filter(self, name, filters):
        """
        Adds a filter to filtlist. Adds name to filtdict.
        
        Parameters
        ----------
        name : string
          Name of given filter
        filters : tuple
          List of filter tuples, see example.
        
        Examples
        --------
        >>> filters = (('P600', 2880, 0), ('P300', 1000, 200))
        
        """
        
        assert (name not in self.filtdict), 'Filter name already exists.'
        
        mask = np.array([True] * self.data[0].shape[0])   # Initialize mask 
        
        for i in range(len(filters)):

            ind = self.inddict[filters[i][0]]   # Get index in self.data
            
            if len(filters[i]) == 3:
                mx = filters[i][1]
                mn = filters[i][2]
                
                mask &= self.data[ind] <= mx
                mask &= self.data[ind] >= mn
                
            elif len(filters[i]) == 2:
                crit = filters[i][1]
                
                mask &= self.data[ind] == crit
                
            else:
                raise Exception('Wrong input for filter tuple.')
       
        self.filtdict[name] = len(self.filtlist)
        self.filtlist.append(mask)
        
        
    def _mask_conv(self, mask):
        """
        Converts mask to iterable lists. To be used in plotting functions, etc.
        
        Paremeters
        ----------
        mask : np.array
          Mask to be converted
        
        Returns
        -------
        uniqueloc : list
          List of all unique save locations.
          Use as first iterator.
        idlist : list
          List of list. Trajectory IDs for every element in 
          uniquelist.yspan
          Use as second iterator.
          
        """
        
        uniqueloc = np.unique(self.filename[mask])
        idlist = []
        for i in range(len(uniqueloc)):
            locmask = self.filename == uniqueloc[i]
            print self.trajid[locmask & mask]
            idlist.append((self.trajid[locmask & mask]))

        return uniqueloc, idlist


def minasct(filelist, yspan, tracer):
    """
    Calculate minimum ascent time for all trajectories from NetCDF files.
    
    Parameters
    ----------
    filelist : list
      List of saved NetDCF file locations
    yspan : float
      Ascent criterion in y-direction
    tracer : str 
      COSMO name of y-axis variable
     
    Returns
    -------
    ascdata : tuple
      Tuple containing three np.arrays:
      * Minimum ascent time
      * Index of ascent start
      * Index of ascent stop
    
    """
    
    # Initialize lists, convert to np.array later
    asct = []
    ascstart = []
    ascstop = []
     
    for i in range(len(filelist)):
        mat = nc.Dataset(filelist[i], 'r').variables[tracer][:, :]
        for j in range(mat.shape[1]):
            asctup = minxspan(mat[:, j], yspan)
            asct.append(asctup[0])
            ascstart.append(asctup[1])
            ascstop.append(asctup[2])
    assert (len(asct) == len(ascstart) == len(ascstop)), \
            'Array lengths do not match'
    ascdata = (np.array(asct), np.array(ascstart), np.array(ascstop))
    return ascdata
            

def minxspan(Array, Criterion, mode = 2):
    """
    Returns the min x span for required criterion
    Automatically converts min to max
    
    ALGORITHM HAS TO BE CHECKED, IMPROVED!!!
    """
    a = np.nan
    b = np.nan
    off = np.where(Array != 0)[0][0]
    Array = Array[Array != 0] # Removes Zero values
    if np.amax(Array) - np.amin(Array) < Criterion:
        asc_span = np.nan
        a = np.nan
        b = np.nan
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
    #assert asc_span > 0 or np.isnan(asc_span), "asc_span is 0 or negative"
    if mode == 1:
        return asc_span
    elif mode == 2:
        #print a, Array.shape[0] ,Array.shape[0] - a
        return (asc_span, Array.shape[0] - a + off, Array.shape[0] - b + off)





# Old code ######################################

    #def new_filter_function(name, function, *args):
        #"""
        #Add a new filter to data list from custom function. 
        
        
        #Parameters
        #----------
        #name : string or tuple
          #Name of the filter. To be used in inddict.
          #If more than one return value in function, use tuple of names
          #in correct order.
        #function : function
          #The Python function to be applied as a filter.
          #Has to return a single value for one trajectory array.
          #If function returns more than one single value, suffixes have to 
          #be specified. The arguments for the function have to be given 
          #in *args
          
        #ADD EXAMPLE!!!
        #"""
        
        ## Add entries to dictionary
        
        #if name == str:
            #self.inddict[name] = len(self.data)
        #elif name == tuple or name == list:
            #self.inddict.update(dict(zip(name), 
                                #range(len(self.data), 
                                      #len(self.data) + len(name))))
        #else:
            #raise TypeError('Invalid input type for name')
        
        ## Evaluate variables to be extracted from NetCDF files
        
        #filevars = []
        #for item in args:
            #if item in nc.Dataset(CaseSpecs.filelist[i], 'r').variables:
                #filevars.append(item)
        #print('Following variables are used:', filevars)
        #assert len(filevars) > 0, "No correct variables chosen"
        
        #datam = []
        #for i in range(len(CaseSpecs.filelist)):
            #rootgrp = nc.Dataset(CaseSpecs.filelist[i], 'r')
            #for n in range(len(filevars)):
                #datam.append(rootgrp.variables[filevars[n]][:, :])
            
            #for j in range(len(rootgrp.dimensions['id'])):
                ## function call, return output
                #if name == str:
                    #pass
    

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


#def IndData(FileList, IndList, MaxFile="default"):
    #"""
    #Gives Matrix of data for all indexed trajectories
    #"""
    #f = open(FileList[0], "r")
    #Spec = cPickle.load(f)
    #f.close()
    #NrTrj = 0
    #if MaxFile == "default":
        #MaxFile = min(len(FileList), len(IndList))
    #for i in range(MaxFile):
        #NrTrj += len(IndList[i])
    #FilterM = np.zeros((NrTrj, Spec[1]+3, Spec[5])) 
    #ind = 0
    #for i in range(MaxFile):
        #print "Opening file", i
        #RawM = OpenSaveFile(Fdef minxspan(Array, Criterion, mode = 2):
    #"""
    #Returns the min x span for required criterion
    #Automatically converts min to max
    #"""
    #a = np.nan
    #b = np.nan
    #off = np.where(Array != 0)[0][0]
    #Array = Array[Array != 0] # Removes Zero values
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
        #return (asc_span, Array.shape[0] - a + off, Array.shape[0] - b + off)ileList[i])[1]
        #for j in IndList[i]:
            #FilterM[ind, :, :] = RawM[j, :, :]
            #ind += 1    
    #return FilterM


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


#def MinXMatrix(FileList, TraceInd, Criterion, StartInd = False, IndMatrix = False, Flat = False):
    #"""
    #Returns 2D list of fastest ascent times for given criterion
    #if Index Matrix is given, search only for indexed trajectories
    #if Flat == True, returns one long array. To be used for histograms!
    #Note: StartInd not compatible with Flat == False!!!
    #"""
    #AscMatrix = []
    #if IndMatrix == False:
        #for i in range(len(FileList)):
            #if i % 10 == 0:
                #print "Going through file %i of %i" % (i,len(FileList))
            #if StartInd:
                #AscArray = [[], []]
            #else:
                #AscArray = []
            #M = OpenSaveFile(FileList[i])[1]
            #for j in range(M.shape[0]):
                #if StartInd:
                    #AscArray[1].append(StartPos(M[j])[0])
                    #AscArray[0].append(MinXSpan(M[j, TraceInd, :], Criterion))
                #else:
                    #AscArray.append(MinXSpan(M[j, TraceInd, :], Criterion))
            #AscMatrix.append(AscArray)
                                
    #elif type(IndMatrix) == list:
        #for i in range(min(len(FileList), len(IndMatrix))):
            #if i % 10 == 0:
                #print "Going through file %i of %i" % (i,min(len(FileList), len(IndMatrix)))
            #if StartInd:
                #AscArray = [[], []]
            #else:
                #AscArray = def minxspan(Array, Criterion, mode = 2):
    #"""
    #Returns the min x span for required criterion
    #Automatically converts min to max
    #"""
    #a = np.nan
    #b = np.nan
    #off = np.where(Array != 0)[0][0]
    #Array = Array[Array != 0] # Removes Zero values
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
        #return (asc_span, Array.shape[0] - a + off, Array.shape[0] - b + off)[]
            #M = OpenSaveFile(FileList[i])[1]
            #for j in IndMatrix[i]:
                #if StartInddef minxspan(Array, Criterion, mode = 2):
    #"""
    #Returns the min x span for required criterion
    #Automatically converts min to max
    #"""
    #a = np.nan
    #b = np.nan
    #off = np.where(Array != 0)[0][0]
    #Array = Array[Array != 0] # Removes Zero values
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
        #return (asc_span, Array.shape[0] - a + off, Array.shape[0] - b + off):
                    #AscArray[1].append(StartPos(M[j])[0])
                    #AscArray[0].append(MinXSpan(M[j, TraceInd, :], Criterion))
                #else:
                    #AscArray.append(MinXSpan(M[j, TraceInd, :], Criterion))   
            #AscMatrix.append(AscArray)
        
    #else: 
        #raise ValueError('IndMatrix type not compatible')
    
    #if Flat == False:
        #if StartInd != False:
            #raise Exception ('Options not compatible!')
        #return AscMatrix
    
    #elif Flat == True:
        #print "Flattening Matrix to Array"
        #if StartInd:
            #AscFlat = [[], []]
            #for i in range(len(AscMatrix)):
                #AscFlat[0] += AscMatrix[i][0]
                #AscFlat[1] += AscMatrix[i][1]
            #AscFlat[0] = np.array(AscFlat[0])
            #AscFlat[1] = np.array(AscFlat[1])
            #StartList = StartIndList()
            #nstart = len(StartList)
            #AscFlatInd = []
            #for i in range(nstart):
                #if StartList[i] in list(AscFlat[1]):
                    #AscFlatInd.append(AscFlat[0][AscFlat[1] == StartList[i]])
            #return AscFlatInd
    
        #else:
            #AscFlat = []
            #for i in range(len(AscMatrix)):
                #AscFlat += AscMatrix[i]
            #return AscFlat
    
    
def MinXSpanOLD(Array, Criterion):
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
