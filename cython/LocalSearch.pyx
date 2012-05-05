#cython: profile=True

import numpy as np
cimport numpy as np
import WalshAnalysis as wal
import itertools as it
import tool as tl
import random
import string
import math
import copy
import sys
import os
import time
import pdb
from operator import itemgetter
from cpython cimport bool
cimport cython
from cython.parallel import prange, parallel, threadid

from libcpp.vector cimport vector
from libcpp.set cimport set
#from libcpp.map cimport map
#from libcpp.pair cimport pair
from libc.stdlib cimport malloc
from cython.operator cimport dereference as deref, preincrement as inc

# standard cimport file libc/stdlib.pxd
class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

cdef class Indiv:
    cdef public char* bit
    cdef public double fit

ctypedef struct InfBit:
    vector[int]* arr
    int WI

ctypedef struct InTer:
    set[int]* arr
    set[int]* WI

ctypedef struct Was:
    int* arr
    double w

ctypedef struct ComArr:
    int** arr
    int size

cdef class LocalSearch:
    cdef InTer** Inter
    cdef vector[InfBit]** infectBit
    cdef Was* WAS
    cdef double* sumArr
    cdef double** C
    cdef list improveA
    cdef object func
    cdef object model
    cdef int MaxFit
    cdef int dim
    cdef double threshold
    cdef int fitEval
    #cdef dict lookup
    cdef ComArr** lookup
    cdef Indiv oldindiv, bsf

    def __init__(self,object model,int MaxFit,int dim):
        self.func = model.compFit
        self.model = model
        self.MaxFit = MaxFit
        self.dim = dim
        self.threshold = 1e-15

    cpdef run(self, fitName, minimize, restart, compM):
        if compM == 'walWalkNext':
            if fitName == 'fit':
                return self.runFitSwalkNext(fitName, minimize, restart)
#            elif fitName == 'mean':
#                return self.runMeanSCwalkNext(fitName, minimize, restart)

    cdef runFitSwalkNext(self,fitName, minimize, restart):
        """ 
        next descent local search running on S
        """
        cdef: 
            int initC = 1
            float updateC = 0

            float start = 0
            float descT = 0
            float pertT = 0
            float updatePertT = 0
            float updateT = 0
            int walkLen = 10

        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.bsf = Indiv()
        self.bsf.bit = self.oldindiv.bit 
        self.bsf.fit = self.oldindiv.fit
#        print 'BEGIN: InitWal'
        self.initWal()
#        print 'END  : InitWal'
        self.genImproveS(minimize)
        self.model.WA = []

        self.fitEval = 0

        initT = os.times()[0] - start

#        print 'BEGIN: search'

        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.nextDesc()
#            print improveN, bestI
#            print self.improveA
            descT = descT + os.times()[0] - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = os.times()[0]
                    diff = self.walk(fitName, minimize, walkLen)
                    pertT = pertT + os.times()[0] - start

                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i]
                        self.update(i)
                        start = os.times()[0]
                        self.updateDeep(i)
                        updateT = updateT + os.times()[0] - start
                        self.updateWAS(i)
                        #print 'BEGIN: updatePertImprS'
                        self.updatePertImprS(i, minimize)
                        #print self.improveA
                        #print 'END  : updatePertImprS'

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI]
#                print 'BEGIN: update', bestI
                self.update(bestI)
#                print 'END  : update'
#                print 'BEGIN: updateWAS'
                start = os.times()[0]
                self.updateDeep(bestI)
                updateT = updateT + os.times()[0] - start
                self.updateWAS(bestI)
#                print 'END  : updateWAS'
#                print 'BEGIN: updateImprS'
                self.updateImprS(bestI, minimize)
#                print self.improveA
#                print 'END  : updateImprS'
                self.fitEval = self.fitEval + 1
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'
#            print 'END  : search'
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}


    cdef walk(self, fitName, minimize, length):
        # update the bsf solution
        if fitName == 'fit' and minimize == True :
            if self.bsf.fit > self.oldindiv.fit:
                self.bsf = Indiv()
                self.bsf.bit = self.oldindiv.bit 
                self.bsf.fit = self.oldindiv.fit
        elif fitName == 'fit' and minimize == False :
            if self.bsf.fit < self.oldindiv.fit:
                self.bsf = Indiv()
                self.bsf.bit = self.oldindiv.bit 
                self.bsf.fit = self.oldindiv.fit
        elif minimize == True :
            if self.bsf.fitG > self.oldindiv.fitG:
                self.bsf = Indiv()
                self.bsf.bit = self.oldindiv.bit 
                self.bsf.fit = self.oldindiv.fit
        elif minimize == False :
            if self.bsf.fitG < self.oldindiv.fitG:
                self.bsf = Indiv()
                self.bsf.bit = self.oldindiv.bit 
                self.bsf.fit = self.oldindiv.fit

        flipBits = random.sample(xrange(self.dim), length)
        for i in flipBits:
            if self.oldindiv.bit[i] == '1':
                self.oldindiv.bit[i] = '0'
            else:
                self.oldindiv.bit[i] = '1'
        return flipBits



    cdef evalPop(self, Indiv indiv):
        indiv.fit = self.func(indiv.bit)
        self.fitEval = self.fitEval + 1
        return indiv

    cdef genImproveS(self,minimize):
        """
        generate the index of best neigh according to sumArr only (surrogate of fitness)
        """
        # check improving move
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
                self.improveA.append(i) 

    cdef nextDesc(self):
        """
        find the next improving move by the similar update trick
        """
        if not self.improveA:
            return False, None

        cdef int I, bestI


        # randomly pick an improving move, which takes only constant time 
        I = random.choice(xrange(len(self.improveA)))
        bestI = self.improveA[I]

        return True, bestI

    cdef initWal(self):
        """ 
        1. 
        compute the sum array for the first time, according to the initial solution
        2. 
        initialize a full matrix, Coincidence, C_ij, which is a symmetric one, and diagonal is empty (as we don't want to recompute SumArr), empty whenever i >= j.
        Two options for implementing C matrix
            a. pre-compute C matrix 
            b. construct it on the fly using dict()
        3. 
        initialize interaction list for each variable
        4.
        initialize a dict() of interaction structure, where interactive bits and the index of WAS (walsh including sign)
        """
        cdef int i,j,k
        cdef InTer* inter 
        cdef vector[InfBit]* vectPtr
        cdef InfBit* strPtr
        cdef ComArr* comb

        cdef Was* was

        #self.sumArr = np.zeros(self.dim)
        self.sumArr = <double*>malloc(self.dim * sizeof(double))
        for i in xrange(self.dim):
            # initialize sumArr
            self.sumArr[i] = 0

        # allocate total space for storing pointers
        self.infectBit = < vector[InfBit]** > malloc(sizeof(void *) * self.dim)
        for i in xrange(self.dim):
            # assign pointers to the allocated space
            vectPtr = new vector[InfBit]()
            self.infectBit[i] = vectPtr


        self.C = <double **>malloc(sizeof(double *) * self.dim)
        for i in xrange(self.dim) :
            self.C[i] = <double *> malloc(sizeof(double) * self.dim)

        #self.WAS = np.tile(Was, len(self.model.w.keys())) # Walsh coefficients with sign, represented in Array
        self.WAS = <Was* > malloc(sizeof(Was)* len(self.model.w.keys()))

#        self.lookup = dict()
        #self.lookup = new map[int,ComArr*]() 
        self.lookup = <ComArr**> malloc(sizeof(ComArr*)*self.dim)
        for i in xrange(self.dim):
            self.lookup[i] = NULL


#        self.Inter = [[] for i in range(self.dim)]
        self.Inter = < InTer** > malloc(sizeof(void *)*self.dim)
        for i in xrange(self.dim):
            self.Inter[i] = NULL
        
        for i in xrange(self.dim):
            for j in xrange(self.dim):
                self.C[i][j] = 0
        for i in xrange(len(self.model.WA)):
            W = int(math.pow(-1, binCount(self.model.WA[i].arr, self.oldindiv.bit))) * self.model.WA[i].w
            #self.WAS[i] = Struct(arr = self.model.WA[i].arr, w = W)
#            self.WAS[i] = Was()
#            self.WAS[i].arr = self.model.WA[i].arr
#            self.WAS[i].w = W
            
            was = <Was *>malloc(sizeof(Was))
            was[0].arr = <int *>malloc(sizeof(int)*len(self.model.WA[i].arr))
            for j in xrange(len(self.model.WA[i].arr)):
                was[0].arr[j] = self.model.WA[i].arr[j]
            was[0].w = W
            self.WAS[i] = was[0]

            comb = self.genComb(len(self.model.WA[i].arr))

            for j in self.model.WA[i].arr:
                self.sumArr[j] = self.sumArr[j] + W
                if len(self.model.WA[i].arr)>1: # for at least order Walsh terms
                    #if not self.Inter[j]: # the entry of i doesn't exist yet
                    if self.Inter[j] == NULL:
#                        self.Inter[j] = InTer()
#                        self.Inter[j].arr = []
#                        self.Inter[j].WI = []
                        inter = <InTer*> malloc(sizeof(InTer))
                        inter[0].arr = new set[int]()
                        inter[0].WI = new set[int]()
                        self.Inter[j] = inter
                        
#                    for k in self.model.WA[i].arr:
#                        if k != j and k not in self.Inter[j].arr:
#                            self.Inter[j].arr.append(k)
#                    if  i not in self.Inter[j].WI:
#                        self.Inter[j].WI.append(i)

                    for k in self.model.WA[i].arr:
                        if k != j :
                            self.Inter[j].arr.insert(k)
                    self.Inter[j].WI.insert(i)


#                print 'BEGIN: self.model.WA[%d].arr' %(i)
                # add list of order >= 3 Walsh terms for the purpose of updating C matrix
                if len(self.model.WA[i].arr) >= 3:
#                    infBit = InfBit()
#                    infBit.arr = self.model.WA[i].arr
#                    infBit.WI = i

                    strPtr = <InfBit *> malloc(sizeof(InfBit))
                    strPtr.WI = i
#                    strPtr.arr = <vector[int]*> malloc(sizeof(vector[int]))
                    #print strPtr.WI, strPtr.arr[0].size(), strPtr.arr.size()
                    strPtr.arr = new vector[int]()
                    for k in self.model.WA[i].arr:
                    #    print 'BEGIN: strPtr.arr.push_back(%d)' %(k)
                        strPtr.arr[0].push_back(k)
                    #    print strPtr.arr[0].size(), strPtr.arr.size()
#                        print 'END  : strPtr.arr.push_back(%d)' %(k)

#                    if not self.infectBit[j]: 
#                        self.infectBit[j] = [infBit]
#                    else :
#                        self.infectBit[j].append(infBit)
                    
                    #print 'BEGIN: self.infectBit[%d].push_back(strPtr[0])' %(j)
                    self.infectBit[j][0].push_back(strPtr[0])
                    #print 'END  : self.infectBit[%d].push_back(strPtr[0])' %(j)
#                print 'END  : self.model.WA[%d].arr' %(i)

#            for l in comb: # for each list in comb
#                j0 = self.model.WA[i].arr[int(l[0])]
#                j1 = self.model.WA[i].arr[int(l[1])]
#                self.C[j0][j1] = self.C[j0][j1] + W

            for l in xrange(comb.size):
                j0 = self.model.WA[i].arr[comb.arr[l][0]]
                j1 = self.model.WA[i].arr[comb.arr[l][1]]
                self.C[j0][j1] = self.C[j0][j1] + W



#    cdef ComArr* genComb(self,int N) nogil:
#        """ 
#        Generate C_N^2 sequence, index are stored, because they are more general, Implemented in an *incremental* fashion.
#        """
#        cdef int c, j, i, counter
#        cdef ComArr* ptr
#        cdef map[int, ComArr*].iterator it
#        cdef pair[int, ComArr*]* p
#
#        it = self.lookup[0].find(N)
#
#        if it != self.lookup.end(): # the key exists before
#            return deref(it).second
#
#        else : # the key not found
#            c = biomial(N, 2)
#            counter = 0
#
#            ptr = <ComArr*> malloc(sizeof(ComArr))
#            ptr.arr = <int**> malloc(sizeof(int*)*c)
#            ptr.size = c
#
#            for i in xrange(c):
#                ptr.arr[i] = <int*> malloc(sizeof(int)*2)
#
#            for i in xrange(N):
#                for j in xrange(i+1,N):
#                    ptr.arr[counter][0] = i
#                    ptr.arr[counter][1] = j
#                    counter = counter + 1
#
#            p = new pair[int,ComArr *](N, ptr)
#            self.lookup.insert(p[0])
#            return ptr


    #cdef ComArr* genComb(self,int N):
    cdef ComArr* genComb(self,int N) nogil:
        """ 
        Generate C_N^2 sequence, index are stored, because they are more general, Implemented in an *incremental* fashion.
        """
        cdef int c, j, i, counter
        cdef ComArr* ptr

        if self.lookup[N] != NULL: # the key exists before
#            print N
#            print "exist"
            return self.lookup[N]
        else : # the key not found
#            print N
#            print "not found"
            c = biomial(N, 2)
            counter = 0

            ptr = <ComArr*> malloc(sizeof(ComArr))
            ptr.arr = <int**> malloc(sizeof(int*)*c)
            ptr.size = c

            for i in xrange(c):
                ptr.arr[i] = <int*> malloc(sizeof(int)*2)

            for i in xrange(N):
                for j in xrange(i+1,N):
                    ptr.arr[counter][0] = i
                    ptr.arr[counter][1] = j
                    counter = counter + 1
            self.lookup[N] = ptr
            return ptr

#    @cython.boundscheck(False)
    cdef update(self, int p):
        """
        By keeping track of coincidence matrix, 
        Cij stands for S_i(y_j) = S_i(x) - C_ij
        partially update the Sum Array and self.WAS, given the bit which is changed
        """
        self.sumArr[p] = - self.sumArr[p]
        #cdef int k, k0, k1, i, ii, len1, len2
        cdef int i,ii
        cdef int len1
        cdef list comb
        cdef set[int].iterator it

        if self.Inter[p]!=NULL:
            """ iterative over self.Inter[p].arr """
            #print 'Not NULL', p
            #for i in prange(len(self.Inter[p].arr), nogil= True):
            it = self.Inter[p].arr.begin()
            #print 'begin'
            while it != self.Inter[p].arr.end():
                ii = deref(it)
                #print ii
                if ii < p:
                    self.sumArr[ii] = self.sumArr[ii] - 2*self.C[ii][p]
                    self.C[ii][p] = - self.C[ii][p]
                else:
                    self.sumArr[ii] = self.sumArr[ii] - 2*self.C[p][ii]
                    self.C[p][ii] = - self.C[p][ii]
                inc(it)
            

    cdef void updateDeep(self, int p):
        cdef int i, k0, k1
        cdef vector[int].iterator it
        cdef vector[int] arr
        cdef InfBit I

        # update the rest of elements in C matrix
        if self.infectBit[p].size() != 0:
            for i in prange(self.infectBit[p].size(), nogil=True):
#            for i in xrange(self.infectBit[p].size()):
                I = self.infectBit[p][0][i]
                arr = I.arr[0]
                it = arr.begin()
                while it != arr.end():
                    if deref(it) == p:
                        arr.erase(it)
                        break
                    inc(it)
#                with gil:
#                    comb = self.genComb(arr.size())
#                    for k in xrange(len(comb)):
#                        k0 = arr[int(comb[k][0])]
#                        k1 = arr[int(comb[k][1])]
#                        self.C[k0][k1] = self.C[k0][k1] - 2 * self.WAS[I.WI].w
                comb = self.genComb(arr.size())
                for k in xrange(comb.size):
                    k0 = arr[int(comb.arr[k][0])]
                    k1 = arr[int(comb.arr[k][1])]
                    self.C[k0][k1] = self.C[k0][k1] - 2 * self.WAS[I.WI].w

    cdef updateImprS(self, int p, bool minimize):
        cdef int i,I
        cdef set[int].iterator it

        self.improveA.remove(p)
        if  self.Inter[p]!=NULL:

            it = self.Inter[p].arr.begin()
            while it != self.Inter[p].arr.end():
                I = deref(it)
                if (minimize == True and self.sumArr[I] > self.threshold) or (minimize == False and self.sumArr[I]< self.threshold ):
                    if I not in self.improveA:
                        self.improveA.append(I)
                elif I in self.improveA:
                    self.improveA.remove(I)
                inc(it)


#            for i in xrange(len(self.Inter[p].arr)): 
#                I = self.Inter[p].arr[i]
#                if (minimize == True and self.sumArr[I] > self.threshold) or (minimize == False and self.sumArr[I]< self.threshold ):
#                    if I not in self.improveA:
#                        self.improveA.append(I)
#                elif I in self.improveA:
#                    self.improveA.remove(I)

    cdef updatePertImprS(self, int p, bool minimize):
        cdef int i,I
        cdef set[int].iterator it

        if  self.Inter[p]!=NULL:
            it = self.Inter[p].arr.begin()
            while it != self.Inter[p].arr.end():
                I = deref(it)
                if (minimize == True and self.sumArr[I] > self.threshold) or (minimize == False and self.sumArr[I]< self.threshold ):
                    if I not in self.improveA:
                        self.improveA.append(I)
                elif I in self.improveA:
                    self.improveA.remove(I)
                inc(it)

#            for i in xrange(len(self.Inter[p].arr)): 
#                I = self.Inter[p].arr[i]
#                if (minimize == True and self.sumArr[I] > self.threshold) or (minimize == False and self.sumArr[I]< self.threshold ):
#                    if I not in self.improveA:
#                        self.improveA.append(I)
#                elif I in self.improveA:
#                    self.improveA.remove(I)

        if (minimize == True and self.sumArr[p] > self.threshold) or (minimize == False and self.sumArr[p]< self.threshold ):
            if p not in self.improveA:
                self.improveA.append(p)
        elif p in self.improveA:
            self.improveA.remove(p)

    cdef updateWAS(self, int p):
        cdef int i, I
        cdef set[int].iterator it

        if self.Inter[p]!=NULL:
            it = self.Inter[p].WI.begin()
            while it != self.Inter[p].WI.end():
                i = deref(it)
                self.WAS[i].w = - self.WAS[i].w
                inc(it)

#            for i in self.Inter[p].WI:
#                self.WAS[i].w = - self.WAS[i].w


    def printW(self):
        """
        print all the original walsh terms
        """
        for k in self.model.w.keys():
            print k, '%.3f' % (self.model.w[k])

    def printWsort(self):
        """
        print all the sorted Walsh terms
        """
#        sorted(self.model.w.values())
        a = sorted(self.model.w.iteritems(), key=itemgetter(1), reverse=True)
        for i in a:
            print i[0], '%.3f' %(i[1])

    def printWAsort(self):
        for i in self.WAsort:
            print i.arr,'\t\t%.3f' %(i.w)

    def printWA(self):
        """
        print all walsh terms with array representation
        """
        for i in range(len(self.model.WA)):
            print self.model.WA[i].arr, self.model.WA[i].w

    def printC(self):
        a = 0
        c = 0
        for i in range(self.dim):
            for j in range(i+1,self.dim):
                a = a + 1
                if self.C[i][j] != 0:
                    print '1',
                    c = c + 1
                else:
                    print '0',
            print 
        print c,a, c/float(a)

#    def printInter(self):
#        for i in self.Inter:
#            print i,self.Inter[i].arr

    def genCombOne(self, odd, order):
        """ generate the number of possible ones, 
            given whether it should be odd or not,
            and the order of Walsh coefficients
        """
        if odd == True:
            if (odd, order) not in self.combLoopup:
                self.combLoopup[odd, order] = 2*np.array(xrange((order+1)/2)) + 1
        else:
            if (odd, order) not in self.combLoopup:
                self.combLoopup[odd, order] = 2*np.array(xrange((order)/2+1))

        return copy.deepcopy(self.combLoopup[odd, order])

    def genPossBit(self, odd, arr):
        comb = self.genCombOne(odd, len(arr))
        iters = []
        for i in comb:
            for j in it.combinations(arr, i):
                iters.append(list(j))
        return iters

cdef Indiv initIndiv(int dim):
    """ initial the search inidividual with random bit string """
    cdef char *randBitStr = <char*>malloc(dim+1 * sizeof(char))
    for j in xrange(dim):
        if random.random()<0.5:
            randBitStr[j] = '0'
        else:
            randBitStr[j] = '1' 

    randBitStr[dim] = '\0'

    cdef Indiv indiv = Indiv()
    indiv.bit = randBitStr
    indiv.fit = 0
    return indiv

cdef int binCount(list arr,str bit):
    """
    count the one bit of union self.model.WA[i].arr and bit
    """
    cdef int s = 0
    cdef int i
    for i in arr:
        if bit[i] == '1':
            s = s + 1
    return s

cdef double sumC(double * a, int d):
    cdef double s = 0

    for i in xrange(d):
        s = s+a[i]

    return s

@cython.cdivision(True)
cdef int biomial(int N, int K) nogil:
    """ compute the combination of N choose K """
    return factorial(N)/( factorial(K) * factorial(N-K) )

cdef int factorial(int N) nogil:
    """ compute N! """
    cdef int c, fact = 1
    for c in xrange(1,N+1):
        fact = fact*c
    return fact
