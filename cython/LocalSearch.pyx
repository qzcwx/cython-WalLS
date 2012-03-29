# cython: profile=True

import numpy as np
cimport numpy as np
#import matplotlib.pyplot as plt
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

# standard cimport file libc/stdlib.pxd

cdef extern from "stdlib.h":
    void free(void* ptr)
    void* malloc(size_t size)
    void* realloc(void* ptr, size_t size)
    # ...

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

cdef class Indiv:
    cdef public char* bit
    cdef public double fit

cdef class LocalSearch:
    #cdef np.ndarray sumArr
    cdef double* sumArr
    cdef double** C
    cdef list improveA
    cdef np.ndarray WAS
    cdef object func
    cdef object model
    cdef int MaxFit
    cdef int dim
    cdef double threshold
    cdef int fitEval

    cdef dict lookup, infectBit, Inter
    cdef Indiv oldindiv, bsf

    def __init__(self,object model,int MaxFit,int dim):
        self.func = model.compFit
        self.model = model
        self.MaxFit = MaxFit
        self.dim = dim
        self.threshold = 1e-15

    def run(self, fitName, minimize, restart, compM):
        if compM == 'walWalkNext':
            if fitName == 'fit':
                return self.runFitSwalkNext(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeanSCwalkNext(fitName, minimize, restart)

    def runFitSwalkNext(self,fitName, minimize, restart):
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
        self.initWal()
        self.genImproveS(minimize)
        self.model.WA = []

        self.fitEval = 0

        initT = os.times()[0] - start

        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.nextDesc()
            descT = descT + os.times()[0] - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = os.times()[0]
                    diff = self.walk(fitName, minimize, walkLen)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i]
                        self.update(i)
                        self.updateWAS(i)
                        self.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start
                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]
                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI]
                self.update(bestI)
                self.updateWAS(bestI)
                self.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}


    def runMeanSCwalkNext(self,fitName, minimize, restart):
        """ 
        next descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0

        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = initIndivNeigh(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sumC(self.sumArr,self.dim))
        self.indiv.fitG = self.oldindiv.fitG
        self.genImproveSC(minimize)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.model.WA = []

        initC = 0
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start
        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.nextDesc()
            descT = descT + os.times()[0] - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = os.times()[0]
                    diff = self.walk(fitName, minimize, walkLen)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i] 
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sumC(self.sumArr, self.dim))

                        self.update(i)
                        self.updateSC(i)
                        self.updateWAS(i)
                        self.updatePertImprSC(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]
                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI] 
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sumC(self.sumArr,self.im))
                self.update(bestI)
                self.updateSC(bestI)
                self.updateWAS(bestI)
                self.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}


    def diffBits(self, a, b):
        diff = []
        for i in xrange(self.dim):
            if a[i] != b[i]:
                diff.append(i)
        return diff

    def walk(self, fitName, minimize, length):
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


    def neighbors(self):
        neighs = []
        for j in xrange(self.dim):
            # flip the jth bit in bit-string
            neighStr = np.copy(self.oldindiv.bit)
            if neighStr[j] == '1':
                neighStr[j] = '0'
            else:
                neighStr[j] = '1'
            neighs.append( neighStr )
        return np.array(neighs)

    def evalPop(self, Indiv indiv):
        indiv.fit = self.func(indiv.bit)
        self.fitEval = self.fitEval + 1
        return indiv

    def evalPopNeigh(self, indiv, fitName, minimize):
        """ evaluate the individual itself """
        indiv.fit = self.func(indiv.bit)
        self.fitEval = self.fitEval + 1
        """ evaluate all neighborhood """
        fitN = np.zeros(self.dim)
        for j in range(self.dim):
            # flip the jth bit in bit-string
            neighStr = np.copy(indiv.bit)
            if neighStr[j] == '1':
                neighStr[j] = '0'
            else:
                neighStr[j] = '1'
            fitN[j] = self.func(neighStr)
        if fitName == 'mean':
            indiv.fitG = np.mean(fitN)
        elif minimize == True :
            indiv.fitG = np.mean(fitN) - np.std(fitN)
        elif minimize == False :
            indiv.fitG = np.mean(fitN) + np.std(fitN)

        return copy.deepcopy(indiv)

    def genImproveS(self,minimize):
        """
        generate the index of best neigh according to sumArr only (surrogate of fitness)
        """
        # check improving move
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
                self.improveA.append(i) 

    def nextDesc(self):
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

    def genImproveSC(self,minimize):
        """
        generate the index of best neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move 
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                self.improveA.append(i)


    def genMeanNext(self,minimize):
        """
        generate the index of next improving neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move 
        improve = False
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                self.improveA.append(i)
                improve = True

        if improve == False:
            return False, None

        bestI = random.choice(self.improveA)
        return True, bestI


    def updateMeanNext(self, p, minimize):
        self.improveA.remove(p)

        if p in self.Inter:
            for i in self.Inter[p].arr:
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None

        bestI = random.choice(self.improveA)
        return True, bestI

    def initWal(self):
        """ 
        1. 
        compute the sum array for the first time, according to the initial solution
        2. 
        initialize a full matrix, Coincidence, C_ij, which is a symmetric one, and diagonal is empty (as we don't want to recompute SumArr), empty whenever i >= j.
        Two options for implementing C matrix
            a. precompute C matrix 
            b. construct it on the fly using dict()
        3. 
        initialize interaction list for each variable
        4.
        initialize a dict of interaction structure, where interactive bits and the index of WAS (walsh including sign)
        """
        cdef int i,j,k
        #self.sumArr = np.zeros(self.dim)
        self.sumArr = <double*>malloc(self.dim * sizeof(double))

        self.C = <double **>malloc(sizeof(double *) * self.dim)

        for i in xrange(self.dim) :
            self.C[i] = <double *> malloc(sizeof(double) * self.dim)

        self.WAS = np.tile(Struct(arr = [], w = 0), len(self.model.w.keys()))# Walsh coefficients with sign, represented in Array
        self.lookup = dict()
        self.infectBit = dict()
        self.Inter = dict()

        for i in xrange(self.dim):
            self.sumArr[i] = 0
        
        for i in xrange(self.dim):
            for j in xrange(self.dim):
                self.C[i][j] = 0

        for i in xrange(len(self.model.WA)):
            W = int(math.pow(-1, binCount(self.model.WA[i].arr, self.oldindiv.bit))) * self.model.WA[i].w
            self.WAS[i] = Struct(arr = self.model.WA[i].arr, w = W)
            comb = self.genComb(len(self.model.WA[i].arr))

            for j in self.model.WA[i].arr:
                self.sumArr[j] = self.sumArr[j] + W
                if len(self.model.WA[i].arr)>1: # for at least order Walsh terms
                    if j not in self.Inter: # the entry of i doesn't exist yet
                        self.Inter[j] = Struct(arr=[], WI=[])

                    for k in self.model.WA[i].arr:
                        if k != j and k not in self.Inter[j].arr:
                            self.Inter[j].arr.append(k)
                    if i not in self.Inter[j].WI:
                        self.Inter[j].WI.append(i)

                # add list of order >= 3 Walsh terms for the purpose of updating C matrix
                if len(self.model.WA[i].arr) >= 3:
                    if j not in self.infectBit: 
                        self.infectBit[j] = [Struct(arr=self.model.WA[i].arr, WI=i)]
                    else :
                        self.infectBit[j].append(Struct(arr=self.model.WA[i].arr, WI=i))

            for l in comb: # for each list in comb
                j0 = self.model.WA[i].arr[int(l[0])]
                j1 = self.model.WA[i].arr[int(l[1])]
                self.C[j0][j1] = self.C[j0][j1] + W


    def initSC(self):
        # compute the SC array
        self.SC = np.zeros(self.dim)
        self.Z = np.zeros(self.dim)
        self.orderC = np.zeros((self.dim,self.dim))

        for p in range(self.dim):

            phi = np.zeros(self.model.k+1)
            if p in self.Inter:
                for i in self.Inter[p].WI:
                    order = len(self.WAS[i].arr)
                    phi[order-1] = phi[order-1] + self.WAS[i].w

            self.Z[p] = self.sumArr[p]
            for i in range(1, self.model.k+1):
                if phi[i] != 0:
                    self.Z[p] = self.Z[p] + i * phi[i]

            self.SC[p] = self.sumArr[p] - 2/float(self.dim) * self.Z[p]

        for i in range(len(self.WAS)):
            lenArr = len(self.WAS[i].arr)
            comb = self.genComb(lenArr)
            for j in comb:
                j0 = self.WAS[i].arr[int(j[0])]
                j1 = self.WAS[i].arr[int(j[1])]
                self.orderC[j0,j1] = self.orderC[j0,j1] + lenArr * self.WAS[i].w

    def genComb(self,N):
        """ 
        Generate C_k^0 sequence, index are stored, because they are more general, Implemented in an *incremental* fashion.
        """
        if N in self.lookup.keys(): # the key exists before
            return self.lookup[N]
        else : 
            comb =  []
            c = 0
            for i in range(N):
                for j in [ k for k in range(N) if k > i]:
                    arr = np.zeros(2)
                    arr[0] = i
                    arr[1] = j
                    comb.append(arr)
                    c = c + 1    
            self.lookup[N] = comb
            return comb

    def update(self, p):
        """
        By keeping track of coincidence matrix, 
        Cij stands for S_i(y_j) = S_i(x) - C_ij
        partially update the Sum Array and self.WAS, given the bit which is changed
        """
        self.sumArr[p] = - self.sumArr[p]
        cdef int ii, k, k0, k1
        cdef object i
        cdef list arr

        if p in self.Inter:
            for ii in self.Inter[p].arr:
                if ii < p:
                    self.sumArr[ii] = self.sumArr[ii] - 2*self.C[ii][p]
                    self.C[ii][p] = - self.C[ii][p]
                else:
                    self.sumArr[ii] = self.sumArr[ii] - 2*self.C[p][ii]
                    self.C[p][ii] = - self.C[p][ii]

        # update the rest of elements in C matrix
        if p in self.infectBit.keys():
            for i in self.infectBit[p]:
                arr = i.arr[:]
                arr.remove(p)
                comb = self.genComb(len(arr))
                for k in xrange(len(comb)):
                    k0 = arr[int(comb[k][0])]
                    k1 = arr[int(comb[k][1])]
                    self.C[k0][k1] = self.C[k0][k1] - 2 * self.WAS[i.WI].w

    def updateImprS(self, int p, bool minimize):
        cdef int i,I
        self.improveA.remove(p)
        if p in self.Inter:
            for i in xrange(len(self.Inter[p].arr)): 
                I = self.Inter[p].arr[i]
                if (minimize == True and self.sumArr[I] > self.threshold) or (minimize == False and self.sumArr[I]< self.threshold ):
                    if I not in self.improveA:
                        self.improveA.append(I)
                elif I in self.improveA:
                    self.improveA.remove(I)

    def updatePertImprS(self, int p, bool minimize):
        cdef int i,I
        if p in self.Inter:
            for i in xrange(len(self.Inter[p].arr)): 
                I = self.Inter[p].arr[i]
                if (minimize == True and self.sumArr[I] > self.threshold) or (minimize == False and self.sumArr[I]< self.threshold ):
                    if I not in self.improveA:
                        self.improveA.append(I)
                elif I in self.improveA:
                    self.improveA.remove(I)

        if (minimize == True and self.sumArr[p] > self.threshold) or (minimize == False and self.sumArr[p]< self.threshold ):
            if p not in self.improveA:
                self.improveA.append(p)
        elif p in self.improveA:
            self.improveA.remove(p)

    def updateWAS(self, int p):
        cdef int i, I
        if p in self.Inter:
            #            for i in xrange(len(self.Inter[p].WI)):
#                I = self.Inter[p].WI[i]
#                self.WAS[I].w = - self.WAS[I].w
            for i in self.Inter[p].WI:
                self.WAS[i].w = - self.WAS[i].w

    def updateSC(self, p):
        self.SC[p] = - self.SC[p]
        self.Z[p] = - self.Z[p]

        #update Z array
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if i < p :
                    self.Z[i] = self.Z[i]  - 2* self.orderC[i,p]
                    self.orderC[i,p] = - self.orderC[i,p]
                else :
                    self.Z[i] = self.Z[i]  - 2* self.orderC[p,i]
                    self.orderC[p,i] = - self.orderC[p,i]
                self.SC[i] = self.sumArr[i] - 2/float(self.dim) * self.Z[i]

        if p in self.infectBit.keys():
            for i in self.infectBit[p]:
                arr = i.arr[:]
                arr.remove(p)
                lenArr = len(arr)
                comb = self.genComb(lenArr)
                for k in range(len(comb)):
                    k0 = arr[int(comb[k][0])]
                    k1 = arr[int(comb[k][1])]
                    self.orderC[k0,k1] = self.orderC[k0,k1] - 2 * (lenArr + 1)* self.WAS[i.WI].w

    def updateImprSC(self, p, minimize):
        self.improveA.remove(p)
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

    def updatePertImprSC(self, p, minimize):
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)
        if (minimize == True and self.SC[p] > self.threshold) or (minimize == False and self.SC[p]<self.threshold):
            if p not in self.improveA:
                self.improveA.append(p)
        elif p in self.improveA:
            self.improveA.remove(p)

    def neighWal(self):
        """ 
        generate neighborhoods and compute their fitnesses (mean of neighs) by Walsh Analysis
        """
        neighs = self.neighbors()
        neighIndiv = np.tile(Struct(fit = 0, fitG=0, bit = '0' ), (self.dim))
        for i in range(self.dim):
            neighIndiv[i] = Struct(fit = 0, bit = neighs[i]) 
        return neighIndiv

    def neighFitWal(self):
        """ 
        generate neighborhoods and compute their fitnesses (real one) by Walsh Analysis
        """
        neighs = self.neighbors()
        neighIndiv = np.tile(Struct(fit = 0,  bit = '0' ), (self.dim))
        for i in range(self.dim):
            neighIndiv[i] = Struct(fit = 0, bit = neighs[i]) 
        return neighIndiv

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

    def printInter(self):
        for i in self.Inter:
            print i,self.Inter[i].arr

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

cpdef Indiv initIndiv(int dim):
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

def initIndivNeigh(int dim):
    """ initial the search inidividual with random bit string """
    indiv = Struct(fit = 0, fitG = 0, bit = '0') 
    randBitStr = []
    for j in xrange(dim):
        if random.random()<0.5:
            randBitStr.append('0')
        else:
            randBitStr.append('1')
    indiv = Struct( fit = 0, fitG = 0, bit = randBitStr)
    return indiv

cpdef int binCount(list arr,str bit):
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
