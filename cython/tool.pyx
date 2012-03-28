import LocalSearch as ls
#import matplotlib.pyplot as plt
import random
import numpy as np
import math
import sys

argvCount = 0

def getArgv(argv):
    global argvCount
    arg = argv[argvCount]
    argvCount = argvCount + 1
    return arg

def globalOpt(model):
    """ find the global optimum on the fly """
    n = model.getN()
    for i in range(int(math.pow(2,n))):
        bit = bin(i)
        bit = bit[2:]
        if len(bit) < n:
            bit = (n - len(bit))*'0' + bit
        if i == 0:
            bestFit = model.compFit(bit)
        else:
            curFit = model.compFit(bit)
            if curFit < bestFit:
                bestFit = curFit
                bestBit = bit
    return bestFit, bestBit

def compFit(model):
    n = model.getN()
    fit = np.zeros(math.pow(2,n))
    bitStr = genSeqBits(n)
    for i in range(int(math.pow(2,n))):
       fit[i] = model.compFit(bitStr[i])
    return bitStr, fit


def checkParam(argv):
    if len(argv) == 1:
        print 'Usage: python demo.py [ComputeMethod] [NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [overwrite] [I] [PopSize] [N] [K] [Q]'
        print """
Constrains: 
1) for SAT, I = [1,100]
2) for SAT, N = 100, K/Q = None
3) for NK/NKQ, I = [0,9]
4) for NK, Q = None
5) for NKQ, Q = {2, 4, 8, 16}
6) for LS,  PopSize = 1
7) [ComputeMethod] = {bf (brute force), walWalk (walsh coefficient with random walk), walRest (walsh coefficient with random walk), supm (super move), bitImp (bit impact), walSearch (Search completely relying on Walsh terms), null (run nothing)}. 
    If 'walSearch', rLS, fit. 
        """
        sys.exit()
 

def evalSol(s,model,fitName,minimize):
    if fitName == 'fit':
        return model.compFit(s)
    else :
        fitN = np.zeros(len(s))
        for j in range(len(s)):
            # flip the jth bit in bit-string
            neighStr = np.copy(s)
            if neighStr[j] == '1':
                neighStr[j] = '0'
            else:
                neighStr[j] = '1'
            fitN[j] = model.compFit(neighStr)
        if fitName == 'mean':
            return np.mean(fitN)
        elif minimize == True :
            return np.mean(fitN) - np.std(fitN)
        elif minimize == False :
            return np.mean(fitN) + np.std(fitN)

def genSeqBits(n):
    bitStr = []
    for i in range(int(math.pow(2,n))):
       bit = bin(i)
       bit = bit[2:]
       if len(bit) < n:
           bit = (n - len(bit))*'0' + bit
       bitStr.append(bit)
    return bitStr

def printSep():
    print 80*'*'
