import nkLandscape as nk
import nkqLandscape as nkq
import WalshAnalysis as wal
import MAXSAT as mx
import LocalSearch as ls
import tool as tl
import os
import numpy as np
import random
import math
import time
import pdb
import sys


def main(argv):
    """ consider as a minimization problem """
    tl.checkParam(argv)

    rseed = 0
    nameOfDir = './result/'
    runtimeDir = './runtime/'
    waltimeDir = './walshtime/'
    traceDir = './trace/'
    prefixNK = '../benchmark/NK/'
    prefixNKQ = '../benchmark/NKQ/'

    random.seed(rseed)

    compMeth = tl.getArgv(argv) # bf(brute force) / wal (walsh analysis)
    probName = tl.getArgv(argv)
    algoName = tl.getArgv(argv)
    fitName = tl.getArgv(argv) # fit/mean/std

    #if compMeth == 'wal' and fitName != 'mean':
    #    print 'ERROR: Walsh analysis can only be applied to compute mean'
    #    sys.exit()

    inst = int(tl.getArgv(argv))
    s = tl.getArgv(argv) # get the setting for population size
    n = int(tl.getArgv(argv))
    if probName != 'SAT':
        k = int(tl.getArgv(argv))

    #maxFit = 1 * n
    maxFit = 100 * n
    #maxFit = 0
    runs = 1
    popSize = 50 # always keep popSize to even number
    q = 0

    if algoName.find('LS') != -1:
        popSize = 1

    if probName == 'SAT':
        """ with SAT, we are forced to set n to 100 """

        """ 
        TODO : 
            need to perform multiple runs for each instance 
        """

        model = mx.MAXSAT()
        res = []

        model.setInstance(inst)
        print 'Instance', inst

        if algoName.find('LS') != -1:
            algo = ls.LocalSearch(model.compFit, maxFit, n)

        tAll = np.zeros(runs)

        for i in range(runs):
            start = os.times()[0]
            if algoName == 'LS':
                res.append(algo.run(fitName, minimize = False, restart = False))
            elif algoName == 'rLS':
                res.append(algo.run(fitName, minimize = False, restart = True))
            tAll[i] = os.times()[0] - start

        if probName == 'SAT':
            nameOfF = nameOfDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'.txt'
        f = open(nameOfF, 'w')
        for i in range(len(res)):
            if fitName != 'fit':
                print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
            else:
                print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
        f.close()

        """ store runtime to files """
        if probName == 'SAT':
            nameOfF = runtimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'.txt'

        f = open(nameOfF, 'w')
        for i in range(runs):
            print >>f,"%g" % (tAll[i])
        f.close()

    else:
        res = []

        if probName == 'NK':
            model = nk.NKLandscape(n,k,prefixNK+'NK-N'+str(n)+'-K'+str(k)+'-I'+str(inst))
            #model = nk.NKLandscape(n,k)
        elif probName == 'NKQ':
            q = int(tl.getArgv(argv))
            model = nkq.NKQLandcape(n, k, q, prefixNKQ+'NKQ-N'+str(n)+'-K'+str(k)+'-I'+str(inst)+'-Q'+str(q))
            #model = nkq.NKQLandcape(n, k, q)

        if compMeth == 'walWalk' or compMeth == 'walRest' or compMeth == 'supm' or compMeth == 'bitImp' or compMeth == 'walSearch' or compMeth == 'checkOptWal' or compMeth == 'checkHyper' or compMeth == 'checkHyperRank' or compMeth == 'hyperSearch' or compMeth == 'hyperSqSearch' or compMeth == 'hyperWalSearch' or compMeth == 'walWalkNext' or compMeth == 'walRestNext':
            start = os.times()[0]
            # Walsh analysis
            w = model.WalshCofLinearLinklist()
            walTime = os.times()[0] - start

            start = os.times()[0]
            if compMeth == 'checkHyper' or compMeth == 'checkHyperRank' or compMeth == 'hyperSearch':
                model.genHyperVote()
            elif compMeth == 'hyperSqSearch':
                model.genHyperSqVote()
            elif compMeth == 'hyperWalSearch':
                model.genHyperWalVote()
            hyperTime = os.times()[0] - start

            # count the number of interative bits
            # model.countInterBits()

            """ store runtime to files """
            if probName == 'NKQ':
                nameOfF = waltimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
            elif probName == 'NK':
                nameOfF = waltimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

            f = open(nameOfF, 'w')
            print >>f,"%g\t%g" % (walTime,hyperTime) 
            f.close()

    #    bit,fit = tl.compFit(model)
    #    a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1]) 
    #    print 'opti\n', a[0][0], a[0][1]
    #    print

    #    for i in a:
    #        print i[0], '%.2f' %(i[1])

    #    for i in zip(bit,fit):
    #        print i[0],'%.3f' %(i[1])
            
    #    print 'bit',bit
    #    print 'fit',fit
    #    print 'mean',np.mean(fit)
    #    print 'w', w

    #    numOpt = lo.localOpt(bit, fit)
    #    print numOpt

        if algoName.find('LS') != -1:
            algo = ls.LocalSearch(model, maxFit, n)

        tAll = np.zeros(runs)
        for i in range(runs):
    #        print 'run', i, ':probName', probName, 'algoName', algoName, 'fitName', fitName, 'I', inst, 'n', n, 'k', k 
            start = os.times()[0]
            if algoName == 'LS':
                res.append(algo.run(fitName, minimize = True, restart = False,compM = compMeth ))
            elif algoName == 'rLS':
                res.append(algo.run(fitName, minimize = True, restart = True,compM = compMeth))
            tAll[i] = os.times()[0] - start

    #    trace = res[0]['trace']
    #    for i in trace:
    ##        print 'Eval', i.fitEval, 'fit', i.fit
    #        print 'Eval', i.fitEval, 'fit', i.fit, 'fitG', i.fitG
    #
    #    plt.plot([i.fitEval for i in trace],[i.fit for i in trace],'.-')
    #    plt.plot([i.fitEval for i in trace],[i.fitG for i in trace],'.-')
    #    plt.show()

        """ store results to files """
        if probName == 'NKQ':
            nameOfF = nameOfDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
        elif probName == 'NK':
            nameOfF = nameOfDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

    #    """ print the mean over multiple runs """
    #    r = np.zeros(runs)
    #    for i in range(runs):
    #        r[i] = res[i]['sol']
    #    print np.mean(r)
        f = open(nameOfF, 'w')
        for i in range(runs):
            if fitName != 'fit':
                print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
            else:
                print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
        f.close()

        """ store trace to files """
        if probName == 'NKQ':
            nameOfF = traceDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
        elif probName == 'NK':
            nameOfF = traceDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'
        f = open(nameOfF, 'w')
        for i in range(runs):
              print >>f,"%g\t%g" % (res[i]['initC'], res[i]['updateC'])
        f.close()


        """ store runtime to files """
        if probName == 'NKQ':
            nameOfF = runtimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
        elif probName == 'NK':
            nameOfF = runtimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

        f = open(nameOfF, 'w')
        print >>f,"All\t\tinit\t\tdesc\t\tpert\t\tupdate\t\tupdatePert\t"
        for i in range(runs):
            print >>f,"%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e\t%0.2e" % (tAll[i], res[i]['init'],res[i]['descT'], res[i]['pertT'], res[i]['updateT'], res[i]['updatePertT'])
        f.close()

        print nameOfF, 'Finish'
        print

if __name__ == '__main__':
    main(sys.argv)
