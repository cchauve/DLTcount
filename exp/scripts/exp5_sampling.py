# Sampling histories

from trees   import *
from DLTcount import *

import math
import random

if __name__=="__main__":
    k   = 16    # size (#leaves) of considered trees
    n   = 30    # size of histories
    mk  = 50    # Number of random trees
    mn  = 10000 # Number of samples
    model = sys.argv[1]
    seed   = 101

    if model == 'UDL':
        MODEL  = {'D':True, 'L': True, 'T':False}
        ranked = False
    elif model == 'UDLT':
        MODEL = {'D':True, 'L': True, 'T':True}
        ranked = False
    elif model == 'RDL':
        MODEL  = {'D':True, 'L': True, 'T':False}
        ranked = True
    elif model == 'RDLT':
        MODEL = {'D':True, 'L': True, 'T':True}
        ranked = True
        
    random.seed(seed)

    output = open('../sampling/samples_'+str(k)+'_'+str(n)+'_'+model,'w')
    
    for i in range(0,mk):    
        utree  = randomOrderedBinaryTree(k)
        labelTree(utree)
        output.write('# S_'+str(i)+'\t'+utree.asNewick())
        if ranked:
            (ranking,rtree) = rankTreeRandomly(utree)
            labelTree(rtree)
            output.write('\t'+str(ranking)+'\n')
            stree = rtree
        else:
            output.write('\n')
            stree = utree

        allLeaves = []
        allNodes  = stree.allNodes()
        for node in allNodes:
            if node.isLeaf():
                allLeaves.append(node.getID())
            
        S,H,D,T = fillMatrices(stree,n,MODEL)

        for j in range(0,mn):
            history = randGen(stree,'H',n,S,H,D,T,MODEL)
            nZ = {}
            for q in allLeaves:
                nZ[q] = history.count('Z'+str(q))
            strHistory = str(history)
            nD = strHistory.count('D')
            nL = strHistory.count('L')
            nT = strHistory.count('T')
            for q in allLeaves:                    
                output.write('Z'+str(q)+':'+str(nZ[q])+' ')
            output.write('D:'+str(nD)+' ')
            output.write('L:'+str(nL)+' ')
            output.write('T:'+str(nT)+'\n')
            
    output.close()
