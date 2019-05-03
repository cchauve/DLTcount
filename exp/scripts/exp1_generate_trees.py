# Experiment 1. 
# We generate trees in a range size (kMin=argv[1],kMax=argv[2])
# For each size we generate nbTrees=argv[3] trees
# The first tree is always the caterpillar
# If k is a power of 2 the second tree is the complete binary tree
# The other trees are random binary trees

from trees   import *
from DLTcount import *

import math
import random

def is_power2(num):
    return num != 0 and ((num & (num - 1)) == 0)

if __name__=="__main__":
    kMin       = 3   # min size (#leaves) of considered trees
    kMax       = 25  # max size (#leaves) of considered trees
    nbTrees    = 100 # number of trees per size
    seed       = 20  # random seed
    outputPrefix = '../trees/trees' # output file

    random.seed(seed)
    TREES       = {}
    for k in range(kMin,kMax+1):
        TREES[k] = {}
        caterpillar = buildCaterpillar(k)
        labelTree(caterpillar)
        idxTrees = 0
        TREES[k][idxTrees]= caterpillar
        idxTrees += 1
        if is_power2(k):
            completeTree = buildCompleteTree(int(math.log(k,2)))
            labelTree(completeTree)
            TREES[k][idxTrees] = completeTree
            idxTrees += 1
        for i in range(idxTrees,nbTrees):
            randomTree = randomOrderedBinaryTree(k)
            labelTree(randomTree)
            TREES[k][i] = randomTree
    for k in TREES.keys():
        output = open(outputPrefix+"_"+str(k),"w")
        output.write("#size tree\tindex tree\ttree\n")
        for i in range(0,nbTrees):
            output.write(str(k)+"\t"+str(i)+"\t"+TREES[k][i].asNewick()+"\n")
            output.flush()
        output.close()
