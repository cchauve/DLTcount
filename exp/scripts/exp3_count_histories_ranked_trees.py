# Experiment 2. 
# Counting histories for ranked species trees, where 10 random rankings are generated per tree

from trees   import *
from DLTcount import *
import random

if __name__=="__main__":
    kMin       = 3   # min size (#leaves) of considered trees
    kMax       = 25  # max size (#leaves) of considered trees
    nbTrees    = 100 # number of trees per size
    nMax       = 50  # max history size
    nbRankings = 10  # random ranking per trees
    seed       = 20  # random seed
    inputPrefix = '../trees/trees' 
    outputPrefix = '../ranked/results_ranked' # output file

    random.seed(seed)

    # Reading trees
    TREES = {}
    for k in range(kMin,kMax+1):
        TREES[k] = {}
        inputFile = open(inputPrefix+'_'+str(k)).readlines()
        for line in inputFile:
            if line[0] != "#":
                line1 = line.rstrip().split('\t')
                size   = line1[0]
                treeID = int(line1[1])
                tree   = newick2Tree(line1[2])
                TREES[k][treeID] = tree
                
    # Counting DL and DLT histories
    for k in TREES.keys():
        output = open(outputPrefix+"_"+str(k),"w")
        output.write("#size_tree\ttree_index\tranking_type\tDL/DLT\ttree/ranking\tnumber_of_histories\n")
        output.flush()
        for i in range(0,nbTrees):
            unrankedTree = TREES[k][i]
            labelTree(unrankedTree)
            for r in range(1,nbRankings+1):
                s1 = str(k)+"\t"+str(i)+"\tU\tDL\t"+unrankedTree.asNewick()+" "
                s2 = str(k)+"\t"+str(i)+"\tU\tDLT\t"+unrankedTree.asNewick()+" "
                (ranking,rankedTree) = rankTreeRandomly(unrankedTree)
                rankingStr = printRanking(ranking) 
                s1 += rankingStr+"\t"
                s2 += rankingStr+"\t"
                for n in range(1,nMax+1):
                    S1,H1,D1,T1 = fillMatrices(rankedTree,n,MODEL={'D':True,'L':True,'T':False},)
                    nbHistories1 = int(H1[rankedTree.getID()][n],)
                    s1 += str(nbHistories1)+" "
                    S2,H2,D2,T2 = fillMatrices(rankedTree,n,MODEL={'D':True,'L':True,'T':True})
                    nbHistories2 = int(H2[rankedTree.getID()][n])
                    s2 += str(nbHistories2)+" "
                output.write(s1+"\n")
                output.write(s2+"\n")
                output.flush()
        output.close()
