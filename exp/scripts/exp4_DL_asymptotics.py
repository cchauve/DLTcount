import sys
import importlib

from trees   import *
from DLTcount import *
from DLTasymptotics import *

if __name__=="__main__":
    kMin       = 3   # min size (#leaves) of considered trees
    kMax       = 25  # max size (#leaves) of considered trees
    nbTrees    = 100 # Number of trees per size
    inputPrefix = '../trees/trees' 
    outputPrefix = '../asymptotics/asymptotics_DL' # output file

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
                
    # Asymptotics for DL histories
    for k in TREES.keys():
        output = open(outputPrefix+"_"+str(k),"w")
        output.write("#size_tree\ttree_index\ttree\tgrowth factor\n")
        output.flush()
        for i in range(0,nbTrees):
            tree = TREES[k][i]
            labelTree(tree)
            stree = tree.asNewick()
            (sysEqs,sysVars) = getSystem(tree,MODEL={'D':True,'L':True,'T':False})
            solveSystem(sysEqs,sysVars,str(k)+'_'+str(i))
            tmp_solver = importlib.import_module('tmp_solver_'+str(k)+'_'+str(i), package=None)
            output.write(str(k)+'\t'+str(i)+'\t'+stree+'\t'+str(tmp_solver.SOL)+'\n')
            os.system('rm tmp_solver_'+str(k)+'_'+str(i)+'.py')
    outputFile.close()

