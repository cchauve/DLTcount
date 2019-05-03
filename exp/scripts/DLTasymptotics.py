__author__  = "Cedric Chauve, Michael Wallner"
__date__    = "April 24, 2019"

# Compute the asymptotics growth factor for a given species tree

import sys
import numpy as np
import sympy as syp
import os
from sympy.matrices import Matrix, zeros
from trees import *

# ---------------------------------------------------------------------------
# Given a tree t, creates a system of equation whose solution gives
# the exponential growth factor of the counting sequence for the number
# of histories conditional to t
#
# naive algorithm: create 1 equation per node
# i.e. does not optimize for identical subgrammars due to identical subtrees
#
# Input: 
# t .............. root of species tree: Node class
# MODEL .......... evolutionary model
# Output: 
# mySys .......... String with system of equations
# M .............. Symbolic matrix of the system
def getSystem(t,MODEL={'D':True,'L':True,'T':False}):
    # Tree t nodes
    allNodes = t.allNodes()
    K        = len(allNodes)
    
    # Create symbolic variables needed for the determinant
    symVarStr  = ','.join(['b'+str(i) for i in range(K)])+',z'
    symVarList = syp.symbols(symVarStr)
    symVarDict = {str(v):v for v in symVarList}

    M     = zeros(K) 
    mySys = []    
    
    for u in allNodes:
        i  = u.getID()
        bi = symVarDict['b'+str(i)]
        eq = 0
        
        # Duplication
        if MODEL['D']:
            eq += -bi+bi*bi
            M[i,i] = 2*symVarList[i]-1

        if u.isLeaf():
            eq += symVarDict['z']
        elif u.isUnary():
            c  = u.getChild().getID()
            bc = symVarDict['b'+str(c)]
            eq += bc # Unary speciation
            M[i,c] = symVarList[c]+1
        else:
            l,r = u.getLeft().getID(),u.getRight().getID()
            bl  = symVarDict['b'+str(l)]
            br  = symVarDict['b'+str(r)]
            eq += bl*br # Speciation
            if MODEL['L']:
                eq += bl+br # SpeciationLoss
            M[i,l] = symVarList[r]+1
            M[i,r] = symVarList[l]+1

        # HGT
        if MODEL['T']:
            if u.getTime() >= 0:
                receivers = u.getContemporary()
            else:
                receivers = u.getIncomparable()
            for v in receivers:
               j  = v.getID()
               bj = symVarDict['b'+str(j)]
               eq += bi*bj
               M[i,i] += symVarList[j]
               M[i,j]  = symVarList[i]

        mySys.append(eq)
        
    detM = M.det(method='det_LU')
    mySys.append(detM)
    
    return((mySys,symVarList))

def solveSystem(mySys,symVarList,suffix=''):
    solver_name = 'tmp_solver_'+suffix+'.py'
    solver_file = open(solver_name,'w')

    solver_file.write('import numpy as np\n')
    solver_file.write('import scipy\n')
    solver_file.write('import scipy.optimize as opt\n')

    solver_file.write('\n')
    solver_file.write('def ftree(X):\n')
    VAR = '\t('
    for v in symVarList:
        VAR +=str(v)+','
    VAR += ') = X\n'
    solver_file.write(VAR.replace(',)',')'))
    EQS = '\tf = '+str(mySys)
    solver_file.write(EQS+'\n')
    solver_file.write('\treturn(f)\n')
    solver_file.write('\n')
    RHS = 'RHS = ['
    for i in range(len(symVarList)):
        RHS +='0,'
    RHS += ','
    solver_file.write(RHS.replace(',,',']')+'\n')
    solver_file.write('\n')
    solver_file.write('def tmp_solve(f):\n')
    solver_file.write('\tsol = opt.newton_krylov(f,RHS,f_tol=1e-14)\n')
    solver_file.write('\treturn(1.0/sol['+str(len(symVarList)-1)+'])\n')
    solver_file.write('SOL = tmp_solve(ftree)\n')

# ----------------------------------------------------------------------
USAGE = 'DLTasymptotics <tree> <MODEL>\n'+'tree = random k | rrandom k | caterpillar k | rcaterpillar k | complete h | rcomplete h | newick string\n'+'\trandom k = random binary tree with k leaves\n'+'\trrandom k = randomly ranked random binary tree with k leaves\n'+'\tcaterpillar k = caterpillar with k leaves\n'+'\trcaterpillar k = randomly ranked caterpillar with k leaves\n'+'\tcomplete h = complete binary tree with 2^h leaves\n'+'\trcomplete h = randomly ranked complete binary tree with 2^h leaves\n'+'\tnewick string = string is the Newick representation of a tree\n'+'MODEL = DL | DLT'

def getTree(s1,s2):
    if s1 == 'newick':
        tree = newick2Tree(s2)
        labelTree(tree)
        return(tree)
    
    if not s2.isdigit():
        print(USAGE)
        sys.exit()
    else:
        k = int(s2)                
        if s1 == 'random':
            tree = randomOrderedBinaryTree(k)
        elif s1 == 'caterpillar':
            tree = buildCaterpillar(k)
        elif s1 == 'complete':
            tree = buildCompleteTree(k)
        elif s1 == 'rrandom':
            tree_aux = randomOrderedBinaryTree(k)
            (ranking,tree) = rankTreeRandomly(tree_aux)
        elif s1 == 'rcaterpillar':
            tree_aux = buildCaterpillar(k)
            (ranking,tree) = rankTreeRandomly(tree_aux)
        elif s1 == 'rcomplete':
            tree_aux = buildCompleteTree(k)
            (ranking,tree) = rankTreeRandomly(tree_aux)
        else:
            tree = newick2Tree(s1)
        labelTree(tree)
        return(tree)

if __name__ == "__main__":
    tree_type = sys.argv[1]
    next_arg = 2
    tree = getTree(tree_type,sys.argv[next_arg])
    next_arg += 1
    stree = tree.asNewick()
    
    # Reading the evolutionary model
    M = sys.argv[next_arg]
    next_arg += 1
    if M == 'DL':
        MODEL={'D':True,'L':True,'T':False}
    elif M == 'DLT':
        MODEL={'D':True,'L':True,'T':True}
    else:
        print(USAGE)
        sys.exit()

    # Asymptotic growth factor
    (sysEqs,sysVars) = getSystem(tree,MODEL)
    solveSystem(sysEqs,sysVars)
    from tmp_solver_ import SOL
    print(stree+'\t'+str(SOL))
    os.system('rm tmp_solver_.py')
