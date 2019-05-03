__author__  = "Cedric Chauve, Yann Ponty, Michael Wallner"
__date__    = "April 24, 2019"

# Class and function to manipulate rooted trees

import sys
import random

# -----------------------------------------------------------------------------
# Class defining a rooted tree, erepresented by a node, the root of the tree
class Node:
    def __init__(self, children=[], id=-1, time=-1):
        self.children = children[:]
        for c in self.children:
            c.setParent(self)
        self.id = id
        self.parent = None
        self.time = time  # For ranked trees

    def getID(self):
        return(self.id)
    def setID(self,id):
        self.id =  id
    def getTime(self):
        return(self.time)
    def setTime(self,time):
        self.time =  time
        
    def getHeight(self): # Height = length of the longest path leaf-root, in number of edges
        h = 0
        for c in self.children:
            hc = 1+c.getHeight()
            if hc > h:
                h = hc
        return(h)
    
    def getLength(self): # Length = number of nodes.
        n = 1
        for c in self.children:
            n += c.getLength()
        return(n)

    def getSize(self): # Size = number of leaves
        if len(self.children)==0:
            n=1
        else:
            n=0
            for c in self.children:
                n += c.getSize()
        return(n)
            

    def getSize(self): # Size = number of leaves
        size = 0
        for c in self.children:
            size += c.getSize()
        return(max(1,size))
        
    def setParent(self,u):
        self.parent = u
    def getParent(self):
        return(self.parent)
    def getRoot(self):
        if self.parent is not None:
            return(self.parent.getRoot())
        else:
            return(self)

    def isRoot(self):
        return(self.parent==None)
    def isBinary(self):
        return(len(self.children)==2)
    def isUnary(self):
        return(len(self.children)==1)
    def isLeaf(self):
        return(len(self.children)==0)
        
    def getChildren(self):
        return(self.children)
    def setChildren(self,children):
        self.children = children
    def getLeft(self):
        assert len(self.children)==2
        return(self.children[0])
    def getRight(self):
        assert len(self.children)==2
        return(self.children[1])
    def getChild(self):
        assert len(self.children)==1
        return(self.children[0])

    def allNodes(self):
        res = []
        for c in self.children:
            res += c.allNodes()
        res += [self]
        return(res)

    def getSiblings(self,u):
        return([v for v in self.children if u!=v])

    def getIncomparable(self):
        res = []
        if not self.isRoot():
            p = self.getParent()
            for c in p.getSiblings(self):
                res += c.allNodes()
            res += p.getIncomparable()
        return(res)

    # Get the nodes in a time slice if the tree is ranked
    def getContemporary(self):
        return([v for v in self.getIncomparable() if v.getTime() == self.time])

    # Remove leaves of a tree
    def shave(self):
        if len(self.children)>0:
            nchildren = [c.shave() for c in self.children if c.shave() is not None]
            return(Node(nchildren, self.id, self.time))
        else:
            return(None)

    # Create a copy of the tree
    def copy(self):
        copyChildren = []
        for child in self.getChildren():
            copyChildren.append(child.copy())
        newTree = Node(copyChildren,self.getID(),self.getTime())
        return(newTree)

    # Build the unique unary-binary tree for a ranked tree
    def buildUnaryBinary(self,currentTime=0):
        t = self.getTime()
        if t>currentTime:
            return(Node([self.buildUnaryBinary(currentTime+1)],-1,currentTime))
        else:
            nChildren = []
            for child in self.getChildren():
                nChildren.append(child.buildUnaryBinary(currentTime+1))
            return(Node(nChildren,self.getID(),currentTime))

    # Representation as a string
    #def asString(self):
    #    return "("+",".join([t.asString() for t in self.children])+")"

    # Representation as a string in Newick format
    def asNewick(self,date=False):
        newick = ""
        children   = self.getChildren()
        nbChildren = len(children)
        if nbChildren > 0:
            newick += "("
        i = 0
        for child in children:
            i += 1
            newick += child.asNewick()
            if i < nbChildren:
                newick += ","
        if nbChildren > 0:
            newick += ")"
        newick += str(self.getID())
        if date:
            newick += ":"+str(self.getTime())
        return(newick)

    # Auxiliary functions to build a random binary tree
    def extendLeft(self):
        extensionR = Node()
        extensionL = Node(self.children)
        extensionR.setParent(self)
        extensionL.setParent(self)
        self.children = [extensionL,extensionR]
        return(self.children)
    def extendRight(self):
        extensionL = Node()
        extensionR = Node(self.children)
        extensionR.setParent(self)
        extensionL.setParent(self)
        self.children = [extensionL,extensionR]
        return(self.children)
# -----------------------------------------------------------------------------

# Label the nodes of a tree in post-order traversal
# Returns the label of the root
def labelTree(t,currentID=0):
    for c in t.getChildren():
        currentID = labelTree(c,currentID)
    t.setID(currentID)
    return(currentID+1)

def printTree(t,indent=0):
    print(" "*(indent)+"-> %s (t=%s)"%(t.getID(),t.getTime()))
    for c in t.getChildren():
        printTree(c,indent+1)

# Check if a tree is balanced
def checkBalanced(t): 
    children = t.getChildren()
    if len(children) == 0:
        return(True)
    else:
        balanced = True
        minSize = -1
        maxSize = -1
        for child in children:               
            balanced = balanced and checkBalanced(child)
            sizeChild = child.getSize()
            if minSize == -1 or minSize > sizeChild:
                minSize = sizeChild
            if maxSize < sizeChild:
                maxSize = sizeChild
        balanced = balanced and abs(minSize-maxSize)<=1
        return(balanced)

# Check if a tree is a caterpillar
def checkCaterpillar(t):
    children = t.getChildren()
    n = t.getLength()
    if len(children) == 0:
        return(True)
    elif len(children) != 2:
        return(False)
    else:
        sizeLeft  = t.getLeft().getLength()
        sizeRight = t.getRight().getLength()
        return(((sizeLeft==1 and sizeRight==n-2) or (sizeLeft==n-2 and sizeRight==1)) and checkCaterpillar(t.getLeft()) and checkCaterpillar(t.getRight()))


    # Auxiliary function to rank randomly a tree
def correctedLength(t):
    m = t.getLength()
    return(m-((m+1)/2))

# Check if a tree is balanced
def checkBalanced(t): 
    children = t.getChildren()
    if len(children) == 0:
        return(True)
    else:
        balanced = True
        minSize = -1
        maxSize = -1
        for child in children:               
            balanced = balanced and checkBalanced(child)
            sizeChild = child.getSize()
            if minSize == -1 or minSize > sizeChild:
                minSize = sizeChild
            if maxSize < sizeChild:
                maxSize = sizeChild
        balanced = balanced and abs(minSize-maxSize)<=1
        return(balanced)

# Check if a tree is a caterpillar
def checkCaterpillar(t):
    children = t.getChildren()
    n = t.getLength()
    if len(children) == 0:
        return(True)
    elif len(children) != 2:
        return(False)
    else:
        sizeLeft  = t.getLeft().getLength()
        sizeRight = t.getRight().getLength()
        return(((sizeLeft==1 and sizeRight==n-2) or (sizeLeft==n-2 and sizeRight==1)) and checkCaterpillar(t.getLeft()) and checkCaterpillar(t.getRight()))
    
    
# Random ranking of a tree: returns the ranking and the unary-binary trees given by the ranking
def rankTreeRandomly(t):
    newTree = t.copy()
    # Generating a ranking for the internal nodes of a copy of t
    res     = []
    l       = correctedLength(newTree)
    weights = {newTree:l}
    leaves  = []
    while len(res)<l:
        acc   = sum(weights.values())
        r     = random.random()*acc
        nodes = weights.keys()
        for u in nodes:
            w = weights[u]
            r -= w
            if r<0:
                res.append(u)
                del weights[u]
                for v in u.getChildren():
                    if v.isLeaf():
                        leaves.append(v)
                    else:
                        weights[v] = correctedLength(v)
                break
    # Setting the times/ranks for the nodes and leaves of the copy of t
    ranks = {}
    for i,v in enumerate(res):
        v.setTime(i)
        ranks[v.getID()] = i
    for v in leaves:
        v.setTime(len(res))
        #ranks[v.getID()] = len(res)
    # Inserting unary nodes on the branches of the ranked tree
    newTree = newTree.buildUnaryBinary()
    # Labeling nodes in postorder
    labelTree(newTree)
    return((ranks,newTree))

# Print a ranking
def printRanking(ranking):
    rankingStr = ""
    for r in ranking.keys():
        rankingStr += "("+str(r)+","+str(ranking[r])+"),"
    rankingStr +=","
    return(rankingStr.replace(",,",""))

# Read a ranking
def readRanking(s):
    ranking={}
    s1 = s.rstrip().split('),(')
    for r in s1:
        r1=r.replace('(','').replace(')','').split(',')
        ranking[int(r1[0])]=int(r1[1])
    return(ranking)

# Rank a tree given a ranking
def rankTree(tree,ranking):
    newTree = tree.copy()
    n = newTree.getSize()-1 # Number of internal nodes
    allNodes = newTree.allNodes()
    for v in allNodes:
        if not v.isLeaf():
            v.setTime(ranking[v.getID()])
        else:
            v.setTime(n)
    # Inserting unary nodes on the branches of the ranked tree
    newTree = newTree.buildUnaryBinary()
    # Labeling nodes in postorder
    labelTree(newTree)
    return(newTree)
    
# Build a caterpillar tree with k leaves
def buildCaterpillar(k):
    if k == 1:
        return(Node())
    else:
        return(Node([Node(),buildCaterpillar(k-1)]))

# Build a complete binary tree of depth h (2^h leaves)
def buildCompleteTree(h):
    if h == 0:
        return(Node())
    else:
        return(Node([buildCompleteTree(h-1),buildCompleteTree(h-1)]))

# Build a balanced binary tree with k leaves
def buildBalancedTree(k):
    if k == 1:
        return(Node())
    else:
        k1 = int(k/2)
        k2 = k-k1
        return(Node([buildBalancedTree(k1),buildBalancedTree(k2)]))
    
# Build a random ordered binary tree with k leaves; Remy's algorithm
def randomOrderedBinaryTree(k):
    nodes = [Node()]
    while(len(nodes))<2*k-1:
        v = random.choice(nodes)
        if random.randint(0,1)==0:
            nodes += v.extendLeft()
        else:
            nodes += v.extendRight()
    nroot = nodes[0].getRoot()
    return(nroot)

# Build a random *unordered* binary tree with k leaves; modified RANRUT algorithm
def precomputeBinaryTrees(n):
    A = [0 for i in range(0,n+1)]
    A[1] = 1
    for i in range(2,n+1):
        for m in range(1,i):
            A[i] += A[m]*(i-1-m)*A[i-1-m]
        if (i-1)%2 == 0:
            m = (i-1)/2
            A[i] += m*A[m]
        A[i] /= (i-1)
    return(A)
def randomBinaryTree_aux(n,A):
    res = []
    if n>1:
        r = random.random()*A[n]*(n-1)
        m = 1
        for m in range(1,n):
            r -= A[m]*(n-m-1)*A[n-1-m]
            if r<0:
                leftTree  = randomBinaryTree_aux(m,A)
                rightTree = randomBinaryTree_aux(n-1-m,A)
                res = [leftTree,rightTree]
                break
        if r>=0:
            leftTree  = randomBinaryTree_aux((n-1)/2,A)
            rightTree = leftTree.copy()
            res = [leftTree,rightTree]
    return(Node(res))
def randomBinaryTree(k):
    n=2*k-1
    A = precomputeBinaryTrees(n)
    return(randomBinaryTree_aux(n,A))

def newick2Tree(s):
    # Assumption: s is the newick string of a rooted binary tree
    if s[0] != "(": # Case 1: tree reduced to a leaf
        s1 = s.split(':')[0]
        return(Node([],s1,-1))
    else: # Case 2, tree with a root
        # Assumption: we are on the opening (
        i = 1 # Next character
        # Looking for the left subtree
        j = i;
        if s[j] == "(":
            # We look for the closing parenthesis
            c = 1
            while c != 0:
                j += 1
                if s[j] == "(":
                    c += 1
                elif s[j] == ")":
                    c -= 1
            s1 = s[i:j+1]
            while s[j] != ",":
                j += 1
        else:
            while s[j] != ",":
                j += 1
            s1 = s[i:j]
        # Looking for the right subtree
        i = j+1
        j = i
        if s[j] == "(":
            # We look for the closing parenthesis
            c = 1
            while c != 0:
                j += 1
                if s[j] == "(":
                    c += 1
                elif s[j] == ")":
                    c -= 1
            s2 = s[i:j+1]
        else:
            while s[j] != ")":
                j += 1
            s2 = s[i:j]
        # We build the tree
        return(Node([newick2Tree(s1),newick2Tree(s2)]))

# -----------------------------------------------------------------------------------------
# Default: generate a random tree with k print it

if __name__ == "__main__":
    k = int(sys.argv[1])
    t = randomBinaryTree(k)
    labelTree(t)
    s = t.asNewick()
    print(s)

# -----------------------------------------------------------------------------------------
# Unused

# # Queue of all nodes of a tree in Breadth-First Search order
# def BFSQueue(t):
#     allNodes   = t.allNodes()
#     nodesQueueAux = []
#     nodesQueueAux.append(t.getRoot())
#     nodesQueue = []
#     while len(nodesQueueAux)>0:
#         v = nodesQueueAux[0]
#         for c in v.getChildren():
#             nodesQueueAux.append(c)
#         nodesQueueAux.remove(v)
#         nodesQueue.append(v)
#     return(list(reversed(nodesQueue)))

# # Compute a list of unique triples (vLabel,c1Label,c2Label) defined recursively as follows:
# # leaves are labeled 0
# # an internal node label is a integer uniquely defined by the unordered pair of the labels of its two children
# # so if two internal nodes have pairs of children with the same labels, they both receive the same label
# # returns the triples (node,label child1, label child2) with the children labels in lexicographic order
# # with each triple, indicating a subtree, appears once.
# def computeUniqueSubtrees(t):
#     allNodes      = t.allNodes()
#     indexSubtrees = {v: -1 for v in allNodes} # index for the subtree rooted at each node
#     labelSubtrees           = {}  # indexed by pairs (leftLabel,rightLabel) and pointing to the subtree rootLabel
#     currentLabel            = 0
#     labelSubtrees[(-1,-1) ] = currentLabel # leaves
#     # Height of leaves
#     for v in allNodes:
#         if v.isLeaf():
#             indexSubtrees[v] = currentLabel
#     # Computing the bottom-up traversal queue of the nodes
#     nodesQueue  = BFSQueue(t)
#     # Bottom-up traversal to compute the index of each node
#     for v in nodesQueue:
#         if not v.isLeaf():
#             leftLabel  = indexSubtrees[v.getLeft()]
#             rightLabel = indexSubtrees[v.getRight()]
#             if (leftLabel,rightLabel) in labelSubtrees.keys():
#                 indexSubtrees[v] = labelSubtrees[(leftLabel,rightLabel)]
#             else:
#                 currentLabel += 1
#                 indexSubtrees[v] = currentLabel
#                 labelSubtrees[(leftLabel,rightLabel)] = currentLabel
#                 labelSubtrees[(rightLabel,leftLabel)] = currentLabel
#     # Extracting the unique triples (labelNode,labelChild1,labelChild2)
#     seenLabels  = []
#     triplesList = []
#     for v in allNodes:
#         vLabel = indexSubtrees[v]
#         if not (v.isLeaf() or vLabel in seenLabels):
#             leftLabel  = indexSubtrees[v.getLeft()]
#             rightLabel = indexSubtrees[v.getRight()]
#             if leftLabel>rightLabel:
#                 triplesList.append((vLabel,rightLabel,leftLabel))
#             else:
#                 triplesList.append((vLabel,leftLabel,rightLabel))
#             seenLabels.append(vLabel)
#
#     return(triplesList)
