"""
Transforms a Polyhedron by running k rounds of Sherali-Adams hierarchy. """
import numpy as np
import itertools
#from math import comb
#from scipy.special import comb
#from cvxopt import matrix
#import sage
#from sage.all import matrix
from functools import reduce
import operator as op

def comb(n, r):
    if r > n: return 0    
    else :
        r = min(r, n-r)
        numer = reduce(op.mul, range(n, n-r, -1), 1)
        denom = reduce(op.mul, range(1, r+1), 1)
        return numer / denom

def fill_look(k,n, look = {}):    
    sm = 0
    for i in range(k):
        lk = {}
        sm = sm + comb(n,i+1)
        for (u,v) in itertools.combinations(range(int(sm)),2):
            if(u < n and v > sm - comb(n,i+1) and len(set(flatten(look,(u,v),i+1,n))) >= i+2):
                lk[extend_index_mem(look,(u,v),n,i+1)] = (u,v)
        look[i+1] = lk
    return look

def flatten(look, pair,k,n,acc = []):
    (a,b) = pair
    if a < b:
        a,b = b,a
    if (a < n and b < n) :
        return acc + [a,b]
    else :
        (bb,aa) = look[k-1][a]
        if b == bb:
            return flatten(look,(aa,bb),k-1,n,acc)
        else: 
            return flatten(look,(aa,bb),k-1,n, acc + [b])

def flatten2(pair,k,n):
    (a,b) = pair
    if a < b:
        a,b = b,a
    if (a < n and b < n) :
        return [a,b]
    else :
        return list(set(invert(a,n,k) + [b]))
        

def compute_index(li,n,k):
    sm = 0
    li.sort()
    for i in range(len(li)):
        sm = sm + int(comb(li[i],i+1))
    for i in range(1, len(li)):
        sm = sm + int(comb(n,i))
    return sm  

def sum_comb(k,n):
    return sum(map(lambda i: int(comb(n,i)), range(1,k+1)))

def invert(ix,n,k):
    bnd = sum_comb(k,n)
    if bnd > ix: return invert(ix,n,k-1)
    return [list(li) for li in itertools.combinations(range(n),k+1) if compute_index(list(li),n,k) == ix][0]
    
    
def extend_index(pair,n,k):
    """
    computes a linearized index from set indexed monomial y_{S}
    INPUT:
    pair -- (a,b) to map
    n -- size of N
    k -- round of SA
    """
    (a,b) = pair
    
    if(a == b):
        return a
    if(b > a):
        a,b = b,a
    if(k == 1):
        return a*(a-1)/2 + b + n 
    if(a < n):
        return extend_index(pair,n, 1)
    else:
        li = flatten2(pair,k,n)
        return int(compute_index(li,n,k))

def extend_index_mem(look,pair,n,k):
    """
    memoized version. must call fill_look first
    computes a linearized index from set indexed monomial y_{S}
    INPUT:
    pair -- (a,b) to map
    n -- size of N
    k -- round of SA
    """
    (a,b) = pair
    
    if(a == b):
        return a
    if(b > a):
        a,b = b,a
    if(k == 1):
        return a*(a-1)/2 + b + n 
    if(a < n):
        return extend_index_mem(look,pair,n, 1)
    else:
        li = list(set(flatten(look,pair,k,n)))
        return int(compute_index(li,n,k))

#must pass previous round instance as arguments
# input c,A,b of the form min cx s.t. Ax <= b
def get_SA_instance(rnd,n,oldA,oldb, look = {}):
    (rows,cols) = oldA.shape
    oldN = cols
    newN = sum(map(lambda j:int(comb(n,j)),range(1,rnd+2)))    
    newA = np.array(np.hstack([oldA, np.zeros((rows,newN - cols))]).tolist())

    acc = []
    bacc = []
    for k in range(n):
        Ak = np.zeros((rows,newN))  #y
        nAk = np.zeros((rows,newN)) #1-y
        newB = np.zeros(rows)
        nnewB = np.zeros(rows)
        for r in range(rows):
            for c in range(cols):
                j = int(extend_index((c,k),n,rnd)) if not look else int(extend_index_mem(look,(c,k),n,rnd))
                if(oldA[r,c] != 0):                
                    Ak[r,j] += oldA[r,c]
                    nAk[r,c] += oldA[r,c]
                    nAk[r,j] += -oldA[r,c]            
            Ak[r,k] += -oldb[r] 
            nAk[r,k] += oldb[r]
            nnewB[r] = oldb[r]            
        acc.append(Ak)
        bacc.append(newB)
        acc.append(nAk)
        bacc.append(nnewB)
        
    xxx_min = 0 * np.ones(newN)
    xxx_max = 1 * np.ones(newN)
    GG3 = np.vstack([
        np.hstack([+np.eye(newN)]),
        np.hstack([-np.eye(newN)])])
    hh3 = np.hstack([xxx_max, -xxx_min])
    
    SA = np.array(np.asarray(list(itertools.chain(*acc))))
    Sb = np.asarray(list(itertools.chain(*bacc)))

    return (np.array(np.vstack([newA,SA,GG3])), np.hstack([oldb,Sb,hh3]))

def run_SA(k,n,A,b,memoize = False):
    look = {}
    if(memoize):
        look = fill_look(k,n)
    tA = A
    tb = b
    for i in range(1,k+1):
        (tA,tb) = get_SA_instance(i,n,tA,tb,look)
    return (tA,tb)
    
