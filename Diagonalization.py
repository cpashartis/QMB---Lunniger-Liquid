#Steps:
#1) Build matrices and save to disk using Build_Matrices.py
#2) Diagonalize to get eigenvectors using this code
#3) Make plots using Value_Plots through the eigenvectors.
#Note. This is the memory-consuming process. 

import numpy 
import scipy as sp
from scipy.sparse import *
from scipy import *
from scipy.sparse.linalg import eigsh
import scipy.io
import timeit
import math

#starttime=timeit.default_timer()
#h = hpy()

N=2
gmax=5.0  #Maximum interaction strength
dg=0.1    #Interval between interaction strengths
#Number of values used will be gmax/dg
    
def which(N):
    '''Define lists of accessible cutoffs, for each N.'''
    
    if N==2:
        ncvec=[1,2,5,10]
    elif N==4:
        ncvec=[1,2,3,4,5,6,7,8]
    elif N==8:
        ncvec=[1,2,3,4]
    elif N==16:
        ncvec=[1,2,3]
    elif N==32:
        ncvec=[1,2]
    elif N==64:
        ncvec=[1]
        
    return ncvec

if __name__ == '__main__': #if script run by itself, not imported
    ncvec=which(N)
    
    for m in ncvec:
        
        #starttime = timeit.default_timer()
        num_states=math.factorial(2*m+N)/(math.factorial(2*m)*math.factorial(N))
        T=scipy.io.mmread('N%dM%dT.mtx' % (N,m)).tocsr()
        V=scipy.io.mmread('N%dM%dV.mtx' % (N,m)).tocsr()
        eigvecs=numpy.empty([num_states,gmax/dg])
        eigvals = numpy.empty(gmax/dg)
        gvals=numpy.arange(0.1,gmax+0.1,dg)
        
        #Diagonalize Hamiltonian (T+gV) for each cutoff in ncvec.
        for i,j in enumerate(gvals):
            eigvals[i], temp =eigsh(T+j*V,1,sigma=0,which='LA')
            eigvecs[:,i] = temp[:,0] #final splice for indexing direction
        #Save the eigenstates named to disk. File named e.g. N4M2eigvecs.npy.
        if N>=2:
            numpy.save('N%dM%deigvecs' % (N,m),eigvecs)
            numpy.save('N%dM%deigvals' % (N,m),eigvals)
    
        #stoptime = timeit.default_timer()
        #print stoptime-starttime
    
    #stoptime=timeit.default_timer()
    #print stoptime-starttime
