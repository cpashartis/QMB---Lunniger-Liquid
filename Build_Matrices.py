#This script will construct the kinetic and potential energy matrices
#along with the number matrix and the projection matrices, and save them all
#to disk for future use.

import math
import numpy as np
from scipy.sparse import *
from scipy import *
from scipy.sparse.linalg import eigsh
import scipy.io
    
def KD(a,b):  #Kronecker delta
    if a == b:
        return 1
    else:
        return 0
        
#starttime=timeit.default_timer()

#Define parameters

nc = 8   # mode cut-off level
N = 2   # number of particles
L = 1    # Ring circumference
g = 0.1  # interaction strength
gmin = g
gmax = 5
dg = 0.1
numg = 100
pi = math.pi
hbar = 1.0
mass = 0.5
lam = hbar**2/(2*mass)
alpha = 1 # renyi order
num_states = math.factorial(2*nc+N)/(math.factorial(2*nc)*math.factorial(N))

def nextstate(currentstate,nc):  
    '''Compute the next occupancy state from the current one'''
    
    np.copy(currentstate)
    k = np.max(np.nonzero(currentstate[0:2*nc])[0])
    currentstate[k] = currentstate[k]-1
    currentstate[k+1] = N-np.sum(currentstate[0:k+1])
    for i in range(k+2,2*nc+1):
        currentstate[i] = 0 
    return currentstate
    
def findstates(N,nc):
    """Fixed to find states inside a function, because used in real space"""
    startstate = np.append([N],[0]*2*nc)
    laststate = np.append([0]*2*nc,[N])
    states = np.empty([num_states,2*nc+1],dtype=int)
    counter = 0

    while np.array_equal(startstate,laststate) == False:
        states[counter] = startstate
        startstate = nextstate(startstate,nc)
        counter += 1
    states[-1] = laststate
    
    return states

def c(vec,m):
    '''Creation operator'''
    
    cvecstart = np.copy(vec)
    if m < -nc or m > nc:
        cval = 0
    elif m >= -nc and m <= nc:
        if m > 0:
            if vec[2*m-1]+1 < 0:
                cval = 0
            else:
                cval = np.sqrt((vec[2*m-1])+1)
            cvecstart[2*m-1] = cvecstart[2*m-1]+1
        else:
            if vec[-2*m]+1 < 0:
                cval = 0
            else:
                cval = np.sqrt(vec[-2*m]+1)
            cvecstart[-2*m] = cvecstart[-2*m]+1
    
    return cval,cvecstart

def a(vec,m):
    '''Annihilation operator'''
    
    avecstart = np.copy(vec)
    if m < -nc or m > nc:
        aval = 0
    elif m >= -nc and m <= nc:
        if m > 0:
            if vec[2*m-1] < 0:
                aval = 0 
            else:
                aval = np.sqrt(vec[2*m-1])
            avecstart[2*m-1] = avecstart[2*m-1]-1
        else:
            if vec[-2*m] < 0:
                aval = 0
            else:
                aval = np.sqrt(vec[-2*m])
            avecstart[-2*m] = avecstart[-2*m]-1

    return aval, avecstart
                    
def Vop(s,m,mp,mdp):
    '''Compute matrix element of the two body operator'''
    
    step1 = a(HS[s],m)
    step2 = a(step1[1],mp)
    step3 = c(step2[1],mdp)
    step4 = c(step3[1],m+mp-mdp)
    inner = 1.0/(L*np.sqrt(1+KD(m,mp))*np.sqrt(1+KD(mdp,m+mp-mdp)))*step1[0]*step2[0]*step3[0]*step4[0]
    
    if inner != 0:
        return inner, HS_dict[tuple(step4[1])]

def get_indices(i): 
    '''Find the modes that give non-zero matrix elements'''
    
    index_list = np.nonzero(HS[i])[0]
    for m in range(len(index_list)):
        if index_list[m] % 2 == 0:
            index_list[m] = -index_list[m]/2
        else: 
            index_list[m] = (index_list[m]+1)/2
    return sorted(index_list)

def Build():
    '''Construct T, V, and n matrices'''
    
    vecstart=[]
    rowstart=[]
    colstart=[]
    vec2 = np.empty(num_states)
    vec3 = np.empty(num_states,dtype=int)
    for j in range(0,num_states):
        
        indices=get_indices(j)
        
        KEvals=np.sum([HS[j,2*s-1]*2*pi**2*hbar**2*s**2/(mass*L**2)+HS[j,2*s]*2*pi**2*hbar**2*s**2/(mass*L**2) for s in range(1,nc+1)])

        vec2[j]=KEvals
        vec3[j]=HS[j,0]
        
        for k in indices:
            for l in indices: 
                for p in range(-nc,nc+1):    
                    if -nc<=k+l-p<=nc:
                        try:
                            if p>k+l-p and k>l:
                                vecstart.append(4*Vop(j,k,l,p)[0])
                                rowstart.append(Vop(j,k,l,p)[1])
                                colstart.append(j)
                            elif p==k+l-p and k>l:
                                vecstart.append(2*Vop(j,k,l,p)[0])
                                rowstart.append(Vop(j,k,l,p)[1])
                                colstart.append(j)
                            elif p>k+l-p and k==l:
                                vecstart.append(2*Vop(j,k,l,p)[0])
                                rowstart.append(Vop(j,k,l,p)[1])
                                colstart.append(j)
                            elif l==p and l==k:
                                vecstart.append(Vop(j,k,l,p)[0])
                                rowstart.append(Vop(j,k,l,p)[1])
                                colstart.append(j)
                    
                        except TypeError:
                            pass       
        KEMatrix=scipy.sparse.spdiags(vec2,0,num_states,num_states)
        nMatrix=scipy.sparse.spdiags(vec3,0,num_states,num_states)
    return  KEMatrix, nMatrix, csc_matrix((vecstart,(rowstart,colstart)),shape=(num_states,num_states))

if __name__ == '__main__':
    HS = findstates(N,nc)
    print HS
    HS_dict = {tuple(a):i for i,a in enumerate(HS.tolist())}
        
    deltaliststart = []
    for i in range(N+1):
        deltastart=np.empty(num_states,dtype=int)
        for j in range(num_states):
            if HS[j,0]==i:
                deltastart[j]=1
            else:
                deltastart[j]=0
        delta=scipy.sparse.spdiags(deltastart,0,num_states,num_states)
        deltaliststart.append(delta)
    D=deltaliststart #Projection matrices
    
    #starttime1=timeit.default_timer()
    
    #could put all three together to speed up time
    T,n,V=Build() #condensate occupancy matrix
    #T=Build()[0] #Kinetic energy matrix
    #V=Build()[2] #potential energy matrix
    ######################################
    scipy.io.mmwrite('N%dM%dn.mtx' % (N,nc),n)
    scipy.io.mmwrite('N%dM%dT.mtx' % (N,nc),T)
    scipy.io.mmwrite('N%dM%dV.mtx' % (N,nc),V)
    for i in range(N+1):
        scipy.io.mmwrite('N%dM%dD%d.mtx' % (N,nc,i),D[i])
        
    #H=T+V
    #stoptime1=timeit.default_timer()
    #print stoptime1-starttime1
