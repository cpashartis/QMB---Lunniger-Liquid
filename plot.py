# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 17:26:54 2015

@author: christoforos
Script for running plots of diagonlization code

Currently configured for:
Plotting Eigen energies of the diagonalization
"""

from matplotlib import use       
use('agg')
from matplotlib import rc_file
rc_file('/Users/christoforos/bin/rcplots/original.rc')
import matplotlib.pyplot as plt
from numpy import *
from math import factorial
from Diagonalization import which
from mcint import integrate

#system vars
N = 2
L = 1
hbar = 1
m=.5

def nextstate(currentstate,nc):  
    '''Compute the next occupancy state from the current one'''
    
    copy(currentstate)
    k = max(nonzero(currentstate[0:2*nc])[0])
    currentstate[k] = currentstate[k]-1
    currentstate[k+1] = N-sum(currentstate[0:k+1])
    for i in range(k+2,2*nc+1):
        currentstate[i] = 0 
    return currentstate
    
def findstates(N,nc):
    """Fixed to find states inside a function, because used in real space"""
    num_states = math.factorial(2*nc+N)/(math.factorial(2*nc)*math.factorial(N))
    startstate = append([N],[0]*2*nc)
    laststate = append([0]*2*nc,[N])
    states = empty([num_states,2*nc+1],dtype=int)
    counter = 0
    
    while array_equal(startstate,laststate) == False:
        states[counter] = startstate
        startstate = nextstate(startstate,nc)
        counter += 1
    states[-1] = laststate
    
    return states
    
def RealEigenfunction(second_quant, L, N, nc, variables):
    """This function calculates the realspace wavevector given an N particle system
    and its wavevector in second quantization"""
    
    #First, transform second quantization to momentum space by symmetrizing
    #there are N! different combinations
    #num_states = math.factorial(2*nc+N)/(math.factorial(2*nc)*math.factorial(N))
    generated_states = findstates(N,nc)
    
    #for now two state
    temp =''
    temp2 = '' #for conjugate
    state_counter = 0
    N_2_coef = []
    for state in generated_states:
        #check for l, -l
    #####################
        if N ==2:
            print nonzero(state)
            if nonzero(state)[0][0] == 0 and len(nonzero(state)[0]) == 1:
                print 0
                N_2_coef.append(second_quant[state_counter])
            elif len(nonzero(state)[0])==2:
                if nonzero(state)[0][0] == 1 and nonzero(state)[0][1]==2:
                    print 1
                    N_2_coef.append(second_quant[state_counter])                
    #############################
        temp += '(%.10f)*(' %(second_quant[state_counter])
        temp2 += '(%.10f)*(' %(second_quant[state_counter])
        variable_counter = 0
        for n in nonzero(state)[0]: #multiply by exponential for each mode
            for total_n in range(state[n]): #how many particles per mode
                if n%2==0: #handle negative stuff
                    j=-1
                else:
                    j=1
                temp += 'exp(-1j*%i*pi/%i*%s)*' %( (n+1)/2 *j , L, variables[variable_counter])
                temp2 += 'exp(1j*%i*pi/%i*%s)*' %( (n+1)/2 *j, L, variables[variable_counter])
                variable_counter +=1
        
        temp = temp.rstrip('*')
        temp2 = temp2.rstrip('*')
        #now symmetrize
        if state[n] != N: #check to see if all particles in one state, if not continue
            temp+='+' #symmetrize first half
            temp2+='+' #symmetrize first half
            variable_counter = len(variables) -1 #reverse order for symmetry
            for n in nonzero(state)[0]: #multiply by exponential for each mode
                for total_n in range(state[n]): #how many particles per mode
                    if n%2 ==0: #handle negative stuff
                        j=-1
                    else:
                        j=1
                    temp += 'exp(-1j*%i*pi/%i*%s)*' %( (n+1)/2 *j , L, variables[variable_counter])
                    temp2 += 'exp(1j*%i*pi/%i*%s)*' %( (n+1)/2 *j, L, variables[variable_counter])
                    variable_counter -=1
            temp = temp.rstrip('*')
            temp2 = temp2.rstrip('*')   
            temp +=')*1/sqrt(factorial(%i))'%(N)
            temp2 +=')*1/sqrt(factorial(%i))'%(N)
        else:
            temp += ')'
            temp2 += ')'
            
        temp +='+'
        temp2 +='+'
        state_counter +=1
    temp = temp.rstrip('+') #remove spare + should there be one
    temp2 = temp2.rstrip('+')
    
#    print temp +'\n'
#    print temp2

    return eval('lambda x,y: ' + temp), eval('lambda x,y: ' + temp2), N_2_coef
    
def sampler(L):
    """Give it an array L with the upper interval of integration"""
    while True:
        yield tuple(random.uniform(0,L)) #returns tuple of random numbers of upper range
        

if __name__ == '__main__': #if script run by itself, not imported
    #define 1st order perturbation theory, cuttoff = 1
    E1= lambda g : g/(L) - m*g**2/(hbar**2*pi**2)
    f1 = plt.figure()
    ax1 = f1.add_subplot(111)
    f2 = plt.figure()
    ax2 = f2.add_subplot(111)
    
    #define exact solutions for N=2
    #renyi_prob = 
    for nc in which(N):
        eigvecs = load('N%dM%deigvecs.npy' % (N,nc))
        eigvals = load('N%dM%deigvals.npy' % (N,nc))
        g=linspace(.1,5.0,len(eigvals))
        #ax1.plot(g, a[0]*g**2 + a[1]*g + a[2], 'k' )
        
        ax2.plot(g,E1(g), 'r--', label = 'perturbation')
        ax2.legend(loc = 'lower right')
        ax2.set_xlabel('Coupling Constant, g')
        ax2.set_ylabel('Eigen-Energy')    
        
        #test real stuff
        #RealEigenfunction(eigvecs,L,N)
        #phi_x1y, phi_c_x1y = RealEigenfunction(eigvecs[:,0], L, N, nc,variables = ['x1','y'])
        #phi_xy1, phi_c_x1y = RealEigenfunction(eigvecs[:,0], L, N, nc,variables = ['x','y1'])
        
        for g_i in array([.1,1.0,5.0]):
            phi, phi_c,coefs = RealEigenfunction(eigvecs[:,where(g == g_i)], L, N, nc,variables = ['x','y'])
            #print coefs, '%.10f' %(coefs[0]**2+coefs[1]**2)
            phi_check = lambda x,y: coefs[0] + coefs[1]*1/sqrt(factorial(2))*( exp(-1j*x*pi)*exp(1j*y*pi)+
                exp(-1j*y*pi)*exp(1j*x*pi))
            integrand0 = lambda x: phi(x[0],x[1])*phi_c(x[2],x[3])*phi(x[2],x[3])*phi_c(x[0],x[1])
            integrand0_check = lambda x:phi_check(x[0],x[1])*phi_check(x[2],x[3])\
                *phi_check(x[2],x[3])*phi_check(x[0],x[1])
            integrand1 = lambda x: phi(x[0],x[1])*phi_c(x[0],x[1])*phi(x[2],x[3])*phi_c(x[2],x[3])
            integrand1_check = lambda x: phi_check(x[0],x[1])*phi_check(x[0],x[1])*phi_check(x[2],x[3])*phi_check(x[2],x[3])
            integrand2 = lambda x: phi(x[0],x[1]) * phi_c(x[2],x[1]) * phi(x[2],x[3]) * phi_c(x[0],x[3])
            integrand2_check = lambda x: phi_check(x[0],x[1]) * phi_check(x[2],x[1]) * phi_check(x[2],x[3]) * phi_check(x[0],x[3])
            
            final = []
            A = linspace(0.01,.5,51)
            B = L-A
            for j in A:
                print j
                k = L-j
                #result0_check,error0_check = integrate(integrand0_check, sampler([j,j,j,j]),measure = j**4, n=100)
                result0, error0 = integrate(integrand0, sampler([j,j,j,j]),measure = j**4, n=1000)
                #result1_check,error0_check = integrate(integrand1_check, sampler([k,k,k,k]),measure = k**4, n=100)
                result1, error1 = integrate(integrand1, sampler([k,k,k,k]),measure = k**4, n=1000)
                #result2_check,error0_check = integrate(integrand2_check, sampler([j,k,j,k]), measure = j**2*k**2, n=100)
                result2, error2 = integrate(integrand2, sampler([j,k,j,k]), measure = j**2*k**2, n=1000)
                final.append(result0+result1+2*result2)
                #print final, result0, result1, result2
                #print result0_check, result1_check, result2_check
            ax1.plot(A,-1*log(final), label = '%.2g' %(g_i))
        ax1.legend()
        ax1.set_xlabel('Region A (length)')
        ax1.set_ylabel('Renyi Entropy ' + r'$\alpha = 2$')
        f1.savefig('Renyi Entanglement Entropy %i' %(nc))
        ax1.clear()
        print 'done %i' %(nc)
        
    f1.savefig('Eigen-Energies in K space')
    plt.clf()
    
    
    


#load in eigen energies and eigenvalues