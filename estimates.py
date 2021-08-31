import numpy as np
import math
import scipy.special


def up_burst(a,k,p,m):
    pro = 0
    d = int(a/2) + 1
    
    for i in range(d, 2*k+1):
        for i1 in range(d, min(a,i)+1):
            pro += pow(p, i)*pow(1-p, 2*k-i)*scipy.special.binom(a, i1)*scipy.special.binom(2*k-a, i-i1)
    if(a%2==0):
        d = int(a/2)
        for i in range(d, 2*k+1-d):
            pro += 1/2*pow(p, i)*pow(1-p, 2*k-i)*scipy.special.binom(a, d)*scipy.special.binom(2*k-a, i-d)
        
    return pro

def low_burst(a,k,p,m):
    pro = 0
    d = int(a/2) + 1
    
    for i in range(d, 2*k+1):
        for i1 in range(d, min(a,i)+1):
            pro += pow(p, i)*pow(1-p, 2*k-i+4*m)*scipy.special.binom(a, i1)*scipy.special.binom(2*k-a, i-i1)
    if(a%2==0):
        d = int(a/2)
        for i in range(d, 2*k+1-d):
            pro += 1/2*pow(p, i)*pow(1-p, 2*k-i+4*m)*scipy.special.binom(a, d)*scipy.special.binom(2*k-a, i-d)
        
    return pro