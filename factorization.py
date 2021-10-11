#!/usr/bin/env python
# coding: utf-8

# In[5]:


#import
import numpy as np
import BasicNumTest as bnt


# In[6]:


#Pollard's p-1 factoring algorithm

def Exponent_Factorization(a,r,n):
    k = 0
    temp = r
    while temp%2 == 0:
        temp = temp >> 1
        k += 1
    b = bnt.Exp(a,temp,n)
    if b == 1 or b == n-1:
        return 'failure'
    for i in range(k):
        b1 = bnt.Exp(a,temp,n)
        if b1 == n-1:
            return 'failure'
        elif b1 == 1:
            d, _, _ = bnt.ExtendedEuclidean(b-1,n)
            return d
        b = b1
        
def Pollard(n,a=2,B=500):
    b = a
    for i in range(1,B+1): 
        b = bnt.Exp(b,i,n)
    if b == 1:
        r = 1
        for j in range(1,B+1):
            r = r*j
        return print(Exponent_Factorization(b,r,n))
    else:
        d, _, _ = bnt.ExtendedEuclidean(b-1,n)
        return d

Pollard(10)


# In[9]:


#Quadratic Sieve

def trivial_factoring(n, factor):
    temp = n
    counter = len(factor) * [0]
    for i in range(len(factor)):
        while temp%factor[i]==0:
            temp = temp/factor[i]
            counter[i] += 1
            counter[i] = counter[i]
            n /= factor[i]
    if n == 1:
        return counter
    else:
        return (n, "is not a product of small prime")

def QS(I,J,n,factor):
    for i in range(1,I):
        for j in range(J):
            x = int(np.sqrt(i*n)+j)
            E = trivial_factoring(np.square(x)%n,factor)
            if type(E) == type([]):
                print(x, E)                     


# In[10]:


n = 3837523
factor = [2,3,5,7,11,13,17,19]

QS(100,30,n,factor)


# In[ ]:




