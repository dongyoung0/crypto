#!/usr/bin/env python
# coding: utf-8

# In[329]:


import numpy as np
import BasicNumTest as bnt


# In[330]:


#이산로그
def Primitive_root_gen(p):
    d = 0
    while d != 1:
        a = np.random.randint(2,p-1)
        if Exp(a,int((p-1)/2), p) % p == p-1:
            d = 1
    return a


# In[331]:


def Pre_factor(m,mx=100):
    target = m
    factor = []
    i = 0
    mx = 100
    while i < mx:
        t_factor = [2,3,5,7,11,13,17,19]
        for ft in t_factor:
            while target%ft==0 and target!=ft:
                target = int(target/ft)
                factor.append(ft)
                
        a = bnt.Bigrand(2,target)
        ft = bnt.Pollard(target,a,B=1000)
        
        if ft>1:
            factor.append(ft)
            target = int(target/ft)
        if target < 4:
            i = mx
        i += 1
    factor.append(target)
    return factor

def Prime_factor(p):
    factor = Pre_factor(p)
    prime_factor = []
    while factor != []:
        for v in factor:
            if bnt.MR_PT(int(v))==True:
                prime_factor.append(v)
                factor.remove(v)
    return set(prime_factor)

def Primitive_root_gen(p):
    prime_factor = Prime_factor(p-1)
    while True:
        a = bnt.Bigrand(1,p)
        itr = 0
        for v in prime_factor:
            itr += 1
            if bnt.Exp(a, int((p-1)/v), p) == 1:
                break
        if itr == len(prime_factor):
            return a
        
def Factoring(n):
    factors = Pre_factor(n)
    factors_set = set(factors)
    factors_dict = {}
    for factor in factors_set:
        factors_dict[factor] = factors.count(factor)
    return factors_dict


# In[333]:


Factoring(120)


# In[358]:


#BS_GS
def BSGS(a,b,p):
    N = int(bnt.np.sqrt(p-1)+1)
    GS = bnt.np.zeros(N)
    for k in range(N):
        GS[k] = b*bnt.Exp(a,2*(p-1)-k*N,p)%p
    for j in range(N):
        c = bnt.Exp(a,p-1+j,p)
        for k in range(N):
            if GS[k] == c:
                return (k*N+j)%(p-1)
            
def BSGS_adv(a,b,p):
    N = int(bnt.np.sqrt(p-1))+1
    GS = bnt.np.zeros(N)
    for k in range(N):
        GS[k] = b*bnt.Exp(a,2*(p-1)-k*N,p)%p
    GS.sort()
    for j in range(N):
        c = bnt.Exp(a,p-1+j,p)
        x = np.where(GS==c)
        if np.shape(x)==(1,1):
            if not (j==0) and x[0][0] == 0:
                return((x[0][0]*N+j))%(p-1)
    return 0
        

def CRT_general(a,n):
    temp = 0
    prod = 1
    for i in range(len(n)):
        prod *= n[i]
    for a_i, n_i in zip(a, n):
        p = prod // n_i
        temp += a_i * bnt.Inv(p, n_i) * p
    return temp % prod


# In[359]:


#Pohlig Hellman

def Pohlig_sub(a,b,p,q,r):
    if (p-1)%(q**r) != 0:
        return 0
    else:
        A = bnt.Exp(a,(p-1)//q,p)
        x = r * [0]
        X = 0
        for i in range(r):
            if i == 0:
                b1 = b
            else:
                b1 = b1*bnt.Exp(a,((p-1)//q)-x[i-1],p)%p
            B = bnt.Exp(b1,(p-1)//(q**(i+1)),p)
            temp = BSGS_adv(A,B,p)
            if temp == 0:
                x[i] = q
            else:
                x[i] = temp
            X += x[i] * (q**i)
        return X%(q**r)

def Pohlig_hellman(a,b,p):
    L = Factoring(p-1)
    q = list(L.keys())
    r = list(L.values())
    x = len(q) * [0]
    factor = len(q) * [0]
    for i in range(len(q)):
        x[i] = Pohlig_sub(a,b,p,q[i],r[i])
        factor[i] = q[i] ** r[i]
    sol = CRT_general(x,factor)
    return sol


# In[360]:


a,b,p = 6,7531,8101
sol = Pohlig_hellman(a,b,p)
sol

