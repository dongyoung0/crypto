# -*- coding: utf-8 -*-
# %%


import BasicNumTest as bnt
import numpy as np
import time

def Sub_DLP(a,b,p,q,r):
    if (p-1)%(q**r) !=0:
        print("errored input")
        return 0
    else:
        A = bnt.Exp(a,(p-1)//q,p)
        x = np.zeros(r)
        X = 0
        for i in range(r):
            if i == 0:
                b1 = b
            else:
                b1 = b1*bnt.Exp(a,(p-1)-int(x[i-1]),p)%p
            B = bnt.Exp(b1,(p-1)//(q**(i+1)),p)
            temp = bnt.BSGS_adv(A,B,p) # A is not a primitive root
            if temp == 0:
                x[i] = q
            else:
                x[i] = temp
            X +=x[i]*(q**i)
        return X%(q**r)
def Factoring(p):
    p_factors = bnt.factoring(p)
    r = np.ones(len(p_factors))
    i = 0
    for q in p_factors:
        while p%(q**r[i]) == 0:
            r[i] +=1
        r[i] -=1
        i +=1 
    return p_factors, r

def Pohlig_Hellman(a,b,p):
    p_factors, r = Factoring(p-1)
    x = np.zeros(len(p_factors))
    factors = np.zeros(len(p_factors))
    for i in range(len(p_factors)):
        x[i] = int(Sub_DLP(a,b,p,int(p_factors[i]),int(r[i])))
        
    for j in range(len(p_factors)):
        factors[j] = p_factors[j]**r[j]
    sol = bnt.CRT_general(x,factors)
    return sol


# %%


p = bnt.Prime_Gen(10000) # 2^30크기 보통 사용하는 소수는 2^128정도
a = bnt.Primitive_Gen(p)
x = bnt.Bigrand(1,p-1)
b = bnt.Exp(a,x,p)
print('p:',p)
print('a:',a)
print('b:',b)
print('x:',x)
start = time.time()
print('PH T/F:',x==Pohlig_Hellman(a,b,p),'time:',time.time()-start)
start = time.time()    
print('advanced T/F:',x==bnt.BSGS_adv(a,b,p),'time:',time.time()-start)
# print('naive T/F:',x==BSGS_naive(a,b,p))


# %%


q1 = bnt.Prime_Gen(1000)
q2 = bnt.Prime_Gen(1000)
p = bnt.MR_PT(1+q1**4*q2**4)


# %%


Sub_DLP(413,562,593,2,4)


# %%


592%2**4


# %%


bnt.factoring(13**10*11**15)


# %%


bnt.Pre_factor(13**12*11**15)


# %%


bnt.BSGS_naive(10,1,11)


# %%


bnt.Exp(3,1,11)


# %%


a = np.array([1,2,3])
p = np.array([3,5,7])
bnt.CRT_general(a,p)


# %%


a = 413
b = 562
p =593
q =2
r =4    

A = bnt.Exp(a,(p-1)//q,p)
# print('A',A)
x = np.zeros(r)
# print('x',x)
X = 0
for i in range(r):
    if i == 0:
        b1 = b
    else:
        b1 = b1*bnt.Exp(a,(p-1)-int(x[i-1]),p)%p
    B = bnt.Exp(b1,(p-1)//(q**(i+1)),p)
    print('i',i,B)
#     temp = bnt.BSGS_adv(A,B,p)
    x[i] = temp
    X +=x[i]*(q**i)
print('x',x,'X',X)
# return X%(q**r)


# %%


def BSGS_adv(a,b,p):
    N = int(np.sqrt(p-1))+1
    GS = np.zeros(N)
    for k in range(N):
        GS[k] = b*bnt.Exp(a,2*(p-1)-k*N,p)%p
#     print(GS)    s
    for j in range(N):
        c = bnt.Exp(a,p-1+j,p)
#         print('c',c)
        x = np.where(GS==c)
#         print('x',x)
        if np.shape(x)!=(1,0):
            if not (j==0 and x[0][0]==0):
#                 print('find')
                return(x[0][0]*N+j)%(p-1)  
    return 0

print(bnt.BSGS_adv(592,392,593))


# %%


a = 592 = -1
b = 392 = even
p = 593
N = int(np.sqrt(p-1))+1
GS = np.zeros(N)
for k in range(N):
    GS[k] = b*bnt.Exp(a,2*(p-1)-k*N,p)%p
#     print(GS)    s
for j in range(N):
    c = bnt.Exp(a,p-1+j,p)
#         print('c',c)
    x = np.where(GS==c)
#         print('x',x)
    if np.shape(x)!=(1,0):
        if not (j==0 and x[0][0]==0):
            print('find')
            print((x[0][0]*N+j)%(p-1)  )
print('end')
# return 0


# %%




