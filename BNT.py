import numpy as np

##########################################
# ######   Basic Operations   #############
# #########################################

def Division(a,b):
    mx = np.max([a,b])
    mn = np.min([a,b])
    q = mx//mn
    r = mx%mn
    return mn,q,r
    
def Euclidean(a,b):
    while b!=0:
        a, q,b = Division(a,b)
        print(a,q,b)
    return a

def ExtendedEuclidean(a,b):
    M = [[1,0],[0,1]]#np.eye(2)
    while b!=0:
        a, q,b = Division(a,b)
        M = np.matmul([[0,1],[1,-q]],M)
    return a,M[0,0],M[0,1]

def Inv(a,p):
    if a>p:# make a to be smaller than p
        _, a =divmod(a,p)
    gcd, _, coeff = ExtendedEuclidean(a,p)# coeff1 is for large input, coeff2 is for small input
    if gcd == 1:
        if coeff < 0:
            return coeff+p
        return coeff
    return print("No exist")

def CRT(a,p1,b,p2):
        _, crt = divmod(a*p2*Inv(p2,p1)+b*p1*Inv(p1,p2),p1*p2)
        return crt
def CRT_general(a,p):
    P = 1
    temp = 0    
    n = len(p)
    if len(a)!=n:
        print("errored input")
    for factor in p:
        P *=factor
    for i in range(n):
        temp +=a[i]*(P//p[i])*Inv(P//p[i],p[i])
    return temp%P
        


def Exp(a,x,p):
    if x ==1:
        return a%p
#     print('a,x,p:',a,x,p)
    l = x.bit_length()-2
    msb = 1 << l
    A = a
    while msb>0:
        A = (A*A)%p
        if x & msb == msb:
            A = (A*a)%p
        msb = msb >> 1
    return A

###########################################
# #########   Primality Test  ##############
# ##########################################

def FermatPT(p):
    a = np.random.randint(2,p-1)
    for i in range(80):
        if Exp(a,p-1,p)>1:
            return(p, "is composite")
    return(p," is prime")

def MR_PT(n,l=80):
    if n == 2 or n == 3:
        return True
    if n&1==0:
#         print(n, "is composite") 
        return False
    for j in range(l):
        a = Bigrand(1, n)    
        k,m = decompose(n-1)
        b0 = Exp(a,m,n)
#         print(a,k,m,b0)
        for i in range(k):
            b1 = (b0**2)%n
            if b0!=1 and b0!=n-1 and b1 ==1: 
#                 print(n, "is composite") 
                return False
            b0 = b1
        if b0 > 1:
#             print(n, "is composite")                
            return False
#     print(n, " is prime") 
    return True

def Bigrand(mn, B):
    A = (1<<63)
    if B>A:
#         print("*",B>>63)
        n = (np.random.randint(mn,1+(B>>63))<<63) +np.random.randint(mn,A)        
    else:
        n = np.random.randint(mn,B)
    return n

def Prime_Gen(B):
    n=Bigrand(1, B)
    while MR_PT(n) == False:
        n=Bigrand(1, B)
    return n

###########################################
# #####  Factoring Algorithms ##############
# ##########################################
def decompose(r):
    k = 0
    temp = r
    while(temp&1==0):
        k +=1
        temp = temp>>1
    return k,temp
        

def ExpFact(a,r,n):
    k,m = decompose(r)
    b0 = Exp(a,m,n)
    if b0==1 or b0==n-1:
#         print("False")
        return 1#False
    for i in range(k):
        b1 = (b0**2)%n
        if b1==n-1:
#             print("False")        
            return 1#False
        if b1==1:
            gcd,_,_= ExtendedEuclidean(b0-1,n)
#             print("gcd:",gcd)            
            return gcd
        b0 = b1
    return 1#False#print("end")

def Pollard(n, a=2, B=500):
    b = a
    Bf=1
    for i in range(1,B+1):
        Bf *=i
        b = int((b**i)% n)
    if b == 1:
#         print("Case1: Exponentiation_Factorization_Method")
        return ExpFact(a,Bf,n)
#     print("Case2")    
    d, _, _ = ExtendedEuclidean(b-1,n)
    return d


def trivial_factoring(n):
    factor = [2,3,5,7,11,13,17,19]
    temp = n
    counter = np.matrix([0,0,0,0,0,0,0,0])
    for i in range(8):
        while temp%factor[i]==0:
            temp = temp/factor[i]
            counter += np.eye(8, dtype=int)[i]
    if temp == 1:
        return counter,1
    return counter,0

def QS_finding_squares(I,J,n):
    for i in range(1,I):
        for j in range(J):
            candidate_int = int((np.sqrt(i*n)+j))
            E,true = trivial_factoring((candidate_int**2)%n)
            if true == 1:
                print(candidate_int,E%2)


###########################################
# #####  Primitive Root Finding ##############
# ##########################################

def Pre_factor(m,mx=100):
    target = m
    factor =[]
    i=0
    while i <mx:
        t_factor = [2,3,5,7,11,13,17,19]
        for ft in t_factor:
            while target%ft==0 and target!=ft:
                target = target//ft
                # Add factor only if it's not in the list
                if ft not in factor:
                    factor.append(ft)
        a = Bigrand(1,target)
        ft = Pollard(target,a,B=1000)
        if ft>1:
            # Add factor only if it's not in the list
            if ft not in factor: 
                factor.append(ft)
            target = target//ft
        if target <4: # target == 2 or 3
            i = mx
        i +=1
    # Add factor only if it's not in the list    
    if target not in factor:
        factor.append(target)
    return factor

def factoring(m,mx=100):
    factor = Pre_factor(m,mx)
    prime_factor=[]
    while factor != []:
        for v in factor:
            if MR_PT(int(v))==True:
                # Add factor only if it's not in the list                
                if v not in prime_factor:
                    prime_factor.append(v)
                factor.remove(v)
        for v in factor:
            factor.remove(v)
            factor = factor + Pre_factor(int(v))
    return prime_factor    

def Primitive_Gen(p):
    p_factors = factoring(p-1)
    TF = 0
    while TF == 0:
        TF = 1
        a = Bigrand(1,p)
        for r in p_factors:
            if Exp(a,(p-1)//r,p) == 1:
                TF = 0
    return a


###########################################
# #####  Finding Discrete Logarithms ##############
# ##########################################

    
def BSGS_naive(a,b,p):
    N = int(np.sqrt(p-1))+1
    GS = np.zeros(N)
    for k in range(N):
        GS[k] = b*Exp(a,2*(p-1)-k*N,p)%p
    for j in range(N):
        c = Exp(a,p-1+j,p)
        for k in range(N):
            if GS[k] == c:
                if not (j==0 and k==0):
                    return (k*N+j)%(p-1)
    return 0

def BSGS_adv(a,b,p): # input: In general, a is a primitive root
    N = int(np.sqrt(p-1))+1
    GS = np.zeros(N)
    for k in range(N):
        GS[k] = b*Exp(a,2*(p-1)-k*N,p)%p
#     print(GS)    s
    for j in range(N):
        c = Exp(a,p-1+j,p)
#         print('c',c)
        x = np.where(GS==c)
#         print('x',x)
        if np.shape(x)!=(1,0):
            if not (j==0 and x[0][0]==0):
#                 print('find')
                return(x[0][0]*N+j)%(p-1)  
    return 0  # this case should be understood the order of 'a'
