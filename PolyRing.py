#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import

import numpy as np
import BasicNumTest as bnt


# In[21]:


class RingElement:
    def __init__(self,num,mod): #Prime Field
        self.num = num % mod
        self.mod = mod
#        if bnt.MR_PT(prime) == False:
#            raise ValueError('The second input {} is not a prime'.format(prime))
        if self.num >= self.mod or self.num < 0:
            raise ValueError('The first input {} is not in the range 0 to {}'.format(num,prime-1))
        
    def __eq__(self, other):
        if other is None:
            return False
        return self.num == other.num and self.mod == other.mod
        
        
    def __add__(self, other):
        if self.mod != other.mod:
            raise TypeError('We cannot add two numbers in different rings.')
        num = (self.num + other.num) % self.mod
        return self.__class__(num, self.mod)
    
    def __sub__(self, other):
        if self.mod != other.mod:
            raise TypeError('We cannot subtract two numbers in different rings.')
        num = (self.num - other.num) % self.mod
        return self.__class__(num, self.mod)
    
    def __mul__(self, other):
        if self.mod != other.mod:
            raise TypeError('We cannot multiply two numbers in different rings.')
        num = (self.num * other.num) % self.mod
        return self.__class__(num, self.mod)
        
    def __pow__(self, exponent):
        a = exponent % (self.mod - 1)
        num = bnt.Exp(self.num, a, self.mod)
        return self.__class__(num, self.mod)
    
    def __neg__(self):
        num = self.mod - self.mod
        return self.__class__(num, self.mod)
    
    def __rmul__(self, coeff):
        num = (coeff * self.num) % self.mod
        return self.__class__(num, self.mod)
    
    def __truediv__(self, other):
        num = self.num * bnt.Inv(other.num,self.mod)
        return self.__class__(num, self.mod)
    
    def inverse(self):
        num = bnt.Inv(self.num,self.mod)
        if type(num) == str:
            return None
        return self.__class__(num,self.mod)


# In[6]:


class PolyRing: ## x^N -1  / Quotient space. N은 항상 동일
    def __init__(self,N,coeffs,mod):
        self.N = N
        self.mod = mod
        self.coeffs = coeffs # Are coefficients ring element? integer?
        for i in range(self.N):
            if self.coeffs[i].mod != self.mod:
                raise ValueError('The second input {} is not a list of same ring elements'.format(coeffs))
        if len(coeffs) != N:
            raise ValueError('error')
#         if coeffs[N-1].num == 0:
#             raise ValueError('error')
                
    def __eq__(self,other):
        return self.N == other.N and self.coeffs == other.coeffs and self.mod == other.mod
    
    def __add__(self,other):
        if self.mod != other.mod:
            raise TypeError('We cannot add two polynomials in different rings.')
        c = list()
        
        for i in range(other.N):
            c.append(self.coeffs[i]+other.coeffs[i])

        return self.__class__(self.N,c,self.mod)
        
    def __sub__(self,other):
        if self.mod != other.mod:
            raise TypeError('We cannot add two polynomials in different rings.')
        c = list()
        for i in range(other.N):
            c.append(self.coeffs[i]-other.coeffs[i])
        return self.__class__(self.N,c,self.mod)
        
    def __mul__(self,other):
        if self.mod != other.mod:
            raise TypeError('We cannot add two polynomials in different rings.')
        c = list()
        for i in range(self.N): # output polynomial's coefficient
            temp = RingElement(0,self.mod)
            for j in range (other.N):
                temp = temp + (self.coeffs[j]*other.coeffs[(i-j)%self.N])
            c.append(temp)
        return self.__class__(self.N,c,self.mod)
    
    def __rmul__(self,coeff):
        c = []
        temp = RingElement(0,self.mod)
        for i in range(self.N):
            temp = coeff * self.coeffs[i]
            c.append(temp)
        return self.__class__(self.N,c,self.mod)
            

    def deg(self):
        degree = self.N
        temp = 0
        while temp == 0:
            degree -= 1
            if degree <0:
                return 0
            temp = self.coeffs[degree].num
        return degree
    
    def poly_div(self, other):
        if self.deg() >= other.deg():
            mx_coeffs = self.coeffs
            mx_deg = self.deg()
            mn_coeffs = other.coeffs
            mn_deg = other.deg()
        else:
            mx_coeffs = other.coeffs
            mx_deg = other.deg()
            mn_coeffs = self.coeffs
            mn_deg = self.deg()
        
        q = self.N*[RingElement(0,self.mod)]
        r = mx_coeffs
        
        for i in range(mx_deg - mn_deg+1):
            if mx_coeffs[mx_deg-i].num!=0:
#                 d = mx_coeffs[mx_deg-i].num/mn_coeffs[mn_deg].num
                d = mx_coeffs[mx_deg-i]/mn_coeffs[mn_deg]
#                 q[mn_deg-i] = RingElement(d,self.mod)
                q[mn_deg-i] = d
    #             mx_coeffs[mx_deg-mn_deg-i:mx_deg-i] -= d * mn_coeffs
                for k in range(mn_deg+1):
                    mx_coeffs[mx_deg-mn_deg-i+k] -= d * mn_coeffs[k]

        mn = self.__class__(self.N,mn_coeffs,self.mod)
        quotient = self.__class__(self.N,q,self.mod)
        remainder = self.__class__(self.N,r,self.mod)
        
        return mn, quotient, remainder

    def Euclidean(self,other):
        a = self
        b = other
        while b.coeffs != self.N * [RingElement(0,self.mod)]:
            a, q, b = a.poly_div(b)
            for i in range(self.N):
                print(a.coeffs[i].num, q.coeffs[i].num, b.coeffs[i].num)
        return a
    
    def ExtendedEuclidean(self,other):
        a = self
        b = other
        one_coeffs = [RingElement(1,self.mod)] + (self.N-1) * [RingElement(0,self.mod)]
        one = self.__class__(self.N,one_coeffs,self.mod)
        zero_coeffs = self.N * [RingElement(0,self.mod)]
        zero = self.__class__(self.N,zero_coeffs,self.mod)
        M = [[one,zero],[zero,one]]
        while b.coeffs != self.N * [RingElement(0,self.mod)]:
            a,q,b = a.poly_div(b)
            M = [[M[1][0], M[1][1]],[M[0][0]-q*M[1][0], M[0][1]-q*M[1][1]]]
        return a, M[0][0], M[0][1]


# In[16]:


# f = [a0, a1, ..., an]
# g = [b0, b1, ..., bm]
# f*g = [a0b0, a0b1+a1b0, a0b2+a1b1+a2b0, ..., anbm]
# c_i = a0bi + a1b(i-1) + ... + 

def print_poly(f):
    for i in range(len(f)):
        print(f[i].num)

def PolyRing_to_Poly(ringpoly):
    temp = ringpoly.coeffs
    degree = ringpoly.deg()
    poly = temp[:degree+1]
    return poly
   

def poly_add(f,g):
    if len(f)>len(g):
        c = f
        for i in range(len(g)):
            c[i] +=g [i]
    else:
        c = f + (len(g)-len(f))*[RingElement(0,f[0].mod)]
        for i in range(len(g)):
            c[i] += g[i]
    return c

def poly_sub(f,g):
    if len(f)>len(g):
        c = f
        for i in range(len(g)):
            c[i] -=g [i]
    else:
        c = f + (len(g)-len(f))*[RingElement(0,f[0].mod)]
        for i in range(len(g)):
            c[i] -= g[i]
    return c


def poly_mul(f,g):
    deg_f = len(f)-1
    deg_g = len(g)-1
    c = (deg_f+deg_g+1)*[RingElement(0,f[0].mod)]
    for i in range(deg_f+1):
        for j in range(deg_g+1):
            c[i+j] += f[i]*g[j]
    return c


def poly_div(f,g):
    if len(f) < len(g):
        return g, [RingElement(0,f[0].mod)], f
    q = (len(f)-len(g)+1)*[RingElement(0,f[0].mod)]
    r = f
    d = RingElement(0,f[0].mod)
    for i in range(len(f)-len(g)+1):
        if r[-1-i] == RingElement(0,f[0].mod):
            q[-1-i] = RingElement(0,f[0].mod)
        else:
            d =  r[-1-i] / g[-1]
            q[-1-i] = d
            for k in range(len(g)):
                r[-1-i-k] -= d * g[-1-k]
    while r[-1] == RingElement(0,f[0].mod):
        r = r[:-1]
        if r == []:
            r = [RingElement(0,f[0].mod)]
            return g, q, r
    return g, q, r


def poly_ExtendedEuclidean(f,g):
    F = f
    G = g
    one = RingElement(1,f[0].mod)
    zero = RingElement(0,f[0].mod)
    M = [[[one],[zero]],[[zero],[one]]]
    while G != len(G) * [zero]:
        F, q, G = poly_div(F,G)
#         M = [[M[1][0], M[1][1]],[M[0][0]-q*M[1][0], M[0][1]-q*M[1][1]]]
        m00 = M[1][0]
        m01 = M[1][1]
        m10 = poly_sub(M[0][0],poly_mul(q,M[1][0]))
        m11 = poly_sub(M[0][1],poly_mul(q,M[1][1]))
        M = [[m00,m01],[m10,m11]]
    return F, M[0][0], M[0][1]
    
def poly_inv(f,p):
    gcd, inv, _ = poly_ExtendedEuclidean(f,p)
    if gcd == [RingElement(1,f[0].mod)]:
        return inv
    return None
    


# In[24]:


# f = -1 + x + x^2 - x^4 + x^6 + x^9 - x^10
# f_p = 1 + 2x + 2x^3 + 2x^4 + x^5 + 2x^7 + x^8 + 2x^9 (mod3)

prime = 3
a0 = RingElement(0,prime)
a1 = RingElement(1,prime)
a_1 = RingElement(-1,prime)

f = [a_1, a1, a1, a0, a_1, a0, a1, a0, a0, a1, a_1]
p = [RingElement(-1,f[0].mod)]+10*[RingElement(0,a[0].mod)]+[RingElement(1,f[0].mod)]

#poly_ExtendedEuclidean(a,p)
print_poly(poly_inv(f,p))


# In[25]:


# f = -1 + x + x^2 - x^4 + x^6 + x^9 - x^10
# f_q = 5+9x+6x^2+16x^3+4x^4+15x^5+16x^6+22x^7+20x^8+18x^9+30x^10 (mod32)
# g = x^11 -1
# N = 11
prime = 32
a0 = RingElement(0,prime)
a1 = RingElement(1,prime)
a_1 = RingElement(-1,prime)

a = [a_1, a1, a1, a0, a_1, a0, a1, a0, a0, a1, a_1]
#a = PolyRing_to_Poly(PolyRing(11,a_l,32))
p = [RingElement(-1,a[0].mod)]+10*[RingElement(0,a[0].mod)]+[RingElement(1,a[0].mod)]

poly_ExtendedEuclidean(a,p)


# In[ ]:




