#!/usr/bin/env python
# coding: utf-8

# In[7]:


#ECC
import numpy as np
import BasicNumTest as bnt
from random import randint
import hashlib
import hmac

class FieldElement:

    def __init__(self, num, prime):
        if num >= prime or num < 0:
            error = 'Num {} not in field range 0 to {}'.format(
                num, prime - 1)
            raise ValueError(error)
        self.num = num
        self.prime = prime

    def __eq__(self, other):
        if other is None:
            return False
        return self.num == other.num and self.prime == other.prime

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        if self.prime != other.prime:
            raise TypeError('Cannot add two numbers in different Fields')
        num = (self.num + other.num) % self.prime
        return self.__class__(num, self.prime)

    def __sub__(self, other):
        if self.prime != other.prime:
            raise TypeError('Cannot subtract two numbers in different Fields')
        num = (self.num - other.num) % self.prime
        return self.__class__(num, self.prime)

    def __mul__(self, other):
        if self.prime != other.prime:
            raise TypeError('Cannot multiply two numbers in different Fields')
        num = (self.num * other.num) % self.prime
        return self.__class__(num, self.prime)

    def __pow__(self, exponent):
        n = exponent % (self.prime - 1)
        num = pow(self.num, n, self.prime)
        return self.__class__(num, self.prime)

    def __truediv__(self, other):
        if self.prime != other.prime:
            raise TypeError('Cannot divide two numbers in different Fields')
        num = (self.num * pow(other.num, self.prime - 2, self.prime)) % self.prime
        return self.__class__(num, self.prime)

    def __rmul__(self, coefficient):
        num = (self.num * coefficient) % self.prime
        return self.__class__(num=num, prime=self.prime)

class Point:

    def __init__(self, x, y, a, b):
        self.a = a
        self.b = b
        self.x = x
        self.y = y
        if self.x is None and self.y is None:
            return
        if self.y**2 != self.x**3 + a * x + b:
            raise ValueError('({}, {}) is not on the curve'.format(x, y))
    # end::source1[]

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y             and self.a == other.a and self.b == other.b

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        if self.a != other.a or self.b != other.b:
            raise TypeError('Points {}, {} are not on the same curve'.format(self, other))

        if self.x is None:
            return other

        if other.x is None:
            return self
        
        if self.x == other.x and self.y != other.y:
            return self.__class__(None, None, self.a, self.b)

        if self.x != other.x:
            s = (other.y - self.y) / (other.x - self.x)
            x = s**2 - self.x - other.x
            y = s * (self.x - x) - self.y
            return self.__class__(x, y, self.a, self.b)

        if self == other and self.y == 0 * self.x:
            return self.__class__(None, None, self.a, self.b)

        if self == other:
            s = (3 * self.x**2 + self.a) / (2 * self.y)
            x = s**2 - 2 * self.x
            y = s * (self.x - x) - self.y
            return self.__class__(x, y, self.a, self.b)

    def __rmul__(self, coefficient):
        coef = coefficient
        temp = self  
        result = self.__class__(None, None, self.a, self.b)  # <2>
        while coef:
            if coef & 1:  
                result += temp
            temp += temp
            coef >>= 1  
        return result


def hash256(s):
    '''two rounds of sha256'''
    return hashlib.sha256(hashlib.sha256(s).digest()).digest()

A = 0
B = 7
P = 2**256 - 2**32 - 977

class S256Field(FieldElement):

    def __init__(self, num, prime=None):
        super().__init__(num=num, prime=P)

class S256Point(Point):

    def __init__(self, x, y, a=None, b=None):
        a, b = S256Field(A), S256Field(B)
        if type(x) == int:
            super().__init__(x=S256Field(x), y=S256Field(y), a=a, b=b)
        else:
            super().__init__(x=x, y=y, a=a, b=b)  

    def __repr__(self):
        if self.x is None:
            return 'S256Point(infinity)'
        else:
            return 'S256Point({}, {})'.format(self.x, self.y)

    def __rmul__(self, coefficient):
        coef = coefficient % N  
        return super().__rmul__(coef)

    def verify(self, z, sig):
        s_inv = pow(sig.s, N - 2, N)
        u = z * s_inv % N  
        v = sig.r * s_inv % N  
        total = u * G + v * self  
        return total.x.num == sig.r  


# In[8]:


prime = 223
a = FieldElement(0, prime)
b = FieldElement(7, prime)


# (170,142) + (60,139) = (220, 181)
x1 = FieldElement(170, prime)
y1 = FieldElement(142, prime)
p1 = Point(x1, y1, a, b)
x2 = FieldElement(60, prime)
y2 = FieldElement(139, prime)
p2 = Point(x2, y2, a, b)
sum1 = p1 + p2
print('(170,142) + (60,139) = (%d, %d)' %(sum1.x.num, sum1.y.num))


# In[9]:


#ECDSA


# In[21]:


def ECDSA_setup():
    N = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
    G = S256Point(0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
    0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8)
    sk = randint(1,N-1)
    return N,G,sk

def ECDSA_signature(msg,sk,N):
    z = int.from_bytes(hash256(msg), 'big')
    r = 0
    s = 0
    while r == 0 or s == 0:
        k = randint(1,N-1)
        r = (k*G).x.num % N
        k_inv = pow(k, N-2, N)
        s = k_inv*(z + r*sk) % N
    return r,s

def ECDSA_verification(r,s,sk,msg,N,G):
    z = int.from_bytes(hash256(msg), 'big')
    if r >= N or s >= N:
        print('invalid1') 
    else:
        w = pow(s, N-2, N)
        u1 = z*w % N
        u2 = r*w % N
        Q = sk*G
        Pt = u1*G + u2*Q
        if Pt == Point(None, None, a, b):
            print('invalid2')
        else:
            if (Pt.x.num - r) % N == 0:
                print('valid')
            else:
                print('invalid3')


# In[23]:


msg = b'Elliptic_Curve23'
N,G,sk = ECDSA_setup()
r,s = ECDSA_signature(msg,sk,N)


# In[24]:


ECDSA_verification(r,s,sk,msg,N,G)


# In[ ]:




