

from functools import reduce
import sys

from cypari import pari
from math import gcd


def random_mat(size, a, b):
    return pari.matrix(size, size, [QQ(random.randint(a,b),1) for i in range(size) for j in range(size)])

def QQ(p,q):
    return (pari.PARI_ONE*p)/q

def adjugate(A):
    cofactor_mat = [[((-1)**(i+j))*determinant(minor(A,i,j)) for j in range(len(A))] for i in range(len(A))]
    return transpose(cofactor_mat)

def determinant(A):
    if len(A) == 2:
        return A[0][0]*A[1][1]-A[0][1]*A[1][0]
    elif len(A) > 2:
        det = 0
        i = 0
        for j in range(len(A)):
            det += A[i][j]*((-1)**(i+j))*determinant(minor(A,i,j))
        return det
    else:
        return A[0][0]


def minor(A,i,j):
    B=[row[:j]+row[j+1:] for row in A]
    to_trash = B.pop(i)
    return B


def lcm(a, b):
    return int(a * b) // gcd(a, b)

def LCM(numbers):
     return reduce(lcm, numbers)

def GCD(numbers):
     return reduce(gcd, numbers)

def inverse(A,dim):
    if len(A) > 1:
        l = LCM([A[i][j].denominator() for i in range(dim) for j in range(dim)])
        A_int = [[(A[i][j]*l).numerator() for j in range(dim)] for i in range(dim)]
#        A_int = [[(A[i][j]*l).numerator//(A[i][j]).denominator for j in range(dim)] for i in range(dim)]
        A_int_adjugate = adjugate(A_int)
        d = determinant(A_int)
        inverse = [[QQ(l,d)*A_int_adjugate[i][j] for j in range(dim)] for i in range(dim)]
        return inverse
    else:
        return pari.matrix([[A[0][0]**(-1)]])

def dot(v,w):
    return (v*w.mattranspose())[0]

def gram_schmidt(vectors):
    orthobasis = []
    orthobasis.append(vectors[0])
    for i in range(1,len(vectors)):
        orthobasis.append(vectors[i])
        for j in range(0,i):
            orthobasis[i] -= ((dot(vectors[i],orthobasis[j]))/(dot(orthobasis[j],orthobasis[j])))*orthobasis[j]
    return orthobasis


def projection(v,orthobasis):
    proj = pari.vector(len(v),[0 for _ in range(len(v))])
    for u in orthobasis:
        proj += ((dot(v,u))/(dot(u,u)))*u
    return proj



















