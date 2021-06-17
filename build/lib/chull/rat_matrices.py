try:
    from collections.abc import Sequence
except ImportError:
    from collections import Sequence

from functools import reduce
import sys

if sys.platform == 'linux':
    from sympy import QQ as mpq
    from math import gcd
    def lcm(a, b):
        return a * b // gcd(a, b)

elif sys.platform == 'darwin':
    from gmpy2 import mpq, gcd, lcm



class Matrix(Sequence):
    def __init__(self,rows):
        self._rows = [Vector(rows[i]) for i in range(len(rows))]
        self._num_rows = len(self._rows)
        self._num_cols = len(self._rows[0])
        self._columns = [Vector([self._rows[i][j] for i in range(self._num_rows)]) for j in range(self._num_cols)]
        #super().__init__()

    def __getitem__(self, i):
        return self._rows[i]
    def __len__(self):
        return len(self._rows)

    def __add__(self,other):
        return Matrix([[self._rows[i][j]+other._rows[i][j] for j in range(self._num_cols)] for i in range(self._num_rows)])
    
    def __sub__(self,other):
        return Matrix([[self._rows[i][j]-other._rows[i][j] for j in range(self._num_cols)] for i in range(self._num_rows)])

    def __neg__(self):
        return Matrix([[-self._rows[i][j] for j in range(self._num_cols)] for i in range(self._num_rows)])

    def __abs__(self):
        return Matrix([[abs(self._rows[i][j]) for j in range(self._num_cols)] for i in range(self._num_rows)])

    def __truediv__(self,other):
        return Matrix([[self._rows[i][j]/other for j in range(self._num_cols)] for i in range(self._num_rows)])

    def __mul__(self,other):
        return self.mult(other)

    def __rmul__(self,other):
        return self.mult(other)

    def __pow__(self,integer):
        if integer < 0:
            return self.inverse().__pow__(-integer)
        else:
            A = Matrix(self._rows)
            for i in range(integer-1):
                A = A.__mul__(self)
            return A

    def __repr__(self):
        rep_str = ''
        as_strings = [[str(self.row(i)[j]) for j in range(self._num_cols)] for i in range(self._num_rows)]
        col_widths = [max([len(as_strings[i][j]) for i in range(self._num_rows)]) for j in range(self._num_cols)]
        for i in range(self._num_rows-1):
            rep_str += '[ ' + ''.join([' '*(col_widths[j]-len(as_strings[i][j])) + str(self.row(i)[j])+' ' for j in range(len(self.row(i)))])+']\n'
        rep_str += '[ ' + ''.join([' '*(col_widths[j]-len(as_strings[self._num_rows-1][j])) + str(self.rows()[-1][j])+' ' for j in range(len(self.rows()[-1]))])+']'
        return rep_str


    def __hash__(self):
        return hash(self.rows())

    def __eq__(self,other):
        return ( self.__class__ == other.__class__ ) and ( self.rows() == other.rows() )

    def __ne__(self,other):
        return ( self.__class__ != other.__class__ ) or ( self.rows() != other.rows() )

    def row(self,i):
        return self.rows()[i]

    def col(self,i):
        return self.cols()[i]

    def rows(self):
        return self._rows

    def cols(self):
        return self._columns


    def transpose(self):
        return Matrix(transpose(self))

    def minor(self,i,j):
        return Matrix(minor(self,i,j))

    def inverse(self):
        return Matrix(inverse(self,self._num_rows))

    def determinant(self):
        return determinant(self)

    def mult(self,other):
        if isinstance(other,Matrix) or isinstance(other,Vector):
            return mult(self,other)
        else:
            return Matrix([[self._rows[i][j]*other for j in range(self._num_cols)] for i in range(self._num_rows)])


class Vector(Sequence):
    def __init__(self,coords):
#        self._coords = tuple([mpq(c.numerator,c.denominator) for c in coords])
        self._coords=tuple(coords)
        self.dim = len(self._coords)
        #super().__init__()

    def __getitem__(self, i):
        return self._coords[i]
    def __len__(self):
        return len(self._coords)

    def __add__(self,other):
        return Vector([self._coords[i]+other._coords[i] for i in range(self.dim)])
    
    def __sub__(self,other):
        return Vector([self._coords[i]-other._coords[i] for i in range(self.dim)])

    def __neg__(self):
        return Vector([-self._coords[i] for i in range(self.dim)])

    def __abs__(self):
        return Vector([abs(self._coords[i]) for i in range(self.dim)])

    def __truediv__(other):
        return Vector([self._coords[i]/other for i in range(self.dim)])

    def __mul__(self,other):
        return self.mult(other)

    def __rmul__(self,other):
        return self.mult(other)

    def __repr__(self):
        return str(self._coords)

    def __hash__(self):
        return hash(self._coords)

    def __eq__(self,other):
        return ( self.__class__ == other.__class__ ) and ( self._coords == other._coords )

    def __ne__(self,other):
        return ( self.__class__ != other.__class__ ) or ( self._coords != other._coords )

    def __le__(self,other):
        return self._coords <= other._coords
    def __lt__(self,other):
        return self._coords < other._coords
    def __ge__(self,other):
        return self._coords >= other._coords
    def __gt__(self,other):
        return self._coords > other._coords

    def as_tuple(self):
        return self._coords

    def mult(self,other):
        if isinstance(other,Vector):
            return self.dot_product(other)
        elif isinstance(other,Matrix):
            return mult(self,other)
        else:
            return Vector([other*self._coords[i] for i in range(self.dim)])

    def dot_product(self, other_vec):
        dot = 0
        for i in range(len(self._coords)):
            dot += self._coords[i]*other_vec._coords[i]
        return dot

    def ns(self):
        return self.dot_product(self)



def mult(A,B): #returns A*B, where A and B are vectors/matrices
    is_vec = [0,0]
    if isinstance(B,Vector):
        B = [[B[i]] for i in range(len(B))]
        is_vec[1] = 1
    if isinstance(A,Vector):
        A = [[A[i] for i in range(len(A))]]
        is_vec[0] = 1
    size = (len(A),len(B[0]))
    C = [[0 for j in range(size[1])] for i in range(size[0])]
    for i in range(len(C)):
        for j in range(len(C[0])):
            for k in range(len(B)):
                C[i][j] += A[i][k]*B[k][j]
    if is_vec == [0,0]:
        return Matrix(C)
    elif is_vec == [0,1]:
        return Vector([C[i][0] for i in range(len(C))])
    elif is_vec == [1,0]:
        return Vector(C[0])
    else:
        return C[0][0]


def transpose(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]

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


#def lcm(a, b):
#    return a * b // gcd(a, b)

def LCM(numbers):
     return reduce(lcm, numbers)

def GCD(numbers):
     return reduce(gcd, numbers)

def inverse(A,dim):
    if len(A) > 1:
        l = LCM([A[i][j].denominator for i in range(dim) for j in range(dim)])
        A_int = [[(A[i][j]*l).numerator for j in range(dim)] for i in range(dim)]
#        A_int = [[(A[i][j]*l).numerator//(A[i][j]).denominator for j in range(dim)] for i in range(dim)]
        A_int_adjugate = adjugate(A_int)
        d = determinant(A_int)
        inverse = [[mpq(l,d)*A_int_adjugate[i][j] for j in range(dim)] for i in range(dim)]
        return inverse
    else:
        return Matrix([[A[0][0]**(-1)]])

def dot(v,w):
    dot = 0
    for i in range(len(v)):
        dot += v.transpose()[i]*w.transpose()[i]
    return dot

def gram_schmidt(vectors):
    orthobasis = []
    orthobasis.append(vectors[0])
    for i in range(1,len(vectors)):
        orthobasis.append(vectors[i])
        for j in range(0,i):
            orthobasis[i] -= ((vectors[i]*orthobasis[j])/(orthobasis[j]*orthobasis[j]))*orthobasis[j]
    return orthobasis


def projection(v,orthobasis):
    proj = Vector([0 for _ in range(len(v))])
    for u in orthobasis:
        proj += ((v*u)/(u*u))*u
    return proj



















