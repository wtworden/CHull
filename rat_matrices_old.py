try:
    from math import gcd
    from collections.abc import Sequence
except ImportError:
    from fractions import gcd
    from collections import Sequence

from functools import reduce



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
        if isinstance(other,Rational) or isinstance(other,int):
            return Matrix([[self._rows[i][j]*other for j in range(self._num_cols)] for i in range(self._num_rows)])
        else:
            return mult(self,other)


class Vector(Sequence):
    def __init__(self,coords):
        self._coords = tuple([Rational((c.numerator,c.denominator)) for c in coords])
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

def QQ(rat):
    return Rational(rat)

class Rational():
    def __init__(self,rat): # rat is an integer or a tuple/list of integers.
        try:
            assert rat[1] != 0
            p = rat[0]; q = rat[1]
        except TypeError:
            p = rat; q=1
        #self.numerator = p*(-2*int(q<0)+1)
        self.numerator = (p//gcd(p,q))*(-2*int(q<0)+1)
        #self.denominator = q*(-2*int(q<0)+1)
        self.denominator = (q//gcd(p,q))*(-2*int(q<0)+1)
        self.numer = self.numerator
        self.denom = self.denominator
        self._tuple = (self.numer,self.denom)

    def __add__(self,other):
        return Rational((self.numerator*other.denominator+other.numerator*self.denominator,self.denominator*other.denominator))
    
    def __radd__(self,other):
        return self.__add__(other)

    def __sub__(self,other):
        return Rational((self.numerator*other.denominator-other.numerator*self.denominator,self.denominator*other.denominator))

    def __rsub__(self,other):
        return self.__neg__().__add__(other)

    def __neg__(self):
        return Rational((-self.numerator,self.denominator))

    def __abs__(self):
        return Rational((abs(self.numerator),self.denominator))

    def __lt__(self,other):
        return self.__sub__(other).numer < 0

    def __le__(self,other):
        return self.__sub__(other).numer <= 0
    
    def __ge__(self,other):
        return self.__sub__(other).numer >= 0
    
    def __gt__(self,other):
        return self.__sub__(other).numer > 0

    def __mul__(self,other):
        try:
            return Rational((self.numerator*other.numerator,self.denominator*other.denominator))
        except AttributeError:
            return other*self

    def __rmul__(self,other):
        return self.__mul__(other)

    def __pow__(self,integer):
        if integer < 0:
            return self.inverse().__pow__(-integer)
        else:
            r = Rational((self.numerator,self.denominator))
            for i in range(integer-1):
                r = r.__mul__(self)
            return r

    def __truediv__(self,other):
        return self.__mul__(Rational((other.denominator,other.numerator)))

    def __rtruediv__(self,other):
        return self.inverse().__mul__(other)

    def __repr__(self):
        if self.denominator != 1:
            return '{}/{}'.format(self.numerator,self.denominator)
        else:
            return '{}'.format(self.numerator)

    def __hash__(self):
        return hash(self._tuple)

    def __eq__(self,other):
        return ( self.__class__ == other.__class__ ) and ( self._tuple == other._tuple )

    def __ne__(self,other):
        return ( self.__class__ != other.__class__ ) or ( self._tuple != other._tuple )

    def inverse(self):
        return Rational((self.denominator,self.numerator))





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

def projection(subspace_basis):
    A = Matrix([v for v in subspace_basis]).transpose()
    proj = A*((A.transpose()*A)**(-1))*A.transpose()
    return proj

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


def lcm(a, b):
    return a * b // gcd(a, b)

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
        inverse = [[Rational((l,d))*A_int_adjugate[i][j] for j in range(dim)] for i in range(dim)]
        return inverse
    else:
        return Matrix([[A[0][0]**(-1)]])

def dot(v,w):
    dot = 0
    for i in range(len(v)):
        dot += v.transpose()[i]*w.transpose()[i]
    return dot























