
import chull.rat_matrices as rat
from random import randint
import chull.convex_hull as ch
from timeit import timeit
import sys
#from sage.all import Polyhedron
if sys.platform == 'linux':
    from sympy import QQ as mpq
elif sys.platform == 'darwin':
    from gmpy2 import mpq


def in_sage():
    try:
        import sage
        return True
    except ModuleNotFoundError:
        return False

def tuples_to_vecs(sage_tuples):
    vectors = [rat.Vector([mpq(c.numer(),c.denom()) for c in p]) for p in sage_tuples]
    return vectors


def get_test_points(dim,sage=False):
    if not sage:
        points = [rat.Vector([mpq(14,9),mpq(19,17),mpq(12,19),mpq(-2)]), rat.Vector([mpq(-14,9),mpq(-19,17),mpq(-12,19),mpq(2)]),rat.Vector([mpq(6),mpq(1,5),mpq(-3,2),mpq(-5,18)]),rat.Vector([mpq(-6),mpq(-1,5),mpq(3,2),mpq(5,18)]),rat.Vector([mpq(-17,19),mpq(-5,9),mpq(-9,17),mpq(1,19)]),rat.Vector([mpq(17,19),mpq(5,9),mpq(9,17),mpq(-1,19)]),rat.Vector([mpq(7),mpq(5,12),mpq(11,20),mpq(5)]),rat.Vector([mpq(-7),mpq(-5,12),mpq(-11,20),mpq(-5)]),rat.Vector([mpq(-19,14),mpq(2,19),mpq(-2),mpq(-1,4)]),rat.Vector([mpq(19,14),mpq(-2,19),mpq(2),mpq(1,4)]),rat.Vector([mpq(-9),mpq(18,7),mpq(10,7),mpq(3)]),rat.Vector([mpq(9),mpq(-18,7),mpq(-10,7),mpq(-3)]),rat.Vector([mpq(19,18),mpq(1,4),mpq(-19,11),mpq(-5,12)]),rat.Vector([mpq(-19,18),mpq(-1,4),mpq(19,11),mpq(5,12)]),rat.Vector([mpq(-5,11),mpq(-8,17),mpq(7,13),mpq(-18,7)]),rat.Vector([mpq(5,11),mpq(8,17),mpq(-7,13),mpq(18,7)]),rat.Vector([mpq(17,5),mpq(5,3),mpq(0),mpq(17,12)]),rat.Vector([mpq(-17,5),mpq(-5,3),mpq(0),mpq(-17,12)]), rat.Vector([mpq(-3),mpq(-20),mpq(15,13),mpq(19,9)]), rat.Vector([mpq(3),mpq(20),mpq(-15,13),mpq(-19,9)]), rat.Vector([mpq(-4,3),mpq(-4,5),mpq(-1,19),mpq(-5,4)]), rat.Vector([mpq(4,3),mpq(4,5),mpq(1,19),mpq(5,4)]), rat.Vector([mpq(-7),mpq(-17,18),mpq(-5,14),mpq(-17,19)]), rat.Vector([mpq(7),mpq(17,18),mpq(5,14),mpq(17,19)]), rat.Vector([mpq(-13,5),mpq(5,7),mpq(-3,2),mpq(-14,15)]), rat.Vector([mpq(13,5),mpq(-5,7),mpq(3,2),mpq(14,15)]), rat.Vector([mpq(-8,15),mpq(-17,9),mpq(-6,19),mpq(-5)]), rat.Vector([mpq(8,15),mpq(17,9),mpq(6,19),mpq(5)]), rat.Vector([mpq(6),mpq(20),mpq(-20,11),mpq(12,19)]), rat.Vector([mpq(-6),mpq(-20),mpq(20,11),mpq(-12,19)])]
        if dim == 4:
            return points
        elif dim == 3:
            return [rat.Vector(v[:-1]) for v in points]
    else:
        points = [vector([QQ((14,9)),QQ((19,17)),QQ((12,19)),QQ((-2))]), vector([QQ((-14,9)),QQ((-19,17)),QQ((-12,19)),QQ((2))]),vector([QQ((6)),QQ((1,5)),QQ((-3,2)),QQ((-5,18))]),vector([QQ((-6)),QQ((-1,5)),QQ((3,2)),QQ((5,18))]),vector([QQ((-17,19)),QQ((-5,9)),QQ((-9,17)),QQ((1,19))]),vector([QQ((17,19)),QQ((5,9)),QQ((9,17)),QQ((-1,19))]),vector([QQ((7)),QQ((5,12)),QQ((11,20)),QQ((5))]),vector([QQ((-7)),QQ((-5,12)),QQ((-11,20)),QQ((-5))]),vector([QQ((-19,14)),QQ((2,19)),QQ((-2)),QQ((-1,4))]),vector([QQ((19,14)),QQ((-2,19)),QQ((2)),QQ((1,4))]),vector([QQ((-9)),QQ((18,7)),QQ((10,7)),QQ((3))]),vector([QQ((9)),QQ((-18,7)),QQ((-10,7)),QQ((-3))]),vector([QQ((19,18)),QQ((1,4)),QQ((-19,11)),QQ((-5,12))]),vector([QQ((-19,18)),QQ((-1,4)),QQ((19,11)),QQ((5,12))]),vector([QQ((-5,11)),QQ((-8,17)),QQ((7,13)),QQ((-18,7))]),vector([QQ((5,11)),QQ((8,17)),QQ((-7,13)),QQ((18,7))]),vector([QQ((17,5)),QQ((5,3)),QQ((0)),QQ((17,12))]),vector([QQ((-17,5)),QQ((-5,3)),QQ((0)),QQ((-17,12))]), vector([QQ((-3)),QQ((-20)),QQ((15,13)),QQ((19,9))]), vector([QQ((3)),QQ((20)),QQ((-15,13)),QQ((-19,9))]), vector([QQ((-4,3)),QQ((-4,5)),QQ((-1,19)),QQ((-5,4))]), vector([QQ((4,3)),QQ((4,5)),QQ((1,19)),QQ((5,4))]), vector([QQ((-7)),QQ((-17,18)),QQ((-5,14)),QQ((-17,19))]), vector([QQ((7)),QQ((17,18)),QQ((5,14)),QQ((17,19))]), vector([QQ((-13,5)),QQ((5,7)),QQ((-3,2)),QQ((-14,15))]), vector([QQ((13,5)),QQ((-5,7)),QQ((3,2)),QQ((14,15))]), vector([QQ((-8,15)),QQ((-17,9)),QQ((-6,19)),QQ((-5))]), vector([QQ((8,15)),QQ((17,9)),QQ((6,19)),QQ((5))]), vector([QQ((6)),QQ((20)),QQ((-20,11)),QQ((12,19))]), vector([QQ((-6)),QQ((-20)),QQ((20,11)),QQ((-12,19))])]
        if dim == 4:
            return points
        elif dim == 3:
            return [vector(v[:-1]) for v in points]

def get_rand_points(num_points, width, dim):
    sage_points = []
    rat_points = []
    for i in range(num_points):
        sage_point = [0 for i in range(dim)]
        rat_point = [0 for i in range(dim)]
        for j in range(dim):
            p = randint(-width,width)
            q = randint(1,width)
            rat_point[j] = mpq(p,q)
            if in_sage():
                sage_point[j] = QQ((p,q))
        if in_sage():
            sage_points.append(vector(sage_point))
        rat_points.append(rat.Vector(rat_point))
    all_sage_points=[]
    all_rat_points=[]
    for p in sage_points:
        all_sage_points.append(p)
        all_sage_points.append(-p)
    for p in rat_points:
        all_rat_points.append(p)
        all_rat_points.append(-p)
    return all_sage_points, all_rat_points

def check_hulls(num_points, width, dim, rand=True):
    if not rand:
        sage_points, rat_points = get_test_points(4,True),get_test_points(4,False)
    else:
        sage_points, rat_points = get_rand_points(num_points, width, dim)
    sage_hull = Polytope(sage_points)
    rat_hull = ch.quick_hull(rat_points, dim)
    
    sage_verts = sorted([tuple(v.vector()) for v in sage_hull.vertices()])
    rat_verts = [v.coords().as_tuple() for v in rat_hull.facets(0)]
    rat_verts = sorted([tuple([QQ((v[i].numerator,v[i].denominator)) for i in range(dim)]) for v in rat_verts])
    if sage_verts != rat_verts:
        return False, sage_hull,rat_hull

    for d in range(1,rat_hull.dim):
        sage_facets = sorted([sorted([tuple(v.vector()) for v in f.vertices()]) for f in sage_hull.faces(d)])
        rat_facets = [[v.coords().as_tuple() for v in f.vertices()] for f in rat_hull.facets(d)]
        rat_facets = sorted([sorted([tuple([QQ((v[i].numerator,v[i].denominator)) for i in range(dim)]) for v in rat_facets[i]]) for i in range(len(rat_facets))])
        if sage_facets != rat_facets:
            return False, sage_hull,rat_hull
        else:
            continue
    return True, sage_hull,rat_hull

def get_bad_points():
    points=[rat.Vector([mpq(-10,11), mpq(22,13), mpq(-27,10), mpq(-3,14)]),rat.Vector([mpq(10,11), mpq(-22,13), mpq(27,10), mpq(3,14)]),rat.Vector([mpq(5,11), mpq(0,1), mpq(-1,18), mpq(24,17)]),rat.Vector([mpq(-5,11), mpq(0,1), mpq(1,18), mpq(-24,17)]),rat.Vector([mpq(2,15), mpq(-7,26), mpq(11,5), mpq(3,1)]),rat.Vector([mpq(-2,15), mpq(7,26), mpq(-11,5), mpq(-3,1)]),rat.Vector([mpq(12,25), mpq(17,6), mpq(5,3), mpq(-1,3)]),rat.Vector([mpq(-12,25), mpq(-17,6), mpq(-5,3), mpq(1,3)]),rat.Vector([mpq(19,23), mpq(-3,19), mpq(-1,19), mpq(-13,7)]),rat.Vector([mpq(-19,23), mpq(3,19), mpq(1,19), mpq(13,7)]),rat.Vector([mpq(0,1), mpq(-3,13), mpq(6,5), mpq(3,2)]),rat.Vector([mpq(0,1), mpq(3,13), mpq(-6,5), mpq(-3,2)]),rat.Vector([mpq(0,1), mpq(15,11), mpq(1,7), mpq(27,5)]),rat.Vector([mpq(0,1), mpq(-15,11), mpq(-1,7), mpq(-27,5)]),rat.Vector([mpq(-15,16), mpq(5,3), mpq(5,26), mpq(-29,30)]),rat.Vector([mpq(15,16), mpq(-5,3), mpq(-5,26), mpq(29,30)]),rat.Vector([mpq(23,28), mpq(3,29), mpq(2,3), mpq(-2,5)]),rat.Vector([mpq(-23,28), mpq(-3,29), mpq(-2,3), mpq(2,5)]),rat.Vector([mpq(3,20), mpq(25,7), mpq(1,15), mpq(3,1)]),rat.Vector([mpq(-3,20), mpq(-25,7), mpq(-1,15), mpq(-3,1)]),rat.Vector([mpq(-22,29), mpq(-5,2), mpq(11,25), mpq(-5,1)]),rat.Vector([mpq(22,29), mpq(5,2), mpq(-11,25), mpq(5,1)]),rat.Vector([mpq(-23,1), mpq(1,1), mpq(8,1), mpq(8,7)]),rat.Vector([mpq(23,1), mpq(-1,1), mpq(-8,1), mpq(-8,7)]),rat.Vector([mpq(12,19), mpq(-11,24), mpq(19,11), mpq(-7,5)]),rat.Vector([mpq(-12,19), mpq(11,24), mpq(-19,11), mpq(7,5)]),rat.Vector([mpq(14,1), mpq(13,20), mpq(3,5), mpq(-25,9)]),rat.Vector([mpq(-14,1), mpq(-13,20), mpq(-3,5), mpq(25,9)]),rat.Vector([mpq(-7,6), mpq(4,1), mpq(-25,22), mpq(5,3)]),rat.Vector([mpq(7,6), mpq(-4,1), mpq(25,22), mpq(-5,3)]),rat.Vector([mpq(2,1), mpq(13,5), mpq(18,17), mpq(-12,25)]),rat.Vector([mpq(-2,1), mpq(-13,5), mpq(-18,17), mpq(12,25)]),rat.Vector([mpq(-6,1), mpq(-3,1), mpq(-3,2), mpq(-3,1)]),rat.Vector([mpq(6,1), mpq(3,1), mpq(3,2), mpq(3,1)]),rat.Vector([mpq(-25,4), mpq(6,19), mpq(24,29), mpq(-13,6)]),rat.Vector([mpq(25,4), mpq(-6,19), mpq(-24,29), mpq(13,6)]),rat.Vector([mpq(2,29), mpq(19,29), mpq(5,9), mpq(23,27)]),rat.Vector([mpq(-2,29), mpq(-19,29), mpq(-5,9), mpq(-23,27)]),rat.Vector([mpq(23,29), mpq(0,1), mpq(-27,14), mpq(1,1)]),rat.Vector([mpq(-23,29), mpq(0,1), mpq(27,14), mpq(-1,1)]),rat.Vector([mpq(-15,8), mpq(1,5), mpq(4,3), mpq(1,3)]),rat.Vector([mpq(15,8), mpq(-1,5), mpq(-4,3), mpq(-1,3)]),rat.Vector([mpq(-1,4), mpq(21,26), mpq(11,5), mpq(11,2)]),rat.Vector([mpq(1,4), mpq(-21,26), mpq(-11,5), mpq(-11,2)]),rat.Vector([mpq(28,11), mpq(6,5), mpq(9,5), mpq(-28,25)]),rat.Vector([mpq(-28,11), mpq(-6,5), mpq(-9,5), mpq(28,25)]),rat.Vector([mpq(11,10), mpq(-5,4), mpq(15,29), mpq(9,14)]),rat.Vector([mpq(-11,10), mpq(5,4), mpq(-15,29), mpq(-9,14)]),rat.Vector([mpq(-2,13), mpq(-15,16), mpq(-11,10), mpq(3,4)]),rat.Vector([mpq(2,13), mpq(15,16), mpq(11,10), mpq(-3,4)]),rat.Vector([mpq(-29,10), mpq(-3,10), mpq(17,4), mpq(-21,4)]),rat.Vector([mpq(29,10), mpq(3,10), mpq(-17,4), mpq(21,4)]),rat.Vector([mpq(-13,1), mpq(1,12), mpq(29,30), mpq(-15,14)]),rat.Vector([mpq(13,1), mpq(-1,12), mpq(-29,30), mpq(15,14)]),rat.Vector([mpq(3,5), mpq(-5,26), mpq(-9,5), mpq(-11,8)]),rat.Vector([mpq(-3,5), mpq(5,26), mpq(9,5), mpq(11,8)]),rat.Vector([mpq(7,15), mpq(16,11), mpq(-18,17), mpq(8,15)]),rat.Vector([mpq(-7,15), mpq(-16,11), mpq(18,17), mpq(-8,15)]),rat.Vector([mpq(7,8), mpq(25,6), mpq(-3,4), mpq(17,26)]),rat.Vector([mpq(-7,8), mpq(-25,6), mpq(3,4), mpq(-17,26)])]
    return points


bad1=[rat.Vector([mpq(-29,10), mpq(-3,10), mpq(17,4), mpq(-21,4)]),rat.Vector([mpq(7,8), mpq(25,6), mpq(-3,4), mpq(17,26)]),rat.Vector([mpq(-6,1), mpq(-3,1), mpq(-3,2), mpq(-3,1)]),rat.Vector([mpq(23,1), mpq(-1,1), mpq(-8,1), mpq(-8,7)])]

bad2=[rat.Vector([mpq(6,1), mpq(3,1), mpq(3,2), mpq(3,1)]),rat.Vector([mpq(-7,8), mpq(-25,6), mpq(3,4), mpq(-17,26)]),rat.Vector([mpq(29,10), mpq(3,10), mpq(-17,4), mpq(21,4)]),rat.Vector([mpq(-23,1), mpq(1,1), mpq(8,1), mpq(8,7)])]

badpt=rat.Vector([mpq(13,1), mpq(-1,12), mpq(-29,30), mpq(15,14)])


n11=rat.Vector([mpq(-54106257635030423679486538767180,1), mpq(291735174496489160207867837364015,1), mpq(-270188334335771729232390447992160,1), mpq(-59277926981508636440089535571450,1)])

n22=rat.Vector([mpq(54106257635030423679486538767180,1), mpq(-291735174496489160207867837364015,1), mpq(270188334335771729232390447992160,1), mpq(59277926981508636440089535571450,1)])

#sage: check_hulls(30,30,5)
#[(-4, -2, 29/20, 6/13, -24), (-23/7, -20, -1/5, 11, 2/11), (-25/8, 6/11, 10, 29/3, -13/3), (-3, -16/5, -17/4, -20/13, 8/11), (-3, 0, 14/25, 5/12, 2/11), (-28/19, 5/14, -29/7, -5/7, 6), (-7/5, 29/7, -29/8, -3/7, 15/16), (-22/17, -15, 19/25, -20/13, -5/11), (-23/19, 6, 23/20, -7/3, 12), (-22/19, 21/10, -9, 23/16, 17/29), (-15/26, -11/19, -3/7, 18, 0), (-2/7, -4/3, 13/4, 16, -22/17), (2/7, 4/3, -13/4, -16, 22/17), (15/26, 11/19, 3/7, -18, 0), (22/19, -21/10, 9, -23/16, -17/29), (23/19, -6, -23/20, 7/3, -12), (22/17, 15, -19/25, 20/13, 5/11), (7/5, -29/7, 29/8, 3/7, -15/16), (28/19, -5/14, 29/7, 5/7, -6), (3, 0, -14/25, -5/12, -2/11), (3, 16/5, 17/4, 20/13, -8/11), (25/8, -6/11, -10, -29/3, 13/3), (23/7, 20, 1/5, -11, -2/11), (4, 2, -29/20, -6/13, 24)]
#[(-4, -2, 29/20, 6/13, -24), (-23/7, -20, -1/5, 11, 2/11), (-25/8, 6/11, 10, 29/3, -13/3), (-3, -16/5, -17/4, -20/13, 8/11), (-3, 0, 14/25, 5/12, 2/11), (-28/19, 5/14, -29/7, -5/7, 6), (-7/5, 29/7, -29/8, -3/7, 15/16), (-22/17, -15, 19/25, -20/13, -5/11), (-23/19, 6, 23/20, -7/3, 12), (-22/19, 21/10, -9, 23/16, 17/29), (-15/26, -11/19, -3/7, 18, 0), (-2/7, -4/3, 13/4, 16, -22/17), (0, -29/17, -7/3, 16/19, -10), (0, 29/17, 7/3, -16/19, 10), (2/7, 4/3, -13/4, -16, 22/17), (15/26, 11/19, 3/7, -18, 0), (22/19, -21/10, 9, -23/16, -17/29), (23/19, -6, -23/20, 7/3, -12), (22/17, 15, -19/25, 20/13, 5/11), (7/5, -29/7, 29/8, 3/7, -15/16), (28/19, -5/14, 29/7, 5/7, -6), (3, 0, -14/25, -5/12, -2/11), (3, 16/5, 17/4, 20/13, -8/11), (25/8, -6/11, -10, -29/3, 13/3), (23/7, 20, 1/5, -11, -2/11), (4, 2, -29/20, -6/13, 24)]
#False


## sage says (0, -29/17, -7/3, 16/19, -10), (0, 29/17, 7/3, -16/19, 10) should not be in there



bads=[rat.Vector([mpq(-4), mpq(-2), mpq(29,20), mpq(6,13), mpq(-24)]), rat.Vector([mpq(-23,7), mpq(-20), mpq(-1,5), mpq(11), mpq(2,11)]), rat.Vector([mpq(-25,8), mpq(6,11), mpq(10), mpq(29,3), mpq(-13,3)]), rat.Vector([mpq(-3), mpq(-16,5), mpq(-17,4), mpq(-20,13), mpq(8,11)]), rat.Vector([mpq(-3), mpq(0), mpq(14,25), mpq(5,12), mpq(2,11)]), rat.Vector([mpq(-28,19), mpq(5,14), mpq(-29,7), mpq(-5,7), mpq(6)]), rat.Vector([mpq(-7,5), mpq(29,7), mpq(-29,8), mpq(-3,7), mpq(15,16)]), rat.Vector([mpq(-22,17), mpq(-15), mpq(19,25), mpq(-20,13), mpq(-5,11)]), rat.Vector([mpq(-23,19), mpq(6), mpq(23,20), mpq(-7,3), mpq(12)]), rat.Vector([mpq(-22,19), mpq(21,10), mpq(-9), mpq(23,16), mpq(17,29)]), rat.Vector([mpq(-15,26), mpq(-11,19), mpq(-3,7), mpq(18), mpq(0)]), rat.Vector([mpq(-2,7), mpq(-4,3), mpq(13,4), mpq(16), mpq(-22,17)]), rat.Vector([mpq(0), mpq(-29,17), mpq(-7,3), mpq(16,19), mpq(-10)]), rat.Vector([mpq(0), mpq(29,17), mpq(7,3), mpq(-16,19), mpq(10)]), rat.Vector([mpq(2,7), mpq(4,3), mpq(-13,4), mpq(-16), mpq(22,17)]), rat.Vector([mpq(15,26), mpq(11,19), mpq(3,7), mpq(-18), mpq(0)]), rat.Vector([mpq(22,19), mpq(-21,10), mpq(9), mpq(-23,16), mpq(-17,29)]), rat.Vector([mpq(23,19), mpq(-6), mpq(-23,20), mpq(7,3), mpq(-12)]), rat.Vector([mpq(22,17), mpq(15), mpq(-19,25), mpq(20,13), mpq(5,11)]), rat.Vector([mpq(7,5), mpq(-29,7), mpq(29,8), mpq(3,7), mpq(-15,16)]), rat.Vector([mpq(28,19), mpq(-5,14), mpq(29,7), mpq(5,7), mpq(-6)]), rat.Vector([mpq(3), mpq(0), mpq(-14,25), mpq(-5,12), mpq(-2,11)]), rat.Vector([mpq(3), mpq(16,5), mpq(17,4), mpq(20,13), mpq(-8,11)]), rat.Vector([mpq(25,8), mpq(-6,11), mpq(-10), mpq(-29,3), mpq(13,3)]), rat.Vector([mpq(23,7), mpq(20), mpq(1,5), mpq(-11), mpq(-2,11)]), rat.Vector([mpq(4), mpq(2), mpq(-29,20), mpq(-6,13), mpq(24)])]