

from itertools import combinations
import chull.rat_mats_pari as rat
import sys
from cypari import pari
import cypari

def convert_to_QQ(rational):
    try:
        return rat.QQ(rational.numerator(), rational.denominator())
    except:
        return rat.QQ(rational.numerator, rational.denominator)


def convert_to_vec(tuple_or_list):
    p = tuple_or_list
    vector = pari.vector(len(p), [convert_to_QQ(c) for c in p])
    return vector

class cached_property(object):
    """
    Descriptor (non-data) for building an attribute on-demand on first use.
    """
    def __init__(self, factory):
        """
        <factory> is called such: factory(instance) to build the attribute.
        """
        self._attr_name = factory.__name__
        self._factory = factory

    def __get__(self, instance, owner):
        # Build the attribute.
        attr = self._factory(instance)

        # Cache the value; hide ourselves.
        setattr(instance, self._attr_name, attr)

        return attr

def row_flatten(list_of_lists):
    return [list_of_lists[i][j] for i in range(len(list_of_lists)) for j in range(len(list_of_lists[0]))]
def col_flatten(list_of_lists):
    return [list_of_lists[i][j] for j in range(len(list_of_lists[0])) for i in range(len(list_of_lists))]

class Polyhedron():
    def __init__(self,dim,facets=[]):
        self.dim = dim
        if dim == 1:
            self._facets = {dim-1:tuple(facets)}
        else:
            self._facets = {dim-1:tuple(facets)}
        if self.dim > 1:
            for d in range(0,dim-1):
                if d == 0:
                    self._facets[d] = tuple(set([dface for f in facets for dface in f.facets(d)]))
                else:
                    self._facets[d] = tuple(set([dface for f in facets for dface in f.facets(d)]))
        self.ambient_dim = len(self.vertices()[0].coords())
  
    @cached_property
    def normal(self): # if this polyhedron is co-dim 1 in the ambient space, this give the normal to the polyhedron (not nec. unit length)
        if not self.ambient_dim - self.dim == 1:
            return 'undefined'
        else:

            basepoint = self.vertices()[0] # pick a vertex to translate everything to the origin by
            pre_basis = [v.coords() for v in self.neighbors(basepoint)]
            basepoint = basepoint.coords()
            basis = [v-basepoint for v in pre_basis][:self.dim]
            orthobasis = rat.gram_schmidt(basis)
            normal = basepoint - rat.projection(basepoint,orthobasis)
            l = rat.LCM([c.denominator() for c in normal])
            normal = pari.vector(len(normal),([c*l for c in normal]))

            return normal


    def __repr__(self):
        verts = ''
        for v in self.vertices():
            verts=verts+str(v)+','

        return "Convex hull of vertices: "+verts[:-1]

    def __hash__(self):
        return hash(tuple(sorted(self.vertices())))

    def __eq__(self,other):
        return ( self.__class__ == other.__class__ ) and ( set(self.vertices()) == set(other.vertices()) )

    def __ne__(self,other):
        return ( self.__class__ != other.__class__ ) or ( set(self.vertices()) != set(other.vertices()) )

    def facets(self,d):
            return self._facets[d]

    def vertices(self):
            return self._facets[0]

    def is_subfacet(self,f1,f2):
        if f1 in f2.facets(f1.dim):
                return True
        else:
            return False


    def neighbors(self,vertex):
        neighbors = []
        for e in self.facets(1):
            if self.is_subfacet(vertex,e):
                for v in e.vertices():
                    if v != vertex:
                        neighbors.append(v)
        return neighbors

    @cached_property
    def negative(self):
        if self.dim == 1:
            return Polyhedron(1,[Vertex(-self.vertices()[0].coords()),Vertex(-self.vertices()[1].coords())])
        else:
            return Polyhedron(self.dim,[f.negative for f in self.facets(self.dim-1)])


class Vertex():
    def __init__(self,coords):
        self._coords = coords # should be a pari.vector
        self.dim = 0

    def __repr__(self):
        return str(self._coords)

    def __hash__(self):
        return hash(self._coords)

    def __eq__(self,other):
        return ( self.__class__ == other.__class__ ) and self._coords == other._coords

    def __ne__(self,other):
        return ( self.__class__ != other.__class__ ) or self._coords != other._coords

    def __sub__(self,other):
        return Vertex( self._coords - other._coords )

    def __add__(self,other):
        return Vertex( self._coords + other._coords )

    def __le__(self,other):
        return tuple(self._coords) <= tuple(other._coords)
    def __lt__(self,other):
        return tuple(self._coords) < tuple(other._coords)    
    def __ge__(self,other):
        return tuple(self._coords) >= tuple(other._coords)   
    def __gt__(self,other):
        return tuple(self._coords) > tuple(other._coords)

    def coords(self):
        return self._coords


def simplex(vector_points,dim):
    points = vector_points
    assert len(points) == dim+1
    if dim == 1:
        P = Polyhedron(1,[Vertex(points[0]),Vertex(points[1])])
        return P
    else:
        facets = [simplex(f,dim-1) for f in list(combinations(points,dim))]
        P = Polyhedron(dim,facets)
        return P


def cone_on_face(face,vertex):
    if face.dim == 0:
        return Polyhedron(1,[face,vertex])
    else:
        facets = [face]
        for f in face.facets(face.dim-1):
            facets.append(cone_on_face(f,vertex))
        return Polyhedron(face.dim+1,facets)

def extend_face(face,vertex):
    subfaces = []
    for subface in face.facets(face.dim-1):
        dist = rat.QQ(sum([ns(v.coords()-vertex.coords()) for v in face.vertices()]),len(face.vertices()))
        subfaces.append((dist,subface))
    closest_subface = sorted(subfaces,key=lambda x: x[0])[0][1]
    subfaces = [subface for _,subface in subfaces[1:]] # get rid of the close face, and forget the distances
    for subsubface in closest_subface.facets(closest_subface.dim-1):
        new_face = cone_on_face(subsubface,vertex)
        subfaces.append(new_face)
    return Polyhedron(face.dim,subfaces)

def ns(v):
    return (v*v.mattranspose())[0]

def quick_hull(vector_points, dim):
    points = vector_points

    if not isinstance(points[0],cypari._pari.Gen):
        points = [convert_to_vec(point) for point in points]


    # choose dim+1 points for the initial hull. We might as well pick the ones of largest norm.
    v_max = sorted([(ns(v),v) for v in points],key=lambda x:x[0])[-1][1]
    B_points = [v_max]
#    points = [p for p in points if p != rat.Vector([0,0,0])]
#    v = points[0]
    B_points = [v_max]
#    points.remove(v); points.remove(-v)
    while len(B_points) < dim:
        orthobasis = rat.gram_schmidt([v for v in B_points])
        perp_norms = dict([(ns((v-rat.projection(v,orthobasis))),v) for v in points])
        max_perp = perp_norms[max(perp_norms.keys())]
        #max_perp = lin_ind[0]

        B_points.append(max_perp)
    F = simplex(B_points,dim-1)
    faces = set([])
    for set1,set2 in partitions([v.coords() for v in F.vertices()]):
        list1 = list(set1); list2 = list(set2)
        faces.add(simplex([-v for v in list1]+[v for v in list2],dim-1))
        faces.add(simplex([v for v in list1]+[-v for v in list2],dim-1))
    above = {} # this dictionary will associate with each face the vertices that are above it 
                 # (i.e., those having positive inner product with the outward normal of the face)
    for F in faces:
        above[F] = outside_of(F,points) 
    points = set([p for F in faces for p,strictly_above in above[F] if strictly_above]) # we can get rid of points not above at least one face
    while len(points) > 0: # as long as there are still points outside of the hull, expand faces to points that lie above them
        for F in faces:
            if len(above[F]) > 0:
                p = above[F][-1][0] # pick the point furthest from F (i.e., the last in the list)
                break
        faces,points,above = expand(faces,p,points,above)

    hull = Polyhedron(dim,list(faces))
    return hull

def expand(faces,point,points,above):
    visible = set([(face,i) for face in faces for i in [0,1] if (point,i) in above[face]]) # faces visible from point
    horizons = set([]) # facets of dimension dim-2 that are on the boundary of the union of visible faces as seen from point
    perfect_horizons = set([]) # facets of dimension dim-1 that are co-planar with point
    subfaces = set([])
    for face, strictly_above in visible:
        if not strictly_above:
            for subface in face.facets(face.dim-1):
                if subface not in subfaces:
                    subfaces.add(subface)
                else:
                    horizons.discard(subface)
            perfect_horizons.add(face)
        else:
            for subface in face.facets(face.dim-1):
                if subface not in subfaces:
                    subfaces.add(subface)
                    horizons.add(subface)
                else:
                    horizons.discard(subface)
    #points_above_visible = set([])
    for face, _ in visible:
        neg_face = face.negative
        faces.remove(face) #remove visible faces as they are now covered by the new faces we'll create below
        faces.remove(neg_face)
        #for p, _ in above[face]:
        #    points_above_visible.add(p)
        del above[face] # also remove these from the above dictionary
        del above[neg_face] # also remove these from the above dictionary
    #points_above_visible.discard(point)
    for h in horizons:
        new_face = cone_on_face(h,Vertex(point)) # new faces come from coning over the horizons to point
        neg_new_face = cone_on_face(h.negative,Vertex(-point))
        #neg_new_face = new_face.negative
        above[new_face] = outside_of(new_face,points) # populate above[new_face] with points that lie above  new_face
        above[neg_new_face] = [(-p,strictly_above) for p,strictly_above in above[new_face]]
        faces.add(new_face) # add to faces each new face
        faces.add(neg_new_face)
    for F in perfect_horizons:
        new_face = extend_face(F, Vertex(point))
        neg_new_face = new_face.negative
        above[new_face] = outside_of(new_face,points) # populate above[new_face] with points that lie above  new_face
        above[neg_new_face] = [(-p,strictly_above) for p,strictly_above in above[new_face]]
        faces.add(new_face) # add to faces each new face
        faces.add(neg_new_face)
        
    points = set([p for F in faces for p,strictly_above in above[F] if strictly_above])
    return faces, points, above




def outside_of(face,points):
    outside = []
    basepoint = face.vertices()[0].coords() # pick a vertex to translate everything to the origin by
    for p in points:
        if p not in [v.coords() for v in face.vertices()]:
            dot = rat.dot((p-basepoint),face.normal)
            if dot >= 0:
                outside.append((p,dot))

#    outside.sort(key=lambda x: x[1])
    return [(p,int(bool(dot))) for p,dot in outside]


def partitions(collection):
    collection = list(collection)
    if len(collection) == 1:
        return [[ collection,[] ]]

    parts = []
    first = collection[0]
    for [set1, set2] in partitions(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        parts.append([set1,set2+[first]])
        parts.append([set1+[first],set2])
    return parts








