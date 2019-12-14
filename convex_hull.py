

from itertools import combinations
import rat_matrices as rat
from gmpy2 import mpq


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

class Polyhedron():
    def __init__(self,dim,facets=[]):
        self.dim = dim
        self._facets = {dim-1:tuple(facets)}
        if self.dim > 1:
            for d in range(0,dim-1):
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
            l = rat.LCM([c.denominator for c in normal])
            normal = rat.Vector([c*l for c in normal])

            return normal


    def __repr__(self):
        verts = ''
        for v in self.vertices():
            verts=verts+str(v)+','

        return "Convex hull of vertices: "+verts[:-1]

    def __hash__(self):
        return hash(self.vertices())

    def __eq__(self,other):
        return ( self.__class__ == other.__class__ ) and ( self.vertices() == other.vertices() )

    def __ne__(self,other):
        return ( self.__class__ != other.__class__ ) or ( self.vertices() != other.vertices() )

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

    def negative(self, facet):
        d = facet.dim
        for F in self.facets(d):
            if set([v.coords() for v in F.vertices()]) == set([-v.coords() for v in facet.vertices()]):
                return F
            else:
                continue


class Vertex():
    def __init__(self,coords):
        self._coords = rat.Vector(coords)
        self.dim = 0

    def __repr__(self):
        return str(self._coords)

    def __hash__(self):
        return hash(self._coords)

    def __eq__(self,other):
        return ( self.__class__ == other.__class__ ) and self._coords == other._coords

    def __ne__(self,other):
        return ( self.__class__ != other.__class__ ) or self._coords != other._coords

    def coords(self):
        return self._coords




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
        dist = mpq(sum([v.ns() for v in face.vertices()]),len(face.vertices()))
        subfaces.append((dist,subface))
    closest_subface = sorted(subfaces)[0]
    subfaces = [subface for _,subface in subfaces[1:]] # get rid of the close face, and forget the distances
    for subsubface in closest_subface.facets(closest_subface.dim-1):
        new_face = cone_on_face(subsubface,vertex)
        subfaces.append(new_face)
    return Polyhedron(subfaces)


def quick_hull(vector_points, dim):
    points = vector_points
    # choose dim+1 points for the initial hull. We might as well pick the ones of largest norm.
    v_max = sorted([(v.ns(),v) for v in points],key=lambda x:x[0])[-1][1]
    B_points = [v_max]
#    points = [p for p in points if p != rat.Vector([0,0,0])]
#    v = points[0]
    B_points = [v_max]
#    points.remove(v); points.remove(-v)
    while len(B_points) < dim:
        orthobasis = rat.gram_schmidt([v for v in B_points])
        perp_norms = dict([((v-rat.projection(v,orthobasis)).ns(),v) for v in points])
        max_perp = perp_norms[max(perp_norms.keys())]
        #max_perp = lin_ind[0]

        B_points.append(max_perp)
    F = simplex(B_points,dim-1)
    facets = []
    for set1,set2 in partitions([v.coords() for v in F.vertices()]):
        facets.append(simplex([-v for v in set1]+[v for v in set2],dim-1))
        facets.append(simplex([-v for v in set2]+[v for v in set1],dim-1))
    hull = Polyhedron(dim,facets)
    above = {} # this dictionary will associate with each face the vertices that are above it 
                 # (i.e., those having positive inner product with the outward normal of the face)
    for F in hull.facets(dim-1):
        above[F] = outside_of(F,points) 
    points = set([p for F in hull.facets(hull.dim-1) for p,strictly_above in above[F] if strictly_above]) # we can get rid of points not above at least one face
    while len(points) > 0: # as long as there are still points outside of the hull, expand faces to points that lie above them
        for F in hull.facets(dim-1):
            if len(above[F]) > 0:
                p = above[F][-1][0] # pick the point furthest from F (i.e., the last in the list)
                break
        hull,points,above = expand(hull,p,points,above)


    return hull

def expand(hull,point,points,above):
    faces = set(hull.facets(hull.dim-1)) # all faces of current hull
    visible = set([(face,i) for face in faces for i in [0,1] if (point,i) in above[face]]) # faces visible from point
    horizons = set([]) # facets of dimension dim-2 that are on the boundary of the union of visible faces as seen from point
    perfect_horizons = set([]) # facets of dimension dim-1 that are co-planar with point
    for face, strictly_above in visible:
        if not strictly_above:
            perfect_horizons.add(face)
        else:
            for subface in face.facets(hull.dim-2):
                on_faces = [f for f,_ in visible if hull.is_subfacet(subface,f)]
                if len(on_faces) == 1: # only add to horizons if it is not shared by multiple faces (so it's not in the    interior)
                    horizons.add(subface)
    #points_above_visible = set([])
    for face, _ in visible:
        neg_face = hull.negative(face)
        faces.discard(face) #remove visible faces as they are now covered by the new faces we'll create below
        faces.discard(neg_face)
        #for p, _ in above[face]:
        #    points_above_visible.add(p)
        del above[face] # also remove these from the above dictionary
        del above[neg_face] # also remove these from the above dictionary
    #points_above_visible.discard(point)
    for h in horizons:
        new_face = cone_on_face(h,Vertex(point)) # new faces come from coning over the horizons to point
        neg_new_face = cone_on_face(hull.negative(h),Vertex(-point)) # new faces come from coning over the horizons to point
        above[new_face] = outside_of(new_face,points) # populate above[new_face] with points that lie above  new_face
        above[neg_new_face] = [(-p,strictly_above) for p,strictly_above in above[new_face]]
        faces.add(new_face) # add to faces each new face
        faces.add(neg_new_face)
    for F in perfect_horizons:
        new_face = extend_face(F, Vertex(point))
        neg_new_face = extend_face(hull.negative(F), Vertex(-point))
        above[new_face] = outside_of(new_face,points) # populate above[new_face] with points that lie above  new_face
        above[neg_new_face] = [(-p,strictly_above) for p,strictly_above in above[new_face]]
        faces.add(new_face) # add to faces each new face
        faces.add(neg_new_face)
        
    hull = Polyhedron(hull.dim,faces) # construct the new hull
    points = set([p for F in hull.facets(hull.dim-1) for p,strictly_above in above[F] if strictly_above])
    return hull, points, above


def simplex(vector_points,dim):
    points = vector_points
    assert len(points) == dim+1
    if dim == 1:
        P = Polyhedron(1,[Vertex(points[0].as_tuple()),Vertex(points[1].as_tuple())])
        return P
    else:
        facets = [simplex(f,dim-1) for f in list(combinations(points,dim))]
        P = Polyhedron(dim,facets)
        return P


def outside_of(face,points):
    outside = []
    basepoint = face.vertices()[0].coords() # pick a vertex to translate everything to the origin by
    for p in points:
        if p not in [v.coords() for v in face.vertices()]:
            dot = (p-basepoint)*face.normal
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








