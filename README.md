# CHull

Given a collection of points C in R^n, computes the convex hull of the union C U C', where C' consists of the negatives of the points in C (i.e., the reflection through the origin of C). This is (almost) an implementation of the Quickhull algorithm of Barber, Dobkin, and Huhdanpaa.

The result is a polyhedron symmetric about the origin, which is the smallest centrally symmetric convex polyhedron containing the points in C. 

CHull was developed for use in TNorm, which only requires computation of centrally symmetric polyhedra. It likely wouldn't be too hard to adapt this code so that it computes the true convex hull of C, and I may do that eventually.
