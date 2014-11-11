===========
Localization
===========

Localization package provides tools for multilateration and triangulation
in '2D','3D' and on earth surface. 
The current model of the earth, supported by the package, is called 'Earth1'. 
Earth1 models the earth as an ideal sphere having radius of 6378.1 kilometers. Typical usage of the package is:


    import localization as lx

To initilaize new localization Project use:

P=lx.Project(mode=<mode>,solver=<solver>)

Currently three modes are supported:
1-2D
2-3D
3-Earth1

Also three solvers can be utilized:
1-LSE for least square error
2-LSE_GC for least square error with geometric constraints. Geometric constraints force the solutions to be in the intersection areas of all multilateration circles.
3- CCA for centroid method, i.e., the solution will be the centroid of the intersection area. If no common intersection area exist, the area with maximum overlap is used.

To add anchors to the project use:

P.add_anchor(<name>,<loc>)

where name denote user provided label of the anchor and <loc> is the location of the anchor provided in tuple, e.g., (120,60).

To add target use:

t,label=P.add_target()

t is the target object and label is the package provided label for the target.

Distance measurements must be added to target object like:

t.add_measure(<anchore_lable>,<measured_distance>)

Finally running P.solve() will locate all targets. You can access the estimated location of the target t by t.loc.
t.loc is a point object. Point object B has "x","y","z" coordinates available by B.x, B.y, B.z respectively.

Install the package by:
pip install localization

contact: kamal.shadi85@gmail.com


An Example
=========

Basic
-------------

To be completed:

1. README

2. DOC
