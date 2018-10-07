# Localization

### Introduction
Localization package provides tools for multilateration and triangulation in `2D`,`3D` and on `earth surface`.
The current model of the earth, supported by the package, is called `Earth1`.
Earth1 models the earth as an ideal sphere having radius of 6378.1 kilometers.

### Installation
The package has been tested with python2.7 and python3.6. Use pip to install the package:
```
pip install localization
```
### Usage

Typical usage of the package is:

```python
import localization as lx
```
To initilaize new localization Project use:

```python
P=lx.Project(mode=<mode>,solver=<solver>)
```

Currently three modes are supported:
 - 2D
 - 3D
 - Earth1

Also three solvers can be used:
 - LSE for least square error
 - LSE_GC for least square error with geometric constraints. Geometric constraints force the solutions to be in the intersection areas of all multilateration circles.
 - CCA for centroid method, i.e., the solution will be the centroid of the intersection area. If no common intersection area exist, the area with maximum overlap is used.

To add anchors to the project use:

```python
P.add_anchor(<name>,<loc>)
```

where name denote user provided label of the anchor and <loc> is the location of the anchor provided in tuple, e.g., (120,60).

To add target use:

```python
t,label=P.add_target()
```

t is the target object and label is the package provided label for the target.

Distance measurements must be added to target object like:

```python
t.add_measure(<anchore_lable>,<measured_distance>)
```

Finally running P.solve() will locate all targets. You can access the estimated location of the target t by t.loc.
t.loc is a point object. Point object `B` has `x`,`y`,`z` coordinates available by `B.x`, `B.y`, `B.z` respectively.

Here is a sample use of the package for three anchors and one target scenario:

```python
import localization as lx

P=lx.Project(mode='2D',solver='LSE')


P.add_anchor('anchore_A',(0,100))
P.add_anchor('anchore_B',(100,100))
P.add_anchor('anchore_C',(100,0))

t,label=P.add_target()

t.add_measure('anchore_A',50)
t.add_measure('anchore_B',50)
t.add_measure('anchore_C',50)

P.solve()

# Then the target location is:

print(t.loc)
```

---
contact: __kamal.shadi85@gmail.com__
