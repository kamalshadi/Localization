import geometry as gx
from shapely.geometry import Polygon,Point,MultiPolygon

""" 2D shapely polygon models for circles"""

def polygonize(c,d):
	l=len(c)
	P=[Point(0, 0).buffer(1.0)]*l
	for i in range(l):
		pc=c[i]
		r=d[i]
		P[i]=Point(pc.x, pc.y).buffer(r)
	return P
