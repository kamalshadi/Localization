import shapely as shx
from shapely.geometry import Polygon, Point, MultiPolygon
import numpy as num
from itertools import combinations

"""The python functions in this files return the area of maximim intersection
 for a list of shapely polygons"""


class geoError(Exception):
    def __init__(self, value):
        self.tag = value

    def __str__(self):
        return repr(self.tag)


def intersection_matrix(P):
    l = len(P)
    I = num.eye(l)
    for i in range(l):
        for j in range(i + 1, l):
            if P[i].intersects(P[j]):
                I[i][j] = 1
                I[j][i] = 1
    return I


def canInd(P, ni):
    l = len(P)
    ind = range(l)
    if ni < 2:
        return [[xx] for xx in ind]
    if ni >= l:
        return [ind]
    im = intersection_matrix(P)
    can = []
    for w in combinations(ind, ni):
        fg = True
        for i in w:
            for j in w:
                if im[i, j] == 0:
                    fg = False
                    break
            if not fg:
                break
        if fg:
            can.append(list(w))
    return can


def checkCan(P, can):
    l = len(can)
    his = float('inf')
    if len(can) < 1:
        return
    passFlag = False
    for k, w in enumerate(can):
        C = P[w[0]]
        fg = False
        for j in w[1:]:
            if C.intersects(P[j]):
                try:
                    C = C.intersection(P[j])
                except:
                    print('(!)Warning: in checkCan')
            else:
                fg = True
                break
        if fg:
            continue
        passFlag = True
        a = C.area
        if a < his:
            his = a
            temp = C
    if passFlag:
        return [temp]
    else:
        return []


def maxPol(P):
    # This function find the area of maximum intersection with minimum area tie-breaker
    l = len(P)
    ni = l
    while ni > 0:
        can = canInd(P, ni)
        if len(can) == 0:
            ni = ni - 1
            continue
        Pl = checkCan(P, can)
        if len(Pl) > 0:
            P = Pl[0]
            break
        else:
            ni = ni - 1
            continue
    if ni > 0:
        return (P, ni)
    raise geoError('UnKnown')
