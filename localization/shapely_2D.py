import geometry as gx
from math import radians, cos, sin, asin, sqrt, pi
import shapely as shx
from shapely.geometry import Polygon, Point, MultiPolygon

sdr = 1  # the resolution of polygons in degree (Sampling Degree Resolution)

''' This file returns the Shapely polygon representation
	of the small circles 
	on the earth based on anchor-target distance measurements '''

res = 1


def d2rad(d):  # Distance in meters to radian
    return d / (gx.E.R)


def boundCheck(p, r):  # Checking if the measurements cover Poles
    x = p.x
    y = p.y
    bc = ['0'] * 2
    if p.y + r > 90:  # North Pole
        bc[0] = '1'
    if p.y - r < -90:  # South Pole
        bc[1] = '1'
    return ''.join(bc)


def order(v, w, mode=0):  # sorting vector v and w pivoting on v
    a = zip(v, w)
    if mode == 0:
        a.sort()
    else:
        a.sort(reverse=True)
    l = zip(*a)
    v = list(l[0])
    w = list(l[1])
    return [v, w]


def polygonize(c, d):
    # (!) small circle that has both north Pole and southPole
    # c is the list of server locations (in gx.Point format)
    # d is the list of distances(Great Circle Distance)
    U = Polygon([(-180, -90), (180, -90), (180, 90), (-180, 90)])
    l = len(c)
    P = [U] * l
    for i in range(l):
        pc = c[i]
        r = d[i]
        deg = d2rad(r) * 180.0 / pi
        if deg >= 180 - res:  # This measurement covers whole earth
            continue
        bc = boundCheck(pc, deg)
        stld = gx.E.gcd2l(r)
        S = gx.Sphere(gx.E.map(pc, True), stld)
        sc = gx.E.small_circle(S)
        vert1 = [gx.E.map(xx) for xx in sc]
        vert = [(xx.x, xx.y) for xx in vert1]
        if bc == '10' or bc == '01':
            lon, lat = zip(*vert)
            lon, lat = order(lon, lat)
            c1 = (180, lat[-1])
            c2 = (180, 90)
            c3 = (-180, 90)
            c4 = (-180, lat[0])
            cur = Polygon(zip(lon, lat) + [c1, c2, c3, c4])
        else:
            lon, lat = zip(*vert)
            pivot = lon[0]
            ii = -1
            j = -1
            fg = True
            for k, xx in enumerate(lon[1:]):
                if fg:
                    if abs(xx - pivot) > min(90, 10 * sdr):  # Passing (+/-)180 meridian
                        ii = k
                        fg = False
                    pivot = xx
                else:
                    if abs(xx - pivot) > min(90, 10 * sdr):  # Passing (+/-)180 meridian for the second time
                        j = k
                        break
                    pivot = xx
                    # Note that a small circle can not pass meridian more than twice
            if ii == -1 and j == -1:
                cur1 = Polygon(zip(lon, lat))
                if bc == '11':  # Containing both Poles
                    cur = U - cur1
                else:
                    cur = cur1
            elif ii == -1 or j == -1:
                raise gx.geoError('Algorithmic bug in polygonizing for earth model')
            else:
                lon1 = lon[ii + 1:j + 1]
                lat1 = lat[ii + 1:j + 1]
                lon2 = lon[j + 1:] + lon[0:ii + 1]
                lat2 = lat[j + 1:] + lat[0:ii + 1]
                if lon1[0] > lon2[0]:
                    P1 = Polygon([(180, lat1[0])] + zip(lon1, lat1) + [(180, lat1[-1])])
                    P2 = Polygon([(-180, lat2[0])] + zip(lon2, lat2) + [(-180, lat2[-1])])
                else:
                    P1 = Polygon([(-180, lat1[0])] + zip(lon1, lat1) + [(-180, lat1[-1])])
                    P2 = Polygon([(180, lat2[0])] + zip(lon2, lat2) + [(180, lat2[-1])])
                cur1 = P1.union(P2)
                if bc == '11':  # Containing both Poles
                    cur = U - cur1
                else:
                    cur = cur1
        P[i] = cur
    return P
