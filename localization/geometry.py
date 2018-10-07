import math as math
import numpy as num

# ~ from myBasic import pickColor

# Global
res = 0.0001
s_res = math.pi / 180.0


def order(v, w, mode=0):
    a = zip(v, w)
    if mode == 0:
        a.sort()
    else:
        a.sort(reverse=True)
    l = zip(*a)
    v = list(l[0])
    w = list(l[1])
    return [v, w]


def cvHull(pL):
    if len(pL) < 3:
        return pL
    from scipy.spatial import ConvexHull
    from scipy.spatial.qhull import QhullError
    S = [xx.std() for xx in pL]
    try:
        hull = ConvexHull(S)
    except QhullError:
        xL = [xx.x for xx in pL]
        xx, pL = order(xL, pL)
        return [pL[0], pL[-1]]
    P = [pL[i] for i in hull.vertices]
    P.append(pL[hull.vertices[0]])
    return P


def solve2p(a, b, c):
    d = float(b ** 2 - 4 * a * c)
    if d < 0:
        return []
    else:
        d = d ** .5
        x1 = (-b + d) / (2 * a)
        x2 = (-b - d) / (2 * a)
    return [x1, x2]


def angleMap(a):
    a1 = a % (2 * math.pi)
    if a1 > math.pi:
        return -2 * math.pi + a1
    else:
        return a1


def uni_obj(ol):
    l = len(ol)
    if l < 2:
        return ol
    out = [ol[0]]
    for i in range(1, l):
        cO = ol[i]
        fg = True
        for w in out:
            if w == cO:
                fg = False
        if fg:
            out.append(cO)
    return out


def circle_intersect(c1, c2):
    d = c1.c.dist(c2.c)
    x1 = (c2.r ** 2 - c1.r ** 2 - d ** 2) / (-2 * d)
    y1 = (c1.r ** 2 - x1 ** 2) ** .5
    y2 = -(c1.r ** 2 - x1 ** 2) ** .5
    v = vec(c1.c, c2.c)
    a = v.angle()
    p1 = point(x1, y1)
    p2 = point(x1, y2)
    p1 = p1.rot(a) + c1.c
    p2 = p2.rot(a) + c1.c
    return [p1, p2]


def circles_intersect(sA):
    # unique touching point between 3 circles
    sA = uni_obj(sA)
    l = len(sA)
    if l < 3:
        raise geoError('InputError')
    # Checking disjoint case
    for i in range(0, l - 1):
        for j in range(i + 1, l):
            if (sA[i].c.dist(sA[j].c) > (sA[i].r + sA[j].r)):
                raise geoError('Disjoint')
    p1, p2 = circle_intersect(sA[0], sA[1])
    p3, p4 = circle_intersect(sA[2], sA[1])
    if p1.dist(p3) < res or p1.dist(p4) < res:
        return p1
    else:
        return p2


def drawC(cA, ax):
    xmin = float('inf')
    xmax = -xmin
    ymin = float('inf')
    ymax = -ymin
    for q, w in enumerate(cA):
        t1 = w.c.x - w.r
        t2 = w.c.x + w.r
        if t1 < xmin:
            xmin = t1
        if t2 > xmax:
            xmax = t2
        t1 = w.c.y - w.r
        t2 = w.c.y + w.r
        if t1 < ymin:
            ymin = t1
        if t2 > ymax:
            ymax = t2
        circ = pl.Circle((w.c.x, w.c.y), radius=w.r, fc='none', lw=2, ec=pickColor(q))
        ax.add_artist(circ)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])


def order(v, w, mode=0):
    a = zip(v, w)
    if mode == 0:
        a.sort()
    else:
        a.sort(reverse=True)
    l = zip(*a)
    v = list(l[0])
    w = list(l[1])
    return [v, w]


class geoError(Exception):
    def __init__(self, value):
        self.tag = value

    def __str__(self):
        return repr(self.tag)


class point:
    def __init__(self, *argv):
        l = len(argv)
        if l == 1:
            self.dim = len(argv[0])
        else:
            self.dim = l
        if l == 1:
            self.x = argv[0][0]
            self.y = argv[0][1]
            try:
                self.z = argv[0][2]
            except IndexError:
                self.z = 0.0
        else:
            if l == 2:
                z = 0
            elif l == 3:
                z = argv[2]
            else:
                raise geoError('Input')
            self.x = float(argv[0])
            self.y = float(argv[1])
            self.z = z

    def __str__(self):
        return 'p(' + str(self.x) + ',' + str(self.y) + ',' + str(self.z) + ')'

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __sub__(self, other):
        if isinstance(other, point):
            tx = self.x - other.x
            ty = self.y - other.y
            tz = self.z - other.z
        else:
            tx = self.x - other.dx
            ty = self.y - other.dy
            tz = self.z - other.dz
        return point(tx, ty, tz)

    def __add__(self, other):
        if isinstance(other, point):
            tx = self.x + other.x
            ty = self.y + other.y
            tz = self.z + other.z
        else:
            tx = self.x + other.dx
            ty = self.y + other.dy
            tz = self.z + other.dz
        return point(tx, ty, tz)

    def __mul__(self, other):
        return point(other * self.x, other * self.y, other * self.z)

    def __rmul__(self, other):
        return point(other * self.x, other * self.y, other * self.z)

    def __div__(self, other):
        return point(self.x / other, self.y / other, self.z / other)

    def __neg__(self):
        x = -self.x
        y = -self.y
        z = -self.z
        return point(x, y, z)

    def area(self):
        return 0.0

    def dist(self, other):
        return ((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2) ** 0.5

    def std(self):
        if self.dim == 2:
            return [self.x, self.y]
        return [self.x, self.y, self.z]

    def c2s(self):
        R = self.dist(point(0, 0, 0))
        lg = math.atan(self.y / self.x)
        lat = acos(self.z / R)
        return (lg, lat, R)

    # ~ def transform(self,p1,p2):
    # ~ if isinstance(p2,point):
    # ~ v=vec(p1,p2)
    # ~ rot=v.angle()
    # ~ return self.transform(p1,rot)
    # ~ else:
    # ~ temp=self-p1
    # ~ rot=p2
    # ~ px=math.cos(rot)*temp.x+math.sin(rot)*temp.y
    # ~ py=-math.sin(rot)*temp.x+math.cos(rot)*temp.y
    # ~ return point(px,py)
    def transform(self, p, rot):
        px = math.cos(rot) * self.x + math.sin(rot) * self.y
        py = -math.sin(rot) * self.x + math.cos(rot) * self.y
        p_t = point(px, py)
        return p_t - p

    def rot(self, a):
        px = math.cos(a) * self.x - math.sin(a) * self.y
        py = math.sin(a) * self.x + math.cos(a) * self.y
        return point(px, py)

    def angle(self, p):
        v = vec(self, p)
        return v.angle()


class vec:
    def __init__(self, *argv):
        if len(argv) >= 2:
            if isinstance(argv[0], point) and isinstance(argv[1], point):
                p = argv[1] - argv[0]
                self.dx = float(p.x)
                self.dy = float(p.y)
                self.dz = float(p.z)
            else:
                self.dx = float(argv[0])
                self.dy = float(argv[1])
                try:
                    self.dz = float(argv[2])
                except IndexError:
                    self.dz = 0.0
        else:
            p = argv[0]
            self.dx = float(p.x)
            self.dy = float(p.y)
            self.dz = float(p.z)

    def __str__(self):
        return 'vec(' + str(self.dx) + ',' + str(self.dy) + ',' + str(self.dz) + ')'

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __sub__(self, other):
        if isinstance(other, vec):
            tx = self.dx - other.dx
            ty = self.dy - other.dy
            tz = self.dz - other.dz
            return vec(tx, ty, tz)
        else:
            tx = self.dx - other.x
            ty = self.dy - other.y
            tz = self.dz - other.z
            return point(tx, ty, tz)

    def __add__(self, other):
        if isinstance(other, vec):
            tx = self.dx + other.dx
            ty = self.dy + other.dy
            tz = self.dz + other.dz
            return vec(tx, ty, tz)
        else:
            tx = self.dx + other.x
            ty = self.dy + other.y
            tz = self.dz + other.z
            return point(tx, ty, tz)

    def __mul__(self, other):
        if isinstance(other, vec):
            return self.cross(other)
        return vec(other * self.dx, other * self.dy, other * self.dz)

    def __rmul__(self, other):
        if isinstance(other, vec):
            return other.cross(self)
        return vec(other * self.dx, other * self.dy, other * self.dz)

    def __div__(self, other):
        return vec(self.dx / other, self.dy / other, self.dz / other)

    def dot(self, v):
        return self.dx * v.dx + self.dy * v.dy + self.dz * v.dz

    def cross(self, v):
        z = self.dx * v.dy - self.dy * v.dx
        x = self.dy * v.dz - self.dz * v.dy
        y = self.dz * v.dx - self.dx * v.dz
        return vec(x, y, z)

    def mag(self):
        return (self.dx ** 2 + self.dy ** 2 + self.dz ** 2) ** 0.5

    def angle(self, *args):
        x = self.dx
        y = self.dy
        z = self.dz

        if len(args) == 0:
            if self.mag() < res:
                return 0.0
            if x >= 0 and y >= 0:
                try:
                    return math.atan(y / x)
                except ZeroDivisionError:
                    return math.pi / 2
            elif x < 0 and y >= 0:
                return math.pi - math.atan(y / abs(x))
            elif x >= 0 and y < 0:
                try:
                    return -math.atan(abs(y) / x)
                except ZeroDivisionError:
                    return -math.pi / 2
            else:
                return math.atan(abs(y) / abs(x)) - math.pi
        elif len(args) == 1:
            b = args[0]
            try:
                rv = math.acos(self.dot(b) / (self.mag() * b.mag()))
                return rv
            except ZeroDivisionError:
                return 0.0

    def rot(self, a):
        dx = self.dx
        dy = self.dy
        dx_t = self.dx * math.cos(a) - self.dy * math.sin(a)
        dy_t = self.dx * math.sin(a) + self.dy * math.cos(a)
        return vec(dx_t, dy_t)

    def norm(self):
        return self / self.mag()

    def floor(self):
        # return to normal vector in the plane normal to the self for
        # Complete coordinate system
        a = self.dx
        b = self.dy
        c = self.dz
        if self.mag() < res:
            raise geoError('ZeroEntity')
        if (abs(a) < res and abs(b) < res):
            return (vec(0, 0, 1), vec(0, 1, 0), vec(1, 0, 0))
        elif (abs(a) < res and abs(c) < res):
            return (vec(0, 1, 0), vec(1, 0, 0), vec(0, 0, 1))
        elif (abs(c) < res and abs(b) < res):
            return (vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1))
        else:
            ex = self
            if abs(c) > res:
                q = 0
                w = .5
                e = -w * b / c
                ey = vec(q, w, e)
            else:
                e = 0
                q = .5
                w = -q * a / b
                ey = vec(q, w, e)
            ez = ex.cross(ey)
        return (ex / ex.mag(), ey / ey.mag(), ez / ez.mag())


class circle:
    def __init__(self, p, r):
        self.c = p
        self.r = float(r)

    def __str__(self):
        return 'Circle[' + self.c.__str__() + ',' + str(self.r) + ']'

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def touch(self, o, fg=0):
        d = self.c.dist(o.c)
        if fg == 0:
            met = self.r + o.r
        else:
            met = 180.0 * (self.r + o.r) / (E.R * math.pi)
        if d < met + res:
            return True
        else:
            return False

    def side(self, p):
        d = p.dist(self.c)
        r = self.r
        if d < r:
            return 1
        elif d < r + res:
            return 0
        else:
            return -1

    def to_poly(self, step):
        alpha = 0
        pll = []
        while alpha < (2 * math.pi):
            x = self.c + point(self.r * math.cos(alpha), self.r * math.sin(alpha))
            pll.append(x)
            alpha = alpha + step
        pll[-1] = pll[0]
        return pll


class line:
    def __init__(self, x, y):
        if isinstance(x, point) and isinstance(y, point):
            if x == y:
                raise geoError('Degenerate')
            try:
                temp = (x.y - y.y) / (x.x - y.x)
                if temp == 0:
                    self.kind = 'Horizental'
                    self.b = x.y
                else:
                    self.kind = 'Normal'
                    self.b = x.y - temp * x.x
            except ZeroDivisionError:
                temp = float('inf')
                self.kind = 'Vertical'
                self.b = x.x
            self.m = temp
        elif (not isinstance(x, point)) and isinstance(y, point):
            self.m = float(x)
            if x == 0:
                self.kind = 'Horizental'
                self.b = y.y
            elif x == float('inf'):
                self.kind = 'Vertical'
                self.b = y.x
            else:
                self.kind = 'Normal'
                self.b = y.y - x * y.x
        else:
            self.m = float(x)
            if x == 0:
                self.kind = 'Horizental'
            elif x == float('inf'):
                self.kind = 'Vertical'
            else:
                self.kind = 'Normal'
            self.b = float(y)

    def __str__(self):
        if self.kind == 'Vertical':
            return 'Line(x=' + str(self.b) + ')'
        elif self.kind == 'Horizental':
            return 'Line(y=' + str(self.b) + ')'
        else:
            return 'Line(' + str(self.m) + 'x+' + str(self.b) + ')'

    def intersect(self, L):
        if isinstance(L, line):
            if self.m == L.m:
                raise geoError('parellel')
            else:
                if self.kind == 'Normal' and L.kind == 'Normal':
                    x = (self.b - L.b) / (L.m - self.m)
                    y = L.m * x + L.b
                    return point(x, y)
                elif self.kind == 'Horizental' and L.kind == 'Normal':
                    x = (self.b - L.b) / L.m
                    return point(x, self.b)
                elif self.kind == 'Vertical' and L.kind == 'Normal':
                    y = L.m * self.b + L.b
                    return point(self.b, y)
                elif L.kind == 'Horizental' and self.kind == 'Normal':
                    x = (L.b - self.b) / self.m
                    return point(x, L.b)
                elif L.kind == 'Vertical' and self.kind == 'Normal':
                    y = self.m * L.b + self.b
                    return point(L.b, y)
                elif self.kind == 'Horizental' and L.kind == 'Vertical':
                    return point(L.b, self.b)
                elif self.kind == 'Vertical' and L.kind == 'Horizental':
                    return point(self.b, L.b)
                else:
                    raise geoError('Unknown')
        else:
            p = L.c
            x0 = p.x
            y0 = p.y
            r = L.r
            m = self.m
            b = self.b
            c0 = b - y0
            if self.kind == 'Vertical':
                x = self.b
                sol = solve2p(1, -2 * y0, (x - x0) ** 2 - r ** 2 + y0 ** 2)
                w = []
                for yy in sol:
                    w.append(point(x, yy))
                return w
            else:
                sol = solve2p(1 + m ** 2, -2 * x0 + 2 * m * c0, x0 ** 2 + c0 ** 2 - r ** 2)
                w = []
                for x in sol:
                    w.append(point(x, m * x + b))
                return w

    def side(self, p):
        if slef.kind == 'Vertical':
            if abs(self.b - p.y) < res:
                met = 0.0
            elif self.b > p.y:
                met = 1
            else:
                met = -1
        elif self.kind == 'Horizental':
            if abs(self.b - p.x) < res:
                met = 0.0
            elif self.b > p.x:
                met = -1
            else:
                met = 1
        else:
            met = L.m * p.x + L.b - p.y
        if abs(met) < res:
            return 0
        elif met < 0.0:
            return -1
        else:
            return 1

    def side(self, p):
        if self.kind == 'Vertical':
            if p.x == self.b:
                return 0
            elif p.x > self.b:
                return -1
            else:
                return 1
        else:
            met = p.y - (self.m * p.x + self.b)
        if met < 0:
            return -1
        elif met > 0:
            return 1
        else:
            return 0


class Triangle:
    def __init__(self, L):
        x = L[0]
        y = L[1]
        z = L[-1]
        L = line(x, y)
        if L.side(z) == 0:
            raise geoError('Degenerate')
        self.a = x
        self.b = y
        self.c = z

    def __str__(self):
        return 'Triangle{' + self.a.__str__().strip('p') + ',' + self.b.__str__().strip(
            'p') + ',' + self.c.__str__().strip('p') + '}'

    def centroid(self):
        V01 = 0.5 * (self.a + self.b)
        V12 = 0.5 * (self.c + self.b)
        L1 = line(V01, self.c)
        L2 = line(V12, self.a)
        return L1.intersect(L2)

    def A(self):
        v1 = vec(self.a, self.b)
        v2 = vec(self.a, self.c)
        return 0.5 * abs(v1.cross(v2))


class Sphere:
    def __init__(self, c, R):
        self.c = c
        self.R = R

    def dist(self, p):
        return abs(p.dist(self.c) - self.R)

    def s2c(self, lon, lat):
        x = self.R * math.cos(lat) * math.cos(lon)
        y = self.R * math.cos(lat) * math.sin(lon)
        z = self.R * math.sin(lat)
        return self.c + point(x, y, z)

    def gcd(self, lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians
        lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
        c = 2 * math.asin(math.sqrt(a))

        dis = E.R * c
        return dis

    def gcd2l(self, d):
        # great circle distance to straight line distance
        theta = d / self.R
        return (2 * (self.R ** 2) * (1 - math.cos(theta))) ** 0.5

    def side(self, p):
        d = pdist(self.c, p)
        if d > R - res and d < R + res:
            return 0
        if d < R:
            return -1
        else:
            return 1

    def small_circle(self, S):
        R = self.R
        r = S.R
        d = self.c.dist(S.c)
        if d > (R + r):
            raise geoError('Disjoint')
        ex1 = vec(S.c - self.c)
        ex, ey, ez = ex1.floor()
        x = (d ** 2 - r ** 2 + R ** 2) / (2 * d)
        temp = (4 * (d ** 2) * (R ** 2)) - ((d ** 2) - (r ** 2) + (R ** 2)) ** 2
        a = math.sqrt(temp) / (2 * d)
        # sampling the small circle and transforming the coordinates
        n = 2 * int(math.pi / s_res)
        theta = num.linspace(0, 2 * math.pi, n)
        sc = []
        for alpha in theta:
            y = a * math.cos(alpha)
            z = a * math.sin(alpha)
            xyz = point(x, y, z)
            pf = xyz.x * ex + xyz.y * ey + xyz.z * ez + self.c
            sc.append(pf)
        return sc

    def map(self, p, inv=False):
        # cartesian to lat/lon
        # if inv is true lat/lon to cartesian
        R = self.R
        if not inv:
            ed = R * (vec(p - self.c).norm())
            ed = point(ed.dx, ed.dy, ed.dz)
            lon = math.atan2(ed.y, ed.x)
            lat1 = math.acos(abs(ed.z) / R)
            if ed.z > 0:
                lat = math.pi / 2 - lat1
            else:
                lat = -(math.pi / 2 - lat1)
            return point(lon, lat) * 180 / math.pi
        if inv:
            p = p * math.pi / 180
            z = R * math.sin(p.y)
            y = R * math.cos(p.y) * math.sin(p.x)
            x = R * math.cos(p.y) * math.cos(p.x)
            return point(x, y, z)


class Plane:
    def __init__(self, v, p):
        self.v = v
        self.p = p


class Polygon:
    def __init__(self, L):
        # L is the list of points on the polygon
        # L[-1]==L[0] Closed polygon
        if L[-1] == L[0]:
            self.L = L
            self.n = len(L)
        else:
            raise geoError('open')

    def area(self):
        s = 0
        for i in range(self.n - 1):
            p1 = self.L[i]
            p2 = self.L[i + 1]
            s = s + p1.x * p2.y - p2.x * p1.y
        return 0.5 * s

    def centroid(self):
        sx = 0
        sy = 0
        for i in range(self.n - 1):
            p1 = self.L[i]
            p2 = self.L[i + 1]
            sx = sx + (p1.x * p2.y - p2.x * p1.y) * (p1.x + p2.x)
            sy = sy + (p1.x * p2.y - p2.x * p1.y) * (p1.y + p2.y)
        p = point(sx, sy)
        return p / (6.0 * self.area())


class Ray:
    def __init__(self, p, a):
        self.c = p
        self.a = angleMap(a)

    def side(self, p):
        ba = p.angle(self.c)
        da = abs(ba - a) / (2 * math.pi)
        if da % (2 * math.pi) < res:
            return 1
        else:
            return 1

    def to_line(self):
        if (abs(self.a - math.pi / 2) < res) or (abs(self.a + math.pi / 2) < res):
            m = float('inf')
            b = self.c.x
            return line(m, b)
        m = math.tan(self.a)
        return line(m, self.c)

    def intersect(self, c):
        if isinstance(c, circle):
            L = self.to_line()
            p = L.intersect(c)
            if len(p) == 0:
                return
            else:
                a1 = self.c.angle(p[0])
                a2 = self.c.angle(p[1])
                if abs(a1 - a2) < res:
                    print('hey')
                    if p[0].dist(self.c) > p[1].dist(self.c):
                        return [p[1], p[0]]
                    else:
                        return [p[0], p[1]]
                elif abs(a1 - self.a) > abs(a2 - self.a):
                    return [p[1]]
                else:
                    return [p[0]]


class ndisc:
    def __init__(self, cA):
        self.cA = cA
        self.n = len(cA)

    def is_disjoint(self):
        l = len(self.cA)
        for i in range(l):
            for j in range(i + 1, l):
                if not self.cA[j].touch(self.cA[i]):
                    return True
        return False

    def sort(self):
        r = [xx.r for xx in self.cA]
        r, self.cA = order(r, self.cA)

    # ~ def remove_ins(self):
    # ~ self.sort()
    # ~ if cA[
    def side(self, p):
        counter = [0, 0, 0]  # [#inside,#on,#outside]
        for cir in self.cA:
            fg = cir.side(p)
            if fg == 1:
                counter[0] = counter[0] + 1
            elif fg == 0:
                counter[1] = counter[1] + 1
            else:
                counter[2] = counter[2] + 1
        return counter

    def get_x0(self, weighted=True):
        l = self.n
        r = [w.r for w in self.cA]
        c = [w.c for w in self.cA]
        S = sum(r)
        W = [(S - w) / ((l - 1) * S) for w in r]
        p0 = point(0, 0)  # Initialized point
        for i in range(l):
            p0 = p0 + W[i] * c[i]
        return p0

    def intersect(self, obj):
        l = self.n
        d = [0.0] * l
        pll = [point(0, 0)] * l
        if isinstance(obj, Ray):
            for i, cir in enumerate(self.cA):
                pll[i] = obj.intersect(cir)[0]
                d[i] = pll[i].dist(obj.c)
        usd, pn = order(d, pll)
        return pn[0]

    def poly(self, step):
        # step is the radian resolution
        pll = []
        cA = self.cA
        if len(cA) == 0:
            raise geoError('InputError')
        elif len(cA) == 1:
            return cA[0].to_poly(step)

        if self.is_disjoint():
            raise geoError('Disjoint')
        self.sort()
        cA = self.cA
        st = False
        fin = False
        cdi = 0
        alpha = 0
        x0 = self.cA[0].c + point(self.cA[0].r, 0)
        xc = x0
        # Finding the first point on Polygon
        while (alpha < (2 * math.pi)):
            if self.side(xc)[-1] == 0:
                pll.append(xc)
                st = True
                break
            else:
                xx = num.cos(step) * cA[cdi].r
                yy = num.sin(step) * cA[cdi].r
                xc_t = point(xx, yy)
                xc = xc_t.transform(-cA[cdi].c, -alpha)
                alpha = alpha + step
        if not st:
            dA = uni_obj(cA)
            if len(dA) == 1:
                return dA[0].to_poly()
            elif len(dA) == 2:
                p_t = circle_intersect(dA[0], dA[1])
                return [p_t]
            else:
                p_t = circles_intersect(dA)
                return [p_t]

        # pivoting
        pv0 = self.get_x0()
        pv = point(pv0.x, pv0.y)
        v0 = vec(xc, pv)
        vc = vec(xc, pv)
        spin = 0.0
        pv_found = False
        while spin < (2 * math.pi):
            if self.side(pv)[0] == self.n:
                pv_found = True
                break
            else:
                vc = vc / 2
                if vc.mag() < 10 * res:
                    spin = spin + 10 * step
                    vc = v0.rot(spin)
                pv = xc + vc
        # pivoting finished
        if not pv_found:
            return [xc]
        alpha = num.linspace(0, 2 * math.pi, int(2 * math.pi / step))
        alpha[-1] = math.pi * 2
        for a in alpha:
            ray = Ray(pv, a)
            try:
                pc = self.intersect(ray)
                pll.append(pc)
            except geoError:
                print('Unknown Error in ray-ndisc intersection')
                raise geoError('Unknown')
        pll[-1] = pll[0]
        return pll


E = Sphere(point(0, 0, 0), 6378100)
# ~

# ~ if __name__=='__main__':
# ~ vc=point(5439680.78365,-959162.488598,3189050.0)
# ~ ws=point(10,30)
# ~ wc=E.map(ws,inv=True)
# ~ st=wc.dist(vc)
# ~ vs=E.map(vc)
# ~ print( st)
# ~ yy=E.gcd(vs.x,vs.y,ws.x,ws.y)
# ~ print( E.gcd2l(yy))
