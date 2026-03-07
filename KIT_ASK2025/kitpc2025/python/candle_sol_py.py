# Python 3
import sys
import math
from functools import cmp_to_key

# --- Types & constants (C++ 대응) ---
ll = int
ld = float

INF = 1e17
TOL = 1e-9
PI = math.acos(-1.0)
LEN = 150

STRONG = 0
WEAK = 1

BLACK   = 0
RED     = (1 << 0)
GREEN   = (1 << 1)
BLUE    = (1 << 2)
YELLOW  = (RED | GREEN)
MAGENTA = (RED | BLUE)
CYAN    = (GREEN | BLUE)
WHITE   = (RED | GREEN | BLUE)

def sign_l(x: ll) -> int:
    return -1 if x < 0 else (1 if x > 0 else 0)

def sign(x: ld) -> int:
    if x < -TOL: return -1
    if x >  TOL: return  1
    return 0

def zero(x: ld) -> bool:
    return sign(x) == 0

def eq(x: ld, y: ld) -> bool:
    return zero(x - y)

def sq(x):
    return x * x

def norm(th: ld) -> ld:
    while th < 0: th += 2 * PI
    while sign(th - 2 * PI) >= 0: th -= 2 * PI
    return th

def fit(x: ld, lo: ld = 0.0, hi: ld = 1.0) -> ld:
    return min(hi, max(lo, x))

# --- Geometry primitives ---

class Pos:
    __slots__ = ("x", "y")
    def __init__(self, x: ld = 0.0, y: ld = 0.0):
        self.x = float(x); self.y = float(y)
    def __eq__(self, p) -> bool:
        return zero(self.x - p.x) and zero(self.y - p.y)
    def __lt__(self, p) -> bool:
        # zero(x - p.x) ? y < p.y : x < p.x
        if zero(self.x - p.x):
            return self.y < p.y
        return self.x < p.x
    def __add__(self, p):
        return Pos(self.x + p.x, self.y + p.y)
    def __sub__(self, p):
        return Pos(self.x - p.x, self.y - p.y)
    def __mul__(self, a):
        # scalar scale or dot product
        if isinstance(a, Pos):  # dot
            return self.x * a.x + self.y * a.y
        else:                   # scalar
            return Pos(self.x * a, self.y * a)
    def __rmul__(self, a):
        return self.__mul__(a)
    def __truediv__(self, p):
        # cross (as in C++ operator/ with Pos)
        return self.x * p.y - self.y * p.x
    def __invert__(self):
        # ~ : rotate 90deg (unused here)
        return Pos(-self.y, self.x)
    def __neg__(self):
        # ! in C++ gave {-x, -y}; here just use -p
        return Pos(-self.x, -self.y)
    def Euc(self) -> ld:
        return self.x * self.x + self.y * self.y
    def mag(self) -> ld:
        return math.sqrt(self.Euc())
    def unit(self):
        m = self.mag()
        return self if zero(m) else (self * (1.0 / m))
    def __repr__(self):
        return f"Pos({self.x},{self.y})"

O = Pos(0.0, 0.0)
Polygon = list  # list[Pos]

def cross(d1: Pos, d2: Pos, d3: Pos) -> ld:
    return (d2 - d1) / (d3 - d2)

def cross4(d1: Pos, d2: Pos, d3: Pos, d4: Pos) -> ld:
    return (d2 - d1) / (d4 - d3)

def ccw(d1: Pos, d2: Pos, d3: Pos) -> int:
    return sign(cross(d1, d2, d3))

def ccw4(d1: Pos, d2: Pos, d3: Pos, d4: Pos) -> int:
    return sign(cross4(d1, d2, d3, d4))

def dot3(d1: Pos, d2: Pos, d3: Pos) -> ld:
    return (d2 - d1) * (d3 - d2)

def dot4(d1: Pos, d2: Pos, d3: Pos, d4: Pos) -> ld:
    return (d2 - d1) * (d4 - d3)

def on_seg_strong(d1: Pos, d2: Pos, d3: Pos) -> bool:
    return (ccw(d1, d2, d3) == 0) and (sign(dot3(d1, d3, d2)) >= 0)

def on_seg_weak(d1: Pos, d2: Pos, d3: Pos) -> bool:
    return (ccw(d1, d2, d3) == 0) and (sign(dot3(d1, d3, d2)) > 0)

def collinear(d1: Pos, d2: Pos, d3: Pos, d4: Pos) -> bool:
    return (ccw(d1, d2, d3) == 0) and (ccw(d1, d2, d4) == 0)

def intersection_point(p1: Pos, p2: Pos, q1: Pos, q2: Pos) -> Pos:
    a1 = cross(q1, q2, p1)
    a2 = -cross(q1, q2, p2)
    return (p1 * a2 + p2 * a1) * (1.0 / (a1 + a2))

def projection(d1: Pos, d2: Pos, d3: Pos) -> ld:
    return dot4(d1, d2, d1, d3) / (d2 - d1).Euc()

def inside(p0: Pos, p1: Pos, p2: Pos, q: Pos, f: int = STRONG) -> bool:
    # triangle half-plane test as in original
    if ccw(p0, p1, p2) < 0:
        return (ccw(p0, p1, q) >= f) or (ccw(p1, p2, q) >= f)
    return (ccw(p0, p1, q) >= f) and (ccw(p1, p2, q) >= f)

def intersect_seg(s1: Pos, s2: Pos, d1: Pos, d2: Pos) -> bool:
    f1 = ccw(s1, s2, d1) * ccw(s2, s1, d2) > 0
    f2 = ccw(d1, d2, s1) * ccw(d2, d1, s2) > 0
    f3 = on_seg_strong(s1, s2, d1) or on_seg_strong(s1, s2, d2) or on_seg_strong(d1, d2, s1) or on_seg_strong(d1, d2, s2)
    return (f1 and f2) or f3

def area_poly(H: Polygon) -> ld:
    sz = len(H)
    if sz < 3: return 0.0
    a = 0.0
    for i in range(sz):
        a += H[i] / H[(i + 1) % sz]
    return a * 0.5

class Seg:
    __slots__ = ("s", "e")
    def __init__(self, s: Pos = None, e: Pos = None):
        self.s = s if s is not None else Pos()
        self.e = e if e is not None else Pos()
    def dir(self) -> Pos:
        return (self.s - self.e).unit()
    def __lt__(self, other) -> bool:
        return self.s == other.s and self.e < other.e or self.s < other.s
    def __eq__(self, other) -> bool:
        return self.s == other.s and self.e == other.e
    def p(self, rt: ld = 0.5) -> Pos:
        return self.s + (self.e - self.s) * rt
    def green(self, lo: ld = 0.0, hi: ld = 1.0) -> ld:
        d = hi - lo
        ratio = (lo + hi) * 0.5
        m = self.p(ratio)
        return m.y * d * (self.s.x - self.e.x)

def dot_seg(p: Seg, q: Seg) -> ld:
    return dot4(p.s, p.e, q.s, q.e)

def intersect_seg2(u: Seg, v: Seg) -> bool:
    return intersect_seg(u.s, u.e, v.s, v.e)

def intersection_param(s1: Seg, s2: Seg, f: bool = STRONG) -> ld:
    p1, p2, q1, q2 = s1.s, s1.e, s2.s, s2.e
    det = (q2 - q1) / (p2 - p1)
    if zero(det): return -1.0
    a1 = ((q2 - q1) / (q1 - p1)) / det
    a2 = ((p2 - p1) / (p1 - q1)) / -det
    if f == WEAK:
        return fit(a1, 0.0, 1.0)
    if (0.0 < a1 < 1.0) and (-TOL < a2 < 1.0 + TOL):
        return a1
    return -1.0

def inner_check(H: Polygon, q: Pos, d: Pos = Pos(0.0, 0.0)) -> int:
    # convex polygon
    sz = len(H)
    for i in range(sz):
        p1, p2 = H[i], H[(i + 1) % sz]
        if ccw(p1, p2, q) < 0:
            return 0
        if on_seg_strong(p1, p2, q) and not eq(d.x, d.y):
            return 1 if sign((p2 - p1) / d) > 0 else 0
        if on_seg_strong(p1, p2, q):
            return 2
    return 1

def inner_check_concave(H: Polygon, p: Pos, s: Pos, e: Pos) -> bool:
    # ray casting, with tie-break on collinear
    cnt = 0
    sz = len(H)
    for i in range(sz):
        cur = H[i]
        nxt = H[(i + 1) % sz]
        if on_seg_strong(cur, nxt, p):
            assert collinear(cur, nxt, s, e)
            return dot4(cur, nxt, s, e) > 0
        if zero(cur.y - nxt.y):
            continue
        c, n = cur, nxt
        if n.y < c.y: c, n = n, c
        if n.y - TOL < p.y or c.y > p.y:
            continue
        cnt += (ccw(c, n, p) > 0)
    return (cnt & 1) == 1

def get_pos(l: Pos, p: Seg, q: Seg) -> Pos:
    p1, p2, q1, q2 = p.s, p.e, q.s, q.e
    # if both endpoints of q are outside triangle (p2,l,p1)
    if (not inside(p2, l, p1, q1, WEAK)) and (not inside(p2, l, p1, q2, WEAK)):
        if intersect_seg(l, p1, q1, q2) and intersect_seg(l, p2, q1, q2):
            return Pos(0.0, 1.0)
        else:
            return Pos(0.0, 0.0)
    tri = [p1, p2, l]
    in1 = inner_check(tri, q1)
    in2 = inner_check(tri, q2)
    if (not in1) and (not in2):
        return Pos(0.0, 0.0)
    r1 = 0.0; r2 = 1.0
    if in1 and in2:
        r1 = intersection_param(p, Seg(l, q1), WEAK)
        r2 = intersection_param(p, Seg(l, q2), WEAK)
    elif in1:
        r1 = intersection_param(p, Seg(l, q1), WEAK)
    elif in2:
        r2 = intersection_param(p, Seg(l, q2), WEAK)
    else:
        r1 = r2 = 0.0
    if r2 < r1:
        r1, r2 = r2, r1
    return Pos(r1, r2)

class Event:
    __slots__ = ("x", "f")
    def __init__(self, x: ld, f: int):
        self.x = float(x); self.f = int(f)
    def __lt__(self, o):
        return (sign(self.x - o.x) < 0) or (eq(self.x, o.x) and self.f < o.f)
    def __eq__(self, o):
        return eq(self.x, o.x) and self.f == o.f
    def __repr__(self):
        return f"Event({self.x},{self.f})"

# --- Globals mirrored from C++ ---
N = M = K = Q = 0
I = [0] * (1 << 3)
C = [0] * LEN
A = [0.0] * (1 << 3)
Lpos = [Pos() for _ in range(1 << 3)]

P = []            # P[0..K]: list[Polygon]
H = [[] for _ in range(1 << 3)]        # H[c]: Polygon
S = [[Seg() for _ in range(LEN)] for _ in range(1 << 3)]  # S[c][k]
Z = [[] for _ in range(1 << 3)]        # Z[c]: list[Seg]
VE = [[[] for _ in range(1 << 3)] for _ in range(1 << 3)] # VE[c1][c2] -> list[list[Event]]
VX = [[[] for _ in range(1 << 3)] for _ in range(1 << 3)] # VX[c1][c2] -> list[list[float]]

cen = Pos()

def cmpt(p: Seg, q: Seg) -> int:
    # assert ccw(cen, p.s, p.e) and ccw(cen, q.s, q.e) non-zero
    assert ccw(cen, p.s, p.e) != 0
    assert ccw(cen, q.s, q.e) != 0
    u = p.s - cen
    v = q.s - cen
    f1 = O < u
    f2 = O < v
    if f1 != f2:
        return -1 if f1 else 1
    # assert not zero(u / v)
    cr = u / v
    assert not zero(cr)
    return -1 if cr > 0 else 1

def unique_consecutive_pos(arr: list) -> list:
    if not arr: return arr
    out = [arr[0]]
    for p in arr[1:]:
        if not (p == out[-1]):
            out.append(p)
    return out

def unique_consecutive_seg(arr: list) -> list:
    if not arr: return arr
    out = [arr[0]]
    for s in arr[1:]:
        if not (s == out[-1]):
            out.append(s)
    return out

def cut(s: Pos, e: Pos, p0: Pos, p1: Pos, p2: Pos, V: list, tmp: list):
    s0 = Seg(s, e)
    if on_seg_strong(s, e, p1):
        x = projection(s, e, p1)
        tq1 = ccw(p0, p1, p2)
        assert tq1 != 0
        if ccw(s, e, p0) == 0:
            d = sign(dot4(s, e, p0, p1))
            if d > 0:
                V.append(x)
                if tq1 > 0: tmp.append(Event(x, -1))
                else: tmp.append(Event(x, 1)); tmp.append(Event(x, -1))
            else:
                if tq1 > 0:
                    return
                else:
                    V.append(x); tmp.append(Event(x, -1))
        elif ccw(s, e, p2) == 0:
            d = sign(dot4(s, e, p1, p2))
            if d > 0:
                V.append(x)
                if tq1 > 0: tmp.append(Event(x, 1))
                else: tmp.append(Event(x, 1)); tmp.append(Event(x, -1))
            else:
                if tq1 > 0:
                    return
                else:
                    V.append(x); tmp.append(Event(x, 1))
        else:
            tq0 = ccw(s, e, p0); tq2 = ccw(s, e, p2)
            assert tq0 != 0 and tq2 != 0
            if tq0 != tq2:
                tq = ccw4(s, e, p0, p2)
                assert tq != 0
                V.append(x)
                if tq < 0: tmp.append(Event(x, 1))
                else: tmp.append(Event(x, -1))
            elif tq1 < 0:
                x = projection(s, e, p1)
                V.append(x)
                tmp.append(Event(x, 1)); tmp.append(Event(x, -1))
    else:
        if ccw(s, e, p2) == 0:
            return
        sk = Seg(p1, p2)
        x = intersection_param(s0, sk)
        if x < 0: return
        tq = ccw4(s, e, p1, p2)
        assert tq != 0
        V.append(x)
        if tq < 0: tmp.append(Event(x, 1))
        else: tmp.append(Event(x, -1))

def clean(s: Pos, e: Pos, s0: Seg, c1: int, tmp: list, ve: list):
    if not tmp:
        if inner_check_concave(H[c1], s0.p(0.5), s, e):
            ve.append(Event(0.0, 1))
            ve.append(Event(1.0, -1))
        return
    tmp.sort()
    if tmp[-1].f == 1:
        ve.append(Event(1.0, -1))
    if tmp[0].f == -1:
        ve.append(Event(0.0, 1))
    ve.extend(tmp)

def query_one():
    global Lpos, A, I, C, H, Z, S, VE, VX, cen

    # 입력: 각 색 점 3개
    rx, ry = map(float, sys.stdin.readline().split())
    gx, gy = map(float, sys.stdin.readline().split())
    bx, by = map(float, sys.stdin.readline().split())
    Lpos[RED] = Pos(rx, ry)
    Lpos[GREEN] = Pos(gx, gy)
    Lpos[BLUE] = Pos(bx, by)

    A = [0.0] * (1 << 3)
    I = [-1] * (1 << 3)
    C = [0] * LEN
    H = [[] for _ in range(1 << 3)]
    Z = [[] for _ in range(1 << 3)]

    # 색별 처리
    for c in (RED, GREEN, BLUE):
        l = Lpos[c]
        f0 = inner_check(P[0], l)  # 0,1,2
        if f0 != 1:   # strictly inside만 진행
            continue
        if f0:
            C[0] |= c
            I[c] = 0
        fk = 0
        for k in range(1, K + 1):
            fk = inner_check(P[k], l)
            if fk:
                C[0] -= c
                C[k] |= c
                I[c] = k
                break
        if fk:
            continue

        # 각 파란 다각형에 대한 '지지' 세그먼트 (극점)
        for k in range(1, K + 1):
            T = P[k]
            M = len(T)
            pl = T[0]
            pr = T[0]
            for j in range(M):
                if ccw(l, pl, T[j]) > 0: pl = T[j]
                if ccw(l, pr, T[j]) < 0: pr = T[j]
            S[c][k] = Seg(pr, pl)

        # 경계 세그먼트 수집 Z[c]
        for k in range(0, K + 1):
            T = P[k]
            sz = len(T)
            for i in range(sz):
                u = T[i]; v = T[(i + 1) % sz]
                if k == 0:
                    assert ccw(l, u, v) > 0
                if (k != 0) and ccw(l, u, v) >= 0:
                    continue
                w = Seg(u, v) if k == 0 else Seg(v, u)
                VP = [Pos(0.0, 0.0)]
                for kk in range(1, K + 1):
                    if kk == k:
                        continue
                    se = get_pos(l, w, S[c][kk])
                    if not eq(se.x, se.y):
                        VP.append(se)
                VP.append(Pos(1.0, 1.0))
                VP.sort()
                hi = 0.0
                for p in VP:
                    if hi < p.x:
                        s = w.p(hi); e = w.p(p.x)
                        Z[c].append(Seg(s, e))
                        hi = p.y
                    else:
                        hi = max(hi, p.y)

        # 각 세그먼트 각도순 정렬 후 꼭짓점 나열 -> H[c]
        cen = l
        Z[c].sort(key=cmp_to_key(cmpt))
        for se in Z[c]:
            H[c].append(se.s)
            H[c].append(se.e)
        H[c] = unique_consecutive_pos(H[c])
        if H[c] and H[c][0] == H[c][-1]:
            H[c].pop()

    # 이벤트 테이블 초기화
    VE = [[[] for _ in range(1 << 3)] for _ in range(1 << 3)]
    VX = [[[] for _ in range(1 << 3)] for _ in range(1 << 3)]

    # 색 쌍(c1,c2) 이벤트 (양쪽 점이 모두 I==0일 때만)
    for i in range(3):
        c1 = (1 << ((i + 1) % 3))
        c2 = (1 << ((i + 2) % 3))
        if I[c1] == 0 and I[c2] == 0:
            for _ in range(2):
                VE[c1][c2] = []
                VX[c1][c2] = []
                sz1 = len(H[c1])
                sz2 = len(H[c2])
                VE[c1][c2] = [[] for _ in range(sz1)]
                VX[c1][c2] = [[] for _ in range(sz1)]
                for j in range(sz1):
                    tmp = []
                    j1 = (j + 1) % sz1
                    s = H[c1][j]; e = H[c1][j1]
                    s0 = Seg(s, e)
                    for k in range(sz2):
                        k0 = (k - 1 + sz2) % sz2
                        k2 = (k + 1) % sz2
                        p0 = H[c2][k0]; p1 = H[c2][k]; p2 = H[c2][k2]
                        cut(s, e, p0, p1, p2, VX[c1][c2][j], tmp)
                    ve = []
                    clean(s, e, s0, c2, tmp, ve)
                    ve.sort()
                    VX[c1][c2][j].sort()
                    VE[c1][c2][j] = ve
                c1, c2 = c2, c1  # swap once

    # R & G & B (세 점 모두 I==0)
    if I[RED] == 0 and I[GREEN] == 0 and I[BLUE] == 0:
        VS = []
        for i in range(3):
            c0 = (1 << i)
            c1 = (1 << ((i + 1) % 3))
            c2 = (1 << ((i + 2) % 3))
            sz = len(H[c0])
            for j in range(sz):
                # ve = merge(VE[c0][c1][j], VE[c0][c2][j])
                ve = sorted(VE[c0][c1][j] + VE[c0][c2][j])
                # surround with [0,1] sentinels
                ve = [Event(0.0, 1)] + ve + [Event(1.0, -1)]
                V = sorted(VX[c0][c1][j] + VX[c0][c2][j])
                V = [0.0] + V + [1.0]
                # unique V
                V2 = []
                for v in V:
                    if not V2 or not eq(V2[-1], v):
                        V2.append(v)
                V = V2
                j1 = (j + 1) % sz
                s = H[c0][j]; e = H[c0][j1]
                s0 = Seg(s, e)
                vp = []
                szr = len(ve); szx = len(V)
                cnt = 0; k = 0
                for x in range(szx - 1):
                    lo = V[x]; hi = V[x + 1]
                    while k < szr and eq(ve[k].x, lo):
                        cnt += ve[k].f; k += 1
                    if cnt > 2:
                        vp.append(Pos(lo, hi))
                for se in vp:
                    VS.append(Seg(s0.p(se.x), s0.p(se.y)))
        VS.sort()
        VS = unique_consecutive_seg(VS)
        for se in VS:
            A[WHITE] += se.green()

    # R&G, G&B, B&R
    for i in range(3):
        c1 = (1 << ((i + 1) % 3))
        c2 = (1 << ((i + 2) % 3))
        if I[c1] == 0 and I[c2] == 0:
            c = c1 | c2
            VS = []
            for ca in (c1, c2):
                cb = c2 if ca == c1 else c1
                sza = len(H[ca]); szb = len(H[cb])
                for j in range(sza):
                    ve = [Event(0.0, 1)]
                    V  = [0.0]
                    ve.extend(VE[ca][cb][j])
                    V.extend(VX[ca][cb][j])
                    ve.append(Event(1.0, -1))
                    V.append(1.0)
                    # unique V
                    V2 = []
                    for v in sorted(V):
                        if not V2 or not eq(V2[-1], v):
                            V2.append(v)
                    V = V2
                    j1 = (j + 1) % sza
                    s = H[ca][j]; e = H[ca][j1]
                    s0 = Seg(s, e)
                    vp = []
                    szr = len(ve); szx = len(V)
                    cnt = 0; k = 0
                    for x in range(szx - 1):
                        lo = V[x]; hi = V[x + 1]
                        while k < szr and eq(ve[k].x, lo):
                            cnt += ve[k].f; k += 1
                        if cnt > 1:
                            vp.append(Pos(lo, hi))
                    for se in vp:
                        VS.append(Seg(s0.p(se.x), s0.p(se.y)))
            VS.sort()
            VS = unique_consecutive_seg(VS)
            for se in VS:
                A[c] += se.green()

    # 단색: 경계 다각형 면적
    for c in (RED, GREEN, BLUE):
        if I[c] == 0:
            A[c] += area_poly(H[c])

    # 파란 다각형들 면적 누적
    for k in range(1, K + 1):
        A[C[k]] += area_poly(P[k])

    # 포함-배제
    if I[RED] == 0 and I[GREEN] == 0 and I[BLUE] == 0:
        for c in range(1, WHITE):
            A[c] -= A[WHITE]

    for i in range(3):
        c1 = (1 << ((i + 1) % 3))
        c2 = (1 << ((i + 2) % 3))
        if I[c1] == 0 and I[c2] == 0:
            c = c1 | c2
            A[c1] -= A[c]
            A[c2] -= A[c]

    A[BLACK] = area_poly(P[0])
    for c in range(1, (1 << 3)):
        A[BLACK] -= A[c]
    for c in range(0, (1 << 3)):
        if A[c] < 0: A[c] = 0.0

    # 출력 (소수 13자리)
    out = []
    for c in (RED, GREEN, BLUE, YELLOW, MAGENTA, CYAN, WHITE, BLACK):
        out.append(f"{A[c]:.13f}")
    sys.stdout.write("\n".join(out) + "\n")

def solve():
    global N, K, Q, P
    data = sys.stdin.read().strip().split()
    it = iter(data)
    try:
        N = int(next(it))
    except StopIteration:
        return
    # P[0]
    P = [ [] for _ in range(0) ]  # reset
    P0 = []
    for _ in range(N):
        x = float(next(it)); y = float(next(it))
        P0.append(Pos(x, y))
    P = [P0]
    K = int(next(it))
    for _k in range(1, K + 1):
        M = int(next(it))
        poly = []
        for _ in range(M):
            x = float(next(it)); y = float(next(it))
            poly.append(Pos(x, y))
        P.append(poly)
    Q = int(next(it))
    # 남은 토큰을 다시 스트림으로
    rest = list(it)
    # 재구성해서 라인별로 넘기자 (각 쿼리는 3줄 x 2실수)
    # 여기선 편의상 바로 sys.stdin 재사용 대신, 작은 래퍼로 읽음
    from io import StringIO
    buf = []
    for i in range(0, len(rest), 2):
        buf.append(f"{rest[i]} {rest[i+1]}\n")
    sys.stdin = StringIO("".join(buf))
    for _q in range(Q):
        query_one()

if __name__ == "__main__":
    import time
    sys.stdin = open("../tests/candle/in/90.in", "r", encoding="utf-8")
    t0_wall = time.perf_counter()
    t0_cpu  = time.process_time()
    solve()
    dt_wall = time.perf_counter() - t0_wall
    dt_cpu  = time.process_time() - t0_cpu
    print(f"[wall] {dt_wall*1000:.2f} ms  [cpu] {dt_cpu*1000:.2f} ms", file=sys.stderr)

