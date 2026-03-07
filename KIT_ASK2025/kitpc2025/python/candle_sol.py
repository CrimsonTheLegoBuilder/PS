import math
import sys
import math
from functools import cmp_to_key
INPUT = sys.stdin.readline

R: int = 1 << 0
G: int = 1 << 1
B: int = 1 << 2
Y: int = R | G
C: int = G | B
M: int = R | B
W: int = R | G | B

TOL: float = 1e-9
PI: float = math.pi


def sign(x: float) -> int:
    if x < -TOL:
        return -1
    if x > TOL:
        return 1
    return 0


def eq(x: float, y: float) -> bool:
    return sign(x - y) == 0


def fit(x: float, lo: float, hi: float) -> float:
    return max(lo, min(x, hi))


def cross(d1: tuple, d2: tuple) -> float:
    return d1[0] * d2[1] - d1[1] * d2[0]


def dot(d1: tuple, d2: tuple) -> float:
    return d1[0] * d2[0] + d1[1] * d2[1]


def ccw(d1: tuple, d2: tuple, d3: tuple) -> int:
    return sign(cross3(d1, d2, d3))


def cross3(d1: tuple, d2: tuple, d3: tuple) -> float:
    return (d2[0] - d1[0]) * (d3[1] - d2[1]) - (d2[1] - d1[1]) * (d3[0] - d2[0])


def cross4(d1: tuple, d2: tuple, d3: tuple, d4: tuple) -> float:
    return (d2[0] - d1[0]) * (d4[1] - d3[1]) - (d2[1] - d1[1]) * (d4[0] - d3[0])


def dot3(d1: tuple, d2: tuple, d3: tuple) -> float:
    return (d2[0] - d1[0]) * (d3[0] - d2[0]) + (d2[1] - d1[1]) * (d3[1] - d2[1])


def dot4(d1: tuple, d2: tuple, d3: tuple, d4: tuple) -> float:
    return (d2[0] - d1[0]) * (d4[0] - d3[0]) + (d2[1] - d1[1]) * (d4[1] - d3[1])


def mul(d: tuple, n: float) -> tuple:
    return d[0] * n, d[1] * n


def add(d1: tuple, d2: tuple) -> tuple:
    return d1[0] + d2[0], d1[1] + d2[1]


def sub(d1: tuple, d2: tuple) -> tuple:
    return d1[0] - d2[0], d1[1] - d2[1]


def euc(d1: tuple) -> float:
    return d1[0] ** 2 + d1[1] ** 2


def mag(d1: tuple) -> float:
    return euc(d1) ** 0.5


def unit(d1: tuple) -> tuple:
    l = mag(d1)
    return d1[0] / l, d1[1] / l


def on_seg_strong(d1: tuple, d2: tuple, t: tuple) -> bool:
    return cross3(d1, d2, t) == 0 and dot3(d1, t, d2) >= 0


def collinear(d1: tuple, d2: tuple, d3: tuple, d4: tuple) -> bool:
    return ccw(d1, d2, d3) == 0 and ccw(d1, d2, d4) == 0


def inner_check(p: list[tuple], t: tuple, d: tuple = (0, 0)) -> int:
    for i0 in range(len(p)):
        i1: int = (i0 + 1) % len(p)
        p0 = p[i0]
        p1 = p[i1]
        if cross3(p0, p1, t) < 0:
            return 0
        if on_seg_strong(p0, p1, t) and not eq(d[0], d[1]):
            if sign(cross(sub(p1, p0), d)) > 0:
                return 1
            else:
                return 0
        if on_seg_strong(p0, p1, t):
            return 2
    return 1


def intersection(s1: list, s2: list, f: int = 0) -> float:
    p1, p2 = s1
    q1, q2 = s2
    det = cross(sub(q2, q1), sub(p2, p1))
    if sign(det) == 0:
        return -1
    a1: float = cross(sub(q2, q1), sub(q1, p1)) / det
    a2: float = cross(sub(p2, p1), sub(p1, q1)) / det
    if f == 1:
        return a1
    if 0 < a1 < 1 and -TOL < a2 < 1 + TOL:
        return a1
    return -1


def intersect(p1: tuple, p2: tuple, q1: tuple, q2: tuple) -> bool:
    f1: bool = ccw(p1, p2, q1) * ccw(p1, p2, q2) > 0
    f2: bool = ccw(q1, q2, p1) * ccw(q1, q2, p2) > 0
    f3: bool = on_seg_strong(p1, p2, q1) or on_seg_strong(p1, p2, q2) \
    or on_seg_strong(q1, q2, p1) or on_seg_strong(q1, q2, p2)
    return (f1 and f2) or f3


def inside(d1: tuple, d2: tuple, d3: tuple, t: tuple, f: int = 0) -> bool:
    if ccw(d1, d2, d3) < 0:
        return ccw(d1, d2, t) >= f or ccw(d2, d3, t) >= f
    return  ccw(d1, d2, t) >= f and ccw(d2, d3, t) >= f


def pos(se: list, r: float = .5) -> tuple:
    v = sub(se[1], se[0])
    l = mag(v)
    return se[0] + mul(v, r / l)


def get_pos(l0: tuple, s1: list, s2: list) -> tuple:
    p1, p2 = s1
    q1, q2 = s2
    if not inside(p2, l0, p1, q1, 1) and not inside(p2, l0, 1, q2, 1):
        if intersect(l0, p1, q1, q2) and intersect(l0, p2, q1, q2):
            return 0, 1
        else:
            return 0, 0
    T: list[tuple] = [p1, p2, l0]
    in1: int = inner_check(T, q1)
    in2: int = inner_check(T, q2)
    if not in1 and not in2:
        return 0, 0
    r1: float = 0
    r2: float = 0
    if in1 and in2:
        r1 = intersection(s1, [l0, q1], 1)
        r2 = intersection(s1, [l0, q2], 1)
    elif in1:
        r1 = intersection(s1, [l0, q1], 1)
    elif in2:
        r2 = intersection(s1, [l0, q2], 1)
    else:
        r1, r2 = 0, 0
    if r2 < r1:
        r2, r1 = r1, r2
    return r1, r2


def inner_check_concave(p: list[tuple], t: tuple, s: tuple, e: tuple) -> int:
    c: int = 0
    for i0 in range(len(p)):
        i1: int = (i0 + 1) % len(p)
        p0: tuple = p[i0]
        p1: tuple = p[i1]
        if on_seg_strong(p0, p1, t):
            assert collinear(p0, p1, s, e)
            if dot4(p0, p1, s, e) > 0:
                return 1
            else:
                return 0
        if eq(p0[1], p1[1]):
            continue
        if p1[1] < p0[1]:
            p0, p1 = p1, p0
        if p1[1] - TOL < t[1] or p0[1] > t[1]:
            continue
        c += ccw(p0, p1, t) > 0
    return c & 1


if __name__ == "__main__":
    n: int = int(INPUT())
    h0 = [tuple(*map(int, INPUT().strip().split())) for _ in range(n)]
    H = [h0]
    b: int = int(INPUT())
    for i in range(b):
        m: int = int(INPUT())
        h = [tuple(*map(int, INPUT().strip().split())) for _ in range(m)]
        H.append(h)
    q: int = int(INPUT())
    for tc in range(1, q + 1):
        rx, ry, bx, by, gx, gy = map(int, INPUT().strip().split())
        cdl = [(0, 0) for _ in range(8)]
        cdl[R] = (rx, ry)
        cdl[G] = (gx, gy)
        cdl[B] = (bx, by)
        c: list[int] = [0 for _ in range(150)]
        I: list[int] = [-1 for _ in range(8)]
        A: list[float] = [0. for _ in range(8)]
        Z: list[list] = [[] for _ in range(8)]
        S: list[list] = [[] for _ in range(8)]
        for i in [R, G, B]:
            f0: int = inner_check(H[0], cdl[i])
            fk: int = 0
            if f0 != 1:
                continue
            if f0:
                c[0] |= i
                I[i] = 0
            for k in range(1, b + 1):
                fk = inner_check(H[k], cdl[i])
                if fk:
                    c[0] -= i
                    c[k] |= i
                    I[i] = k
                    break
            if fk:
                continue
            for k in range(1, b + 1):
                m = len(H[k])
                pl = H[k][0]
                pr = H[k][0]
                for j in range(m):
                    if ccw(cdl[i], pl, H[k][j]) > 0:
                        pl = H[k][j]
                    if ccw(cdl[i], pl, H[k][j]) < 0:
                        pr = H[k][j]
                S[i].append([pl, pr])
            sz = len(H[0])
            for j in range(sz):
                u, v = H[0][j], H[0][(j + 1) % sz]
                assert ccw(cdl[i], u, v) > 0
                w = [u, v]
                vp: list[tuple] = [(0, 0)]
                for k in range(1, b + 1):
                    se: tuple = get_pos(cdl[i], w, S[i][k - 1])
                    if not eq(se[0], se[1]):
                        vp.append(se)
                vp.append((1, 1))
                vp.sort(key=lambda o: (o[0], o[1]))
                hi: float = 0
                for p in vp:
                    if hi < p[0]:
                        s: tuple = pos(w, hi)
                        e: tuple = pos(w, p[0])
                        Z[i].append([s, e])
                        hi = p[1]
                    else:
                        hi = max(hi, p[1])
            for k in range(0, b + 1):
                sz = len(H[k])
                for j in range(sz):
                    u, v = H[k][j], H[k][(j + 1) % sz]
                    if ccw(cdl[i], u, v) >= 0:
                        continue
                    w = [u, v]



