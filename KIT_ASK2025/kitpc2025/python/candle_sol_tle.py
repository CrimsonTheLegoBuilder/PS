# Python 3
import sys
import math
import time
from functools import cmp_to_key

# ---- consts ----
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

# ---- utils ----
def sign(x: float) -> int:
    if x < -TOL: return -1
    if x >  TOL: return  1
    return 0
def zero(x: float) -> bool: return sign(x) == 0
def eq(a: float, b: float) -> bool: return zero(a - b)
def fit(x: float, lo: float = 0.0, hi: float = 1.0) -> float: return min(hi, max(lo, x))

# ---- geometry ----
class Pos:
    __slots__ = ("x", "y")
    def __init__(self, x=0.0, y=0.0): self.x=float(x); self.y=float(y)
    def __eq__(self, o): return zero(self.x-o.x) and zero(self.y-o.y)
    def __lt__(self, o):
        if zero(self.x-o.x): return self.y < o.y
        return self.x < o.x
    def __add__(self, o): return Pos(self.x+o.x, self.y+o.y)
    def __sub__(self, o): return Pos(self.x-o.x, self.y-o.y)
    def __mul__(self, a):
        if isinstance(a, Pos): return self.x*a.x + self.y*a.y  # dot
        return Pos(self.x*a, self.y*a)
    def __rmul__(self, a): return self.__mul__(a)
    def __truediv__(self, o): return self.x*o.y - self.y*o.x    # cross
    def Euc(self): return self.x*self.x + self.y*self.y
    def pstr(self): return f"({self.x},{self.y})"  # debug

O = Pos(0.0, 0.0)
Polygon = list  # list[Pos]

def cross(d1: Pos, d2: Pos, d3: Pos) -> float:  return (d2-d1) / (d3-d2)
def cross4(d1: Pos, d2: Pos, d3: Pos, d4: Pos) -> float: return (d2-d1) / (d4-d3)
def ccw(d1: Pos, d2: Pos, d3: Pos) -> int: return sign(cross(d1,d2,d3))
def dot4(d1: Pos, d2: Pos, d3: Pos, d4: Pos) -> float: return (d2-d1) * (d4-d3)

def on_seg_strong(d1: Pos, d2: Pos, d3: Pos) -> bool:
    return (ccw(d1,d2,d3)==0) and (sign((d3-d1)*(d2-d3)) >= 0)

def collinear(d1: Pos, d2: Pos, d3: Pos, d4: Pos) -> bool:
    return (ccw(d1,d2,d3)==0) and (ccw(d1,d2,d4)==0)

def intersection_point(p1: Pos, p2: Pos, q1: Pos, q2: Pos) -> Pos:
    a1 = cross(q1,q2,p1); a2 = -cross(q1,q2,p2)
    return (p1*a2 + p2*a1) * (1.0/(a1+a2))

def projection(d1: Pos, d2: Pos, d3: Pos) -> float:
    return dot4(d1,d2,d1,d3) / (d2-d1).Euc()

def inside(p0: Pos, p1: Pos, p2: Pos, q: Pos, f: int = STRONG) -> bool:
    if ccw(p0,p1,p2) < 0: return (ccw(p0,p1,q) >= f) or (ccw(p1,p2,q) >= f)
    return (ccw(p0,p1,q) >= f) and (ccw(p1,p2,q) >= f)

def intersect_seg(s1: Pos, s2: Pos, d1: Pos, d2: Pos) -> bool:
    f1 = ccw(s1,s2,d1) * ccw(s2,s1,d2) > 0
    f2 = ccw(d1,d2,s1) * ccw(d2,d1,s2) > 0
    f3 = on_seg_strong(s1,s2,d1) or on_seg_strong(s1,s2,d2) or on_seg_strong(d1,d2,s1) or on_seg_strong(d1,d2,s2)
    return (f1 and f2) or f3

def area_poly(H: Polygon) -> float:
    sz = len(H)
    if sz < 3: return 0.0
    a = 0.0
    for i in range(sz): a += H[i] / H[(i+1)%sz]
    return a*0.5

def norm_poly(H: Polygon):
    a = area_poly(H)
    if a < 0: H.reverse()

def graham_scan(C: Polygon) -> Polygon:
    H = []
    if len(C) < 3:
        C.sort(); return C
    C[0], mn = min((p,i) for i,p in enumerate(C))[0], min(range(len(C)), key=lambda i:(C[i].x,C[i].y))
    C[0], C[mn] = C[mn], C[0]
    base = C[0]
    C[1:] = sorted(C[1:], key=lambda p: ( -ccw(base,p,Pos(base.x+1,base.y)), (base-p).Euc() ))
    # unique
    out=[C[0]]
    for p in C[1:]:
        if not (p==out[-1]): out.append(p)
    C = out
    for p in C:
        while len(H)>=2 and ccw(H[-2],H[-1],p) <= 0: H.pop()
        H.append(p)
    return H

def polygon_cut(ps: Polygon, b1: Pos, b2: Pos) -> Polygon:
    qs = []
    n = len(ps)
    for i in range(n):
        p1 = ps[i]; p2 = ps[(i+1)%n]
        d1 = ccw(b1,b2,p1); d2 = ccw(b1,b2,p2)
        if d1 >= 0: qs.append(p1)
        if d1*d2 < 0: qs.append(intersection_point(p1,p2,b1,b2))
    return qs

def suth_hodgman(C: Polygon, clip: Polygon) -> Polygon:
    ret = C[:]
    for i in range(len(clip)):
        b1 = clip[i]; b2 = clip[(i+1)%len(clip)]
        ret = polygon_cut(ret,b1,b2)
        if not ret: break
    return ret

class Seg:
    __slots__=("s","e","i")
    def __init__(self, s=None, e=None, i=-1): self.s=s or Pos(); self.e=e or Pos(); self.i=i
    def __lt__(self, o): return (self.s==o.s and self.e<o.e) or (self.s<o.s)
    def __eq__(self, o): return self.s==o.s and self.e==o.e
    def p(self, rt:float=0.5)->Pos: return self.s + (self.e - self.s)*rt
    def green(self, lo:float=0.0, hi:float=1.0)->float:
        d = hi-lo; m = self.p((lo+hi)*0.5)
        return m.y * d * (self.s.x - self.e.x)

def intersection_param(s1: Seg, s2: Seg, f: int = STRONG) -> float:
    p1,p2,q1,q2 = s1.s, s1.e, s2.s, s2.e
    det = (q2-q1) / (p2-p1)
    if zero(det): return -1.0
    a1 = ((q2-q1) / (q1-p1)) / det
    a2 = ((p2-p1) / (p1-q1)) / -det
    if f == WEAK: return fit(a1,0.0,1.0)
    if (0.0 < a1 < 1.0) and (-TOL < a2 < 1.0+TOL): return a1
    return -1.0

def inner_check(H: Polygon, q: Pos, d: Pos=Pos(0.0,0.0)) -> int:
    # convex (on-edge도 1 취급)
    sz = len(H)
    for i in range(sz):
        p1, p2 = H[i], H[(i+1)%sz]
        if ccw(p1,p2,q) < 0: return 0
        if on_seg_strong(p1,p2,q) and not eq(d.x,d.y):
            return 1 if sign((p2-p1)/d)>0 else 0
    return 1

def inner_check_concave(H: Polygon, p: Pos, s: Pos, e: Pos) -> bool:
    cnt = 0; sz = len(H)
    for i in range(sz):
        cur = H[i]; nxt = H[(i+1)%sz]
        if on_seg_strong(cur,nxt,p):
            assert collinear(cur,nxt,s,e)
            return dot4(cur,nxt,s,e) > 0
        if zero(cur.y - nxt.y): continue
        c, n = (cur,nxt) if nxt.y >= cur.y else (nxt,cur)
        if n.y - TOL < p.y or c.y > p.y: continue
        cnt += ccw(c,n,p) > 0
    return (cnt & 1) == 1

def get_pos(l: Pos, p: Seg, q: Seg) -> Pos:
    p1,p2,q1,q2 = p.s,p.e,q.s,q.e
    if (not inside(p2,l,p1,q1,WEAK)) and (not inside(p2,l,p1,q2,WEAK)):
        if intersect_seg(l,p1,q1,q2) and intersect_seg(l,p2,q1,q2): return Pos(0.0,1.0)
        else: return Pos(0.0,0.0)
    tri = [p1,p2,l]
    in1 = inner_check(tri,q1); in2 = inner_check(tri,q2)
    if (not in1) and (not in2): return Pos(0.0,0.0)
    r1=0.0; r2=1.0
    if in1 and in2:
        r1 = intersection_param(p, Seg(l,q1), WEAK)
        r2 = intersection_param(p, Seg(l,q2), WEAK)
    elif in1: r1 = intersection_param(p, Seg(l,q1), WEAK)
    elif in2: r2 = intersection_param(p, Seg(l,q2), WEAK)
    else: r1=r2=0.0
    if r2<r1: r1,r2=r2,r1
    return Pos(r1,r2)

def intersections_on_edge(l: Seg, H: Polygon):
    sz = len(H); ret=[]
    for i in range(sz):
        p0 = H[i]; p1 = H[(i+1)%sz]
        k = Seg(p0,p1)
        if collinear(l.s,l.e,p0,p1):
            for p in (p0,p1):
                ix = projection(l.s,l.e,p)
                ret.append(fit(ix))
        elif ccw(l.s,l.e,p0)*ccw(l.s,l.e,p1) <= 0:
            ix = intersection_param(l,k,WEAK)
            ret.append(fit(ix))
    ret.sort()
    # unique with epsilon
    out=[]
    for v in ret:
        if not out or not eq(out[-1], v): out.append(v)
    return out

# ---- globals ----
N=M=K=Q=0
A=[0.0]*(1<<3)
I=[0]*(1<<3)
C=[0]*LEN
Lpos=[Pos() for _ in range(1<<3)]
P=[[] for _ in range(LEN)]          # polygons input
Hc=[[] for _ in range(LEN)]         # per-color hull-ish poly
Ttri=[[] for _ in range(1<<3)]      # debug triangles (not printed)
Sseg=[[Seg() for _ in range(LEN)] for _ in range(1<<3)]
Z=[[] for _ in range(1<<3)]
VX=[[[] for _ in range(1<<3)] for _ in range(1<<3)]
cen = Pos()

def cmpt(p: Seg, q: Seg) -> int:
    assert ccw(cen,p.s,p.e)!=0
    assert ccw(cen,q.s,q.e)!=0
    u = p.s - cen; v = q.s - cen
    f1 = (O < u); f2 = (O < v)
    if f1 != f2: return -1 if f1 else 1
    cr = u / v
    assert not zero(cr)
    return -1 if cr > 0 else 1

def query_one():
    global A,I,C,Hc,Z,Sseg,VX,cen,Lpos
    rx,ry = map(float, sys.stdin.readline().split())
    gx,gy = map(float, sys.stdin.readline().split())
    bx,by = map(float, sys.stdin.readline().split())
    Lpos[RED]=Pos(rx,ry); Lpos[GREEN]=Pos(gx,gy); Lpos[BLUE]=Pos(bx,by)

    A=[0.0]*(1<<3)
    I=[-1]*(1<<3)
    for i in range(LEN): C[i]=0
    for c in range(1, (1<<3)): Ttri[c].clear()
    for c in range(1, (1<<3)): Hc[c].clear()
    for c in range(1, (1<<3)): Z[c].clear()

    for c in (RED,GREEN,BLUE):
        l = Lpos[c]
        f0 = inner_check(P[0], l)
        if not f0: continue
        if f0: C[0] |= c; I[c]=0

        fk = 0
        for k in range(1, K+1):
            fk = inner_check(P[k], l)
            if fk:
                C[0] -= c; C[k] |= c; I[c]=k
                break
        if fk: continue

        # support segments for each blue polygon
        for k in range(1, K+1):
            Hk = P[k]; M = len(Hk)
            pl = Hk[0]; pr = Hk[0]
            for j in range(M):
                if ccw(l,pl,Hk[j]) > 0: pl = Hk[j]
                if ccw(l,pr,Hk[j]) < 0: pr = Hk[j]
            Sseg[c][k] = Seg(pr, pl)

        # from P0 edges
        sz = len(P[0])
        for i in range(sz):
            u = P[0][i]; v = P[0][(i+1)%sz]
            assert ccw(l,u,v) > 0
            w = Seg(u,v)
            VP = [Pos(0.0,0.0)]
            for k in range(1, K+1):
                se = get_pos(l,w,Sseg[c][k])
                if not eq(se.x, se.y): VP.append(se)
            VP.append(Pos(1.0,1.0))
            VP.sort()
            hi = 0.0
            for p in VP:
                if hi < p.x:
                    s = w.p(hi); e = w.p(p.x)
                    Ttri[c].append([l,s,e])
                    Z[c].append(Seg(s,e))
                    hi = p.y
                else:
                    hi = max(hi, p.y)

        # from blue polygons' edges
        for k in range(1, K+1):
            Hk = P[k]; sz = len(Hk)
            for i in range(sz):
                u = Hk[i]; v = Hk[(i+1)%sz]
                if ccw(l,u,v) >= 0: continue
                w = Seg(v,u)
                VP=[Pos(0.0,0.0)]
                for kk in range(1, K+1):
                    if kk==k: continue
                    se = get_pos(l,w,Sseg[c][kk])
                    if not eq(se.x,se.y): VP.append(se)
                VP.append(Pos(1.0,1.0))
                VP.sort()
                hi = 0.0
                for p in VP:
                    if hi < p.x:
                        s = w.p(hi); e = w.p(p.x)
                        Ttri[c].append([l,s,e])
                        Z[c].append(Seg(s,e))
                        hi = p.y
                    else:
                        hi = max(hi, p.y)

        cen = l
        Z[c].sort(key=cmp_to_key(cmpt))
        for se in Z[c]:
            Hc[c].append(se.s); Hc[c].append(se.e)
        # unique consecutive
        tmp=[]
        for p in Hc[c]:
            if not tmp or not (p==tmp[-1]): tmp.append(p)
        Hc[c] = tmp
        if Hc[c] and Hc[c][0] == Hc[c][-1]: Hc[c].pop()

    # intersections VX for color pairs
    for i in range(3):
        c1 = (1 << ((i+1)%3))
        c2 = (1 << ((i+2)%3))
        if I[c1]==0 and I[c2]==0:
            for _ in range(2):
                VX[c1][c2].clear()
                sz = len(Hc[c1]); VX[c1][c2] = [[] for _ in range(sz)]
                for j in range(sz):
                    V = [0.0, 1.0]
                    se = Seg(Hc[c1][j], Hc[c1][(j+1)%sz])
                    tmpv = intersections_on_edge(se, Hc[c2])
                    V.extend(tmpv)
                    V.sort()
                    # unique eps
                    out=[]
                    for v in V:
                        if not out or not eq(out[-1],v): out.append(v)
                    VX[c1][c2][j] = out
                c1, c2 = c2, c1

    # R&G&B (white)
    if I[RED]==0 and I[GREEN]==0 and I[BLUE]==0:
        VS=[]
        for i in range(3):
            c0 = (1<<i); c1=(1<<((i+1)%3)); c2=(1<<((i+2)%3))
            sz = len(Hc[c0])
            for j in range(sz):
                V = [0.0,1.0]
                V += VX[c0][c1][j]
                V += VX[c0][c2][j]
                V.sort()
                uniq=[]
                for v in V:
                    if not uniq or not eq(uniq[-1],v): uniq.append(v)
                se0 = Seg(Hc[c0][j], Hc[c0][(j+1)%sz], c0)
                for a,b in zip(uniq, uniq[1:]):
                    VS.append(Seg(se0.p(a), se0.p(b), c0))
        VS.sort()
        # unique consecutive segs
        out=[]
        for s in VS:
            if not out or not (s==out[-1]): out.append(s)
        VS = out
        for se in VS:
            f=True
            for i in range(3):
                c = (1<<i)
                if c==se.i: continue
                m = se.p(0.5)
                if not inner_check_concave(Hc[c], m, se.s, se.e):
                    f=False; break
            if f: A[WHITE] += se.green()

    # pairs
    for i in range(3):
        c1 = (1 << ((i+1)%3))
        c2 = (1 << ((i+2)%3))
        if I[c1]==0 and I[c2]==0:
            c = c1|c2
            VS=[]
            for _ in range(2):
                sz=len(Hc[c1])
                for j in range(sz):
                    se0 = Seg(Hc[c1][j], Hc[c1][(j+1)%sz], c1)
                    V = VX[c1][c2][j][:]
                    V.sort()
                    uniq=[]
                    for v in V:
                        if not uniq or not eq(uniq[-1],v): uniq.append(v)
                    for a,b in zip(uniq, uniq[1:]):
                        VS.append(Seg(se0.p(a), se0.p(b), c1))
                c1,c2 = c2,c1
            VS.sort()
            out=[]
            for s in VS:
                if not out or not (s==out[-1]): out.append(s)
            VS=out
            for se in VS:
                f=True
                for cc in (c1,c2):
                    if cc==se.i: continue
                    m=se.p(0.5)
                    if not inner_check_concave(Hc[cc], m, se.s, se.e):
                        f=False; break
                if f: A[c] += se.green()

    # singles
    for c in (RED,GREEN,BLUE):
        if I[c]==0:
            A[c] += area_poly(Hc[c])

    # blue polygons area add
    for k in range(1, K+1):
        A[C[k]] += area_poly(P[k])

    # inclusion-exclusion
    if I[RED]==0 and I[GREEN]==0 and I[BLUE]==0:
        for c in range(1, WHITE):
            A[c] -= A[WHITE]

    for i in range(3):
        c1=(1<<((i+1)%3)); c2=(1<<((i+2)%3))
        if I[c1]==0 and I[c2]==0:
            c=c1|c2
            A[c1] -= A[c]
            A[c2] -= A[c]

    A[BLACK] = area_poly(P[0])
    for c in range(1,(1<<3)): A[BLACK] -= A[c]
    for c in range(0,(1<<3)):
        if A[c] < 0: A[c]=0.0

    # output
    out=[]
    for c in (RED,GREEN,BLUE,YELLOW,MAGENTA,CYAN,WHITE,BLACK):
        out.append(f"{A[c]:.13f}")
    sys.stdout.write("\n".join(out) + "\n")

def solve():
    global N,K,Q,P
    data = sys.stdin.read().strip().split()
    it = iter(data)
    try:
        N = int(next(it))
    except StopIteration:
        return
    P = [[] for _ in range(LEN)]
    P0=[]
    for _ in range(N):
        x=float(next(it)); y=float(next(it)); P0.append(Pos(x,y))
    P[0]=P0; norm_poly(P[0])
    K = int(next(it))
    for i in range(1, K+1):
        M = int(next(it)); poly=[]
        for _ in range(M):
            x=float(next(it)); y=float(next(it)); poly.append(Pos(x,y))
        norm_poly(poly)
        P[i]=poly
    Q = int(next(it))
    # repackage remaining tokens into lines for query_one
    from io import StringIO
    rest = list(it)
    buf=[]
    for i in range(0, len(rest), 2):
        buf.append(f"{rest[i]} {rest[i+1]}\n")
    sys.stdin = StringIO("".join(buf))
    for _ in range(Q):
        query_one()

if __name__ == "__main__":
    # === 파일 경로만 바꾸세요 ===
    sys.stdin  = open("../tests/candle/in/23.in", "r", encoding="utf-8")
    t0 = time.perf_counter()
    solve()
    dt = time.perf_counter() - t0
    print(f"[wall] {dt*1000:.2f} ms", file=sys.stderr)
