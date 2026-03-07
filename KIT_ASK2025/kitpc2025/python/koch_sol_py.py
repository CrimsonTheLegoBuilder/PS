# Python 3
import sys
sys.setrecursionlimit(1_000_000)

LEN = 10_000  # 10000 x 10000

def sign(x: int) -> int:
    return -1 if x < 0 else (1 if x > 0 else 0)

class Pos:
    __slots__ = ("x", "y")
    def __init__(self, x=0, y=0):
        self.x = x; self.y = y
    def __eq__(self, o): return self.x == o.x and self.y == o.y
    def __add__(self, o): return Pos(self.x + o.x, self.y + o.y)
    def __sub__(self, o): return Pos(self.x - o.x, self.y - o.y)
    def __mul__(self, n: int): return Pos(self.x * n, self.y * n)

# 6방향 (원본과 동일)
DIR = [Pos(1,0), Pos(1,1), Pos(-1,1), Pos(-1,0), Pos(-1,-1), Pos(1,-1)]

# 전역 출력 경계
sx = LEN // 2
sy = LEN // 2
ex = 0
ey = 0

# 10000 x 10000 배열: 각 행을 bytearray(0/1)로
B = [bytearray(LEN) for _ in range(LEN)]

def mark(y: int, x: int):
    global sx, sy, ex, ey
    if x < 0 or x >= LEN or y < 0 or y >= LEN:
        return  # 안전장치(정상 입력이면 안 걸림)
    if x < sx: sx = x
    if y < sy: sy = y
    if x > ex: ex = x
    if y > ey: ey = y
    B[y][x] = 1

def draw_segment(ps: Pos, pe: Pos, d: int):
    """리프 세그먼트 그리기. C++ 동작과 동일:
       - 수평(d==0,3): 양 끝점 포함
       - 대각(그 외): while(fx--)로 끝점 미포함
    """
    global sx, sy, ex, ey
    # 경계 갱신(끝점 포함)
    sx = min(sx, ps.x, pe.x)
    ex = max(ex, ps.x, pe.x)
    sy = min(sy, ps.y, pe.y)
    ey = max(ey, ps.y, pe.y)

    if d == 0 or d == 3:
        # 수평
        assert ps.y == pe.y
        y = ps.y
        x0, x1 = sorted((ps.x, pe.x))
        # 끝점 포함
        for x in range(x0, x1 + 1):
            mark(y, x)
    else:
        # 45도 대각(끝점 미포함)
        fx = abs(ps.x - pe.x)
        fy = abs(ps.y - pe.y)
        assert fx == fy
        dx = sign(pe.x - ps.x)
        dy = sign(pe.y - ps.y)
        x, y = ps.x, ps.y
        while fx > 0:
            mark(y, x)
            x += dx; y += dy
            fx -= 1

def recur(ps: Pos, pe: Pos, d: int, s: int, SV, SH):
    if s == 1:
        draw_segment(ps, pe, d)
        return

    # 같은 방향으로 3등분
    pv = DIR[d]
    itx = (SV[s-1] - 1) if pv.y else (SH[s-1] - 1)
    ity = (SV[s-1] - 1) if pv.y else 0
    pv = Pos(pv.x * itx, pv.y * ity)

    p0 = ps
    p1 = ps + pv
    p2 = ps + pv * 2
    p3 = ps + pv * 3
    assert p3 == pe

    recur(p0, p1, d, s - 1, SV, SH)
    recur(p2, p3, d, s - 1, SV, SH)

    # 꺾이는 두 조각
    pv2 = DIR[(d + 5) % 6]
    itx2 = (SV[s-1] - 1) if pv2.y else (SH[s-1] - 1)
    ity2 = (SV[s-1] - 1) if pv2.y else 0
    pv2 = Pos(pv2.x * itx2, pv2.y * ity2)
    pm = p1 + pv2

    recur(p1, pm, (d + 5) % 6, s - 1, SV, SH)
    recur(pm, p2, (d + 1) % 6, s - 1, SV, SH)

def solve():
    global sx, sy, ex, ey, B

    line = sys.stdin.readline()
    if not line:
        return
    N = int(line.strip())

    if N == 1:
        sys.stdout.write("  *\n")
        sys.stdout.write(" * *\n")
        sys.stdout.write("*****\n")
        return

    MAXN = max(12, N + 1)
    SV = [0] * MAXN
    SH = [0] * MAXN
    SV[1], SH[1] = 3, 5
    for i in range(2, MAXN):
        SV[i] = SV[i - 1] * 3 - 2
        SH[i] = SH[i - 1] * 3 - 2

    # 중심 초기화
    sx = LEN >> 1
    sy = LEN >> 1
    ex = 0
    ey = 0
    # (B는 이미 0으로 초기화된 상태)

    # 꼭짓점 계산(원본 동일)
    p0 = Pos(sx, sy) - Pos(SH[N] // 2, SV[N - 1] // 2)
    p1 = p0 + Pos(SH[N] - 1, 0)
    p2 = p0 + Pos(SV[N] - 1, SV[N] - 1)

    recur(p0, p1, 0, N, SV, SH)
    recur(p1, p2, 2, N, SV, SH)
    recur(p2, p0, 4, N, SV, SH)

    # 0->' ', 1->'*' 변환 테이블 (256바이트)
    trans = bytes.maketrans(bytes(range(256)),
                            bytes([(32 if v == 0 else (42 if v == 1 else 63)) for v in range(256)]))
    out = []
    write = out.append

    # sy..ey 각 행에 대해, 마지막 '*' 위치까지 출력 (뒤 공백 제거)
    for y in range(sy, ey + 1):
        row = B[y]
        # ex부터 왼쪽으로 마지막 1 찾기 (없으면 -1)
        i = row.rfind(1, 0, ex + 1)
        if i == -1 or i < sx:
            write("\n")
            continue
        segment = bytes(row[sx:i + 1]).translate(trans)
        write(segment.decode("ascii") + "\n")

    sys.stdout.writelines(out)

if __name__ == "__main__":
    solve()
