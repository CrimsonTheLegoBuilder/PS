import sys
from fractions import Fraction
from functools import cmp_to_key
INPUT = sys.stdin.readline


R: int = 1 << 0
G: int = 1 << 1
B: int = 1 << 2
Y: int = R | G
C: int = G | B
M: int = R | B
W: int = R | G | B


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
        
    