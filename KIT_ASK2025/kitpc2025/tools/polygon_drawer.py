import os
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
# plt.ion()  # 인터랙티브 모드


if __name__ == "__main__":
    sys.stdin = open('D:\\PS\\kitpc2025\\tests\\candle\\in\\18.in', 'r')

    N: int = int(sys.stdin.readline())
    H: list = []
    P: list = [tuple(map(int, sys.stdin.readline().split())) for _ in range(N)]
    H.append(P)
    K: int = int(sys.stdin.readline())
    for _ in range(K):
        M: int = int(sys.stdin.readline())
        P: list = [tuple(map(int, sys.stdin.readline().split())) for _ in range(M)]
        H.append(P)

    color = ['red', 'blue']

    fig, ax = plt.subplots(figsize=(8, 8))
    for i, polygon in enumerate(H):
        x, y = zip(*polygon)
        ax.fill(x, y, edgecolor=color[0 if not i else 1], fill=False)

    ax.set_aspect('equal')
    ax.margins(x=0.2, y=0.2)
    ax.grid(True)

    plt.show(block=True)