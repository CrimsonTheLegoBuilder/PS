import os
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# 색상 매핑
COLOR_MAP = {
    'BLACK': 'gray',
    'RED': 'red',
    'GREEN': 'green',
    'BLUE': 'blue',
    'YELLOW': 'yellow',
    'MAGENTA': 'magenta',
    'CYAN': 'cyan',
    'WHITE': 'white',
}

# 추가로 색칠할 폴리곤들

L = [
      [
    (0.00000000, 0.00000000),
    (10.00000000, 0.00000000),
    (10.00000000, 4.00000000),
    (0.00000000, 4.00000000)
  ]
]
R = [
      [
    (1.00000000, 1.00000000),
    (3.00000000, 1.00000000),
    (3.00000000, 3.00000000),
    (1.00000000, 3.00000000)
  ]
]
G = [
      [
    (4.00000000, 1.00000000),
    (6.00000000, 1.00000000),
    (6.00000000, 3.00000000),
    (4.00000000, 3.00000000)
  ]
]
B = [
      [
    (7.00000000, 1.00000000),
    (9.00000000, 1.00000000),
    (9.00000000, 3.00000000),
    (7.00000000, 3.00000000)
  ]
]
C = []
M = []
Y = []
W = []


# 그릴 순서
draw_order = [
    ('BLACK', L),
    ('RED', R),
    ('GREEN', G),
    ('BLUE', B),
    ('YELLOW', Y),
    ('MAGENTA', M),
    ('CYAN', C),
    ('WHITE', W),
]

if __name__ == "__main__":
    sys.stdin = open('D:\\PS\\kitpc2025\\tests\\candle\\in\\01.in', 'r')

    N: int = int(sys.stdin.readline())
    H: list = []
    P: list = [tuple(map(int, sys.stdin.readline().split())) for _ in range(N)]
    H.append(P)
    K: int = int(sys.stdin.readline())
    for _ in range(K):
        M_: int = int(sys.stdin.readline())
        P: list = [tuple(map(int, sys.stdin.readline().split())) for _ in range(M_)]
        H.append(P)

    fig, ax = plt.subplots(figsize=(8, 8))

    for color_name, polygons in draw_order:
        for poly in polygons:
            if len(poly) >= 3:
                # print(poly)
                x, y = zip(*poly)
                ax.fill(x, y, facecolor=COLOR_MAP[color_name], edgecolor='none', antialiased=False)

    color = ['red', 'blue']
    for i, polygon in enumerate(H):
        x, y = zip(*polygon)
        ax.fill(x, y, edgecolor=color[0 if not i else 1], fill=False, linewidth=3)

    ax.set_aspect('equal')
    ax.margins(x=1, y=1)
    ax.set_xlim(-1, 11)
    ax.set_ylim(-1, 5)
    ax.axis('off')  # 축 숨기기
    plt.show()