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

]
R = [
  # [
  #   (0.00000000, 40.00000000),
  #   (-100.00000000, -20.00000000),
  #   (-60.00000000, -60.00000000),
  # ],
  [
    (0.00000000, 40.00000000),
    (-60.00000000, -60.00000000),
    (-50.58823529, -61.17647059),
  ],
  [
    (0.00000000, 40.00000000),
    (-32.33082707, -63.45864662),
    (-17.55102041, -65.30612245),
  ],
  [
    (0.00000000, 40.00000000),
    (17.31543624, -69.66442953),
    (20.00000000, -70.00000000),
  ],
  # [
  #   (0.00000000, 40.00000000),
  #   (20.00000000, -70.00000000),
  #   (80.00000000, -30.00000000),
  # ],
  # [
  #   (0.00000000, 40.00000000),
  #   (80.00000000, -30.00000000),
  #   (100.00000000, 40.00000000),
  # ],
  # [
  #   (0.00000000, 40.00000000),
  #   (100.00000000, 40.00000000),
  #   (70.00000000, 80.00000000),
  # ],
  # [
  #   (0.00000000, 40.00000000),
  #   (70.00000000, 80.00000000),
  #   (-30.00000000, 60.00000000),
  # ],
  # [
  #   (0.00000000, 40.00000000),
  #   (-30.00000000, 60.00000000),
  #   (-80.00000000, 20.00000000),
  # ],
  # [
  #   (0.00000000, 40.00000000),
  #   (-80.00000000, 20.00000000),
  #   (-100.00000000, -20.00000000),
  # ],
  # [
  #   (0.00000000, 40.00000000),
  #   (-25.00000000, -20.00000000),
  #   (-30.00000000, -40.00000000),
  # ],
  # [
  #   (0.00000000, 40.00000000),
  #   (-40.00000000, -40.00000000),
  #   (-25.00000000, -20.00000000),
  # ],
  # [
  #   (0.00000000, 40.00000000),
  #   (-31.66666667, -44.44444444),
  #   (-25.00000000, -40.00000000),
  # ],
  # [
  #   (0.00000000, 40.00000000),
  #   (-5.00000000, 10.00000000),
  #   (15.00000000, -55.00000000),
  # ]
]
G = [

]
Y = [

]
B = [

]
M = [

]
C = [

]
W = [

]

grey = [
  [
    (0.00000000, 40.00000000),
    (-50.58823529, -61.17647059),
    (-38.50746269, -62.68656716),
  ],
  [
    (0.00000000, 40.00000000),
    (-45.26315789, -61.84210526),
    (-32.33082707, -63.45864662),
  ],
  [
    (0.00000000, 40.00000000),
    (-17.55102041, -65.30612245),
    (17.31543624, -69.66442953),
  ]
]

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
    sys.stdin = open('D:\\PS\\kitpc2025\\tests\\candle_renew\\in\\101.in', 'r')

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
                ax.fill(x, y, facecolor=COLOR_MAP[color_name], edgecolor='none', antialiased=False, zorder=0,alpha=0.5)

    color = ['red', 'blue']
    for i, polygon in enumerate(H):
        x, y = zip(*polygon)
        ax.fill(x, y, edgecolor=color[0 if not i else 1], fill=False, linewidth=1, zorder=1)

    ax.set_aspect('equal')
    ax.margins(x=1, y=1)
    ax.set_xlim(-105, 105)
    ax.set_ylim(-75, 85)
    ax.axis('off')  # 축 숨기기
    plt.show()