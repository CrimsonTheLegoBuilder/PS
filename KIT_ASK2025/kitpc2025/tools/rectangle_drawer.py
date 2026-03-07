import matplotlib.pyplot as plt
import matplotlib.patches as patches

# 웹페이지 크기
W, H, N = map(int, input().split())

# 이미 게재된 광고 목록 (X_lu, Y_lu, X_rd, Y_rd)
ads = []
for _ in range(N):
    x1, y1, x2, y2 = map(int, input().split())
    ads.append((x1, y1, x2, y2))

max_ad = tuple(map(int, input().split()))

fig, ax = plt.subplots()

# 웹페이지 테두리
webpage_border = patches.Rectangle(
    (0, 0), W, H,
    linewidth=2, edgecolor='black', facecolor='none'
)
ax.add_patch(webpage_border)

# 광고 배너들 (반투명한 사각형)
for ad in ads:
    x1, y1, x2, y2 = ad
    ad_rect = patches.Rectangle(
        (x1, y1), x2 - x1, y2 - y1,
        linewidth=1, edgecolor='blue', facecolor='blue', alpha=0.4
    )
    ax.add_patch(ad_rect)

x1, y1, x2, y2 = max_ad
max_ad_rect = patches.Rectangle(
    (x1, y1), x2 - x1, y2 - y1,
    linewidth=2, edgecolor='red', facecolor='none', linestyle='--'
)
ax.add_patch(max_ad_rect)

# 축 설정
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.set_aspect('equal')
plt.gca().invert_yaxis()  # 웹처럼 y축 아래로 증가
plt.title("test case 15")
plt.grid(True)
plt.savefig("awesome_15.png", dpi=300, bbox_inches='tight')
plt.show()

