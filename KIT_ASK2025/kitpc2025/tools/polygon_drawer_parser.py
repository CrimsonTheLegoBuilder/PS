#!/usr/bin/env python3
# draw_polys.py
import argparse, ast, math, sys
from pathlib import Path
import numpy as np
import matplotlib
# 대용량 렌더링 최적화
matplotlib.rcParams['agg.path.chunksize'] = 20000  # 아주 큰 경로 분할 렌더링
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

# 파일에 있는 키들: L,R,G,B,Y,M,C,W  -> 색상
COLOR_MAP = {
    'L': 'gray',     # BLACK로 쓰시던 영역
    'R': 'red',
    'G': 'lime',
    'B': 'blue',
    'Y': 'yellow',
    'M': 'magenta',
    'C': 'cyan',
    'W': 'white',
}

DRAW_ORDER = ['L', 'R', 'G', 'B', 'Y', 'M', 'C', 'W']  # 겹침 순서

def parse_literal_assignments(text: str):
    """
    파일 전체를 AST로 파싱하여 L/R/G/B/Y/M/C/W 변수에 대입된 '리터럴(list/tuple/float/int)'만 읽어온다.
    exec 없이 안전하게 처리.
    """
    tree = ast.parse(text, filename="<polys>", mode="exec")
    out = {k: [] for k in DRAW_ORDER}
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        # 다중 대입도 대비
        for tgt in node.targets:
            if isinstance(tgt, ast.Name) and tgt.id in out:
                # 값이 전부 리터럴이어야 함
                try:
                    val = ast.literal_eval(node.value)
                except Exception as e:
                    raise ValueError(f"{tgt.id} 에 literal이 아닌 값이 있어 파싱 실패: {e}")
                # 기대 형태: [ [ (x,y), ... ], [ ... ], ... ]
                if not isinstance(val, (list, tuple)):
                    raise ValueError(f"{tgt.id}는 리스트여야 합니다.")
                # numpy PolyCollection용으로 float32/64 배열로 정리
                cleaned = []
                for poly in val:
                    if not isinstance(poly, (list, tuple)):
                        continue
                    # 좌표쌍만 받음
                    pts = []
                    for p in poly:
                        if (isinstance(p, (list, tuple)) and len(p) == 2 and
                                isinstance(p[0], (int, float)) and isinstance(p[1], (int, float))):
                            pts.append((float(p[0]), float(p[1])))
                    if len(pts) >= 3:  # 면만 그린다(점/선 스킵)
                        cleaned.append(np.asarray(pts, dtype=np.float32))
                out[tgt.id] = cleaned
    return out

def compute_bounds(groups):
    xs_min, ys_min = math.inf, math.inf
    xs_max, ys_max = -math.inf, -math.inf
    for polys in groups.values():
        for a in polys:
            if a.size == 0: 
                continue
            x = a[:,0]; y = a[:,1]
            xs_min = min(xs_min, float(x.min()))
            ys_min = min(ys_min, float(y.min()))
            xs_max = max(xs_max, float(x.max()))
            ys_max = max(ys_max, float(y.max()))
    if not math.isfinite(xs_min):  # 비어있을 때
        xs_min, ys_min, xs_max, ys_max = -1.0, -1.0, 1.0, 1.0
    return xs_min, ys_min, xs_max, ys_max

def main():
    ap = argparse.ArgumentParser(description="다각형 데이터(L/R/G/B/Y/M/C/W) 파일을 그려 PNG/PDF로 저장")
    ap.add_argument("input", type=Path, help="다각형 데이터 텍스트(파이썬 리터럴 할당 형식)")
    ap.add_argument("-o", "--output", type=Path, default=None, help="저장 경로(.png/.pdf). 미지정 시 화면 표시")
    ap.add_argument("--dpi", type=int, default=200, help="저장 DPI")
    ap.add_argument("--margin", type=float, default=0.05, help="오토 범위 여백 비율(0.05=5%)")
    ap.add_argument("--xlim", nargs=2, type=float, help="x축 범위 수동 지정: xmin xmax")
    ap.add_argument("--ylim", nargs=2, type=float, help="y축 범위 수동 지정: ymin ymax")
    ap.add_argument("--no-aa", action="store_true", help="안티에일리어싱 끄기(더 빠름)")
    args = ap.parse_args()

    text = args.input.read_text(encoding="utf-8")
    groups = parse_literal_assignments(text)

    fig, ax = plt.subplots(figsize=(8, 8))
    # 색상별로 PolyCollection 추가(일괄 드로잉)
    for key in DRAW_ORDER:
        polys = groups.get(key, [])
        if not polys:
            continue
        coll = PolyCollection(
            polys,
            facecolors=COLOR_MAP[key],
            edgecolors='none',
            antialiased=not args.no_aa,
            closed=True
        )
        ax.add_collection(coll)

    # 자동 범위 또는 수동 범위
    if args.xlim and args.ylim:
        ax.set_xlim(*args.xlim)
        ax.set_ylim(*args.ylim)
    else:
        xmin, ymin, xmax, ymax = compute_bounds(groups)
        dx = xmax - xmin; dy = ymax - ymin
        mx = dx * args.margin if dx > 0 else 1.0
        my = dy * args.margin if dy > 0 else 1.0
        ax.set_xlim(xmin - mx, xmax + mx)
        ax.set_ylim(ymin - my, ymax + my)

    ax.set_aspect('equal')
    ax.axis('off')

    if args.output:
        # 포맷은 확장자 보고 자동 결정
        fig.savefig(args.output, dpi=args.dpi, bbox_inches='tight', pad_inches=0.0)
        print(f"saved: {args.output}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
