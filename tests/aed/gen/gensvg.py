# sisendi visualiseerija

import sys

S = 20

fi = open(sys.argv[1], 'rt') if len(sys.argv) > 1 else sys.stdin
fo = open(sys.argv[2], 'wt') if len(sys.argv) > 2 else sys.stdout

n = int(fi.readline().strip())

s = fi.readline().strip().split()
ux = int(s[0])
uy = -int(s[1])

px, py = [], []
for i in range(n):
	s = fi.readline().strip().split()
	px.append(int(s[0]))
	py.append(-int(s[1]))

minx = min(px + [ux]) - 1
maxx = max(px + [ux]) + 1
miny = min(py + [uy]) - 1
maxy = max(py + [uy]) + 1

ppx = sorted(px)
ppy = sorted(py)

fo.write('<?xml version="1.0" encoding="UTF-8"?>\n')
fo.write('<svg xmlns="http://www.w3.org/2000/svg" viewBox="{} {} {} {}">\n'.format(S * (minx), S * (miny), S * (maxx - minx), S * (maxy - miny)))

for i in range(n):
	if ppx[i] != ppx[i - 1] and ppx[i] != ux:
		fo.write('  <line x1="{}" y1="{}" x2="{}" y2="{}" stroke="gray" />\n'.format(S * ppx[i], S * miny, S * ppx[i], S * maxy))
	if ppy[i] != ppy[i - 1] and ppy[i] != uy:
		fo.write('  <line x1="{}" y1="{}" x2="{}" y2="{}" stroke="gray" />\n'.format(S * minx, S * ppy[i], S * maxx, S * ppy[i]))

fo.write('  <line x1="{}" y1="{}" x2="{}" y2="{}" stroke="black" />\n'.format(S * ux, S * miny, S * ux, S * maxy))
fo.write('  <line x1="{}" y1="{}" x2="{}" y2="{}" stroke="black" />\n'.format(S * minx, S * uy, S * maxx, S * uy))

fo.write('  <circle cx="{}" cy="{}" r="5" fill="red" />\n'.format(S * ux, S * uy))

for i in range(n):
	fo.write('  <circle cx="{}" cy="{}" r="3" fill="blue" />\n'.format(S * px[i], S * py[i]))
	fo.write('  <line x1="{}" y1="{}" x2="{}" y2="{}" stroke="blue" />\n'.format(S * px[i - 1], S * py[i - 1], S * px[i], S * py[i]))
	#fo.write('  <text x="{}" y="{}" alignment-baseline="middle" text-anchor="middle">{}</text>\n'.format(S / 2 * (px[i - 1] + px[i]) , S / 2 * (py[i - 1] + py[i]), 1 + i))

fo.write('</svg>\n')

fi.close()
fo.close()
