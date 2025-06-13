# sisendi korrektsuse kontrollija

import sys

fi = open(sys.argv[1], 'rt')
fo = open(sys.argv[2], 'wt')

n = int(fi.readline().strip())
if 3 > n or n > 1000:
	sys.exit('N out of range')
fo.write('{}\n'.format(n))

x, y = [int(_) for _ in fi.readline().strip().split()]
if -1e5 > x or x > 1e5:
	sys.exit('X out of range')
if -1e5 > y or y > 1e5:
	sys.exit('Y out of range')
fo.write('{} {}\n'.format(x, y))
p = (x, y)

pp = []
for i in range(n):
	x, y = [int(_) for _ in fi.readline().strip().split()]
	if -1e5 > x or x > 1e5:
		sys.exit('X{} out of range'.format(1 + i))
	if -1e5 > y or y > 1e5:
		sys.exit('Y{} out of range'.format(1 + i))
	pp.append((x, y))
	fo.write('{} {}\n'.format(x, y))

def cp(p0, p1, p2):
	x0, y0 = p0
	x1, y1 = p1
	x2, y2 = p2
	dx1, dy1 = x1 - x0, y1 - y0
	dx2, dy2 = x2 - x0, y2 - y0
	return dx1 * dy2 - dx2 * dy1

def dp(p0, p1, p2):
	x0, y0 = p0
	x1, y1 = p1
	x2, y2 = p2
	dx1, dy1 = x1 - x0, y1 - y0
	dx2, dy2 = x2 - x0, y2 - y0
	return dx1 * dx2 + dy1 * dy2

for i0 in range(n):
	i1 = (i0 + 1) % n
	if cp(p, pp[i0], pp[i1]) == 0:
		sys.exit('P on P{}-P{}'.format(1 + i0, 1 + i1))
	for j in range(n):
		if j == i0 or j == i1:
			continue
		if cp(pp[j], pp[i0], pp[i1]) == 0 and dp(pp[j], pp[i0], pp[i1]) <= 0:
			sys.exit('P{} on P{}-P{}'.format(1 + j, 1 + i0, 1 + i1))

for i0 in range(n):
	i1 = (i0 + 1) % n
	for j0 in range(i0 - 1):
		j1 = j0 + 1
		if cp(pp[i0], pp[i1], pp[j0]) * cp(pp[i0], pp[i1], pp[j1]) < 0 and \
				cp(pp[j0], pp[j1], pp[i0]) * cp(pp[j0], pp[j1], pp[i1]) < 0:
			sys.exit('P{}-P{} and P{}-P{} intersect'.format(1 + j0, 1 + j1, 1 + i0, 1 + i1))

cvx = True
for i0 in range(n):
	i1, i2 = (i0 + 1) % n, (i0 + 2) % n
	if cp(pp[i0], pp[i1], pp[i2]) <= 0:
		cvx = False

prp = True
for i0 in range(n):
	i1 = (i0 + 1) % n
	x0, y0 = pp[i0]
	x1, y1 = pp[i1]
	dx, dy = x1 - x0, y1 - y0
	if dx != 0 and dy != 0:
		prp = False

fi.close()
fo.close()

if len(sys.argv) > 3 and sys.argv[3] == '-v':
	s = []
	if cvx and n <= 50:
		s.append(1)
	if cvx:
		s.append(2)
	if prp and n <= 50:
		s.append(3)
	if prp:
		s.append(4)
	if n <= 50:
		s.append(5)
	if n <= 200:
		s.append(6)
	s.append(7)
	print('subtasks: {0}'.format(', '.join([str(_) for _ in s])))
