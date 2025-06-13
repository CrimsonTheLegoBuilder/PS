# --- coding: utf-8 ---

with open('aedsis.txt', 'rt') as sis:
	# aiapostide arv, Uno asukoht, aiapostide asukohad
	n = int(sis.readline().strip())
	u = [int(_) for _ in sis.readline().strip().split()]
	p = []
	for i in range(n):
		p.append([int(_) for _ in sis.readline().strip().split()])

# p0->p1 ja q0->q1 vektorkorrutis (täpsemalt selle Z-koordinaat)
def vkorr(p0, p1, q0, q1):
	p0x, p0y = p0
	p1x, p1y = p1
	q0x, q0y = q0
	q1x, q1y = q1
	return (p1x - p0x) * (q1y - q0y) - (p1y - p0y) * (q1x - q0x)

# tegur, millega on vaja korrutada vektorit p0->p1, et selle lõpp
# langeks sirgele q0-q1; teguri tagastame hariliku murruna a/b
def tegur(p0, p1, q0, q1):
	# sirge p0-p1 punktid rahuldavad võrrandit p0 + t * p0->p1
	# sirge q0-q1 punktid rahuldavad võrrandit q0 + u * q0->q1
	# sellest võrrandisüsteemist saame avaldada t = a/b kujul
	a, b = vkorr(p0, q0, q0, q1), vkorr(p0, p1, q0, q1)
	return (a, b) if b >= 0 else (-a, -b)

# tagastab (p0, p1) või (p1, p0) nii, et esimene
# oleks punktist u vaadates vasakul ja teine paremal
def orient(p0, p1):
	if vkorr(u, p0, u, p1) <= 0:
		return (p0, p1)
	else:
		return (p1, p0)

# kas aialõik p[i]-p[i+1] on punktist u nähtav
def nahtav(i):
	# mugavuseks: vp-pp on see lõik, mille nähtavust uurime
	j0, j1 = i, (i + 1) % n
	vp, pp = orient(p[j0], p[j1])
	# vektorid u->vk ja u->pk osutavad vp-pp nähtava osa otsi
	# avp/bvp ja app/bpp on tegurid, millega neid vektoreid
	# tuleks korrutada, et nad ulataks täpselt sirgeni vp-pp
	vk, avp, bvp = vp, 1, 1
	pk, app, bpp = pp, 1, 1
	# lahutame lõigust maha teiste varjus olevad sektorid
	while j1 != i:
		# mugavuseks: vo-po on see lõik, mille kohta nüüd
		# uurime, kas see kitsendab lõigu vp-pp nähtavat osa
		j0, j1 = j1, (j1 + 1) % n
		vo, po = orient(p[j0], p[j1])
		# tegurid, millega tuleks vektoreid u->vk ja u->pk
		# korrutada, et nad ulataks täpselt sirgeni vo-po
		avo, bvo = tegur(u, vk, vo, po)
		apo, bpo = tegur(u, pk, vo, po)
		# kas lõik vo-po lõikab kiiri u->vk ja u->pk
		lv = vkorr(u, vk, u, vo) * vkorr(u, vk, u, po) <= 0 and avo >= 0
		lp = vkorr(u, pk, u, vo) * vkorr(u, pk, u, po) <= 0 and apo >= 0
		if lv and lp:
			# lõikab mõlemat kiirt: läbib sektori kogu laiuses
			# a1/b1 >= a2/b2 asemel kasutame a1*b2 >= a2*b1
			if avp * bvo >= avo * bvp and app * bpo >= apo * bpp:
				# sirgest vp-pp punkti u pool: varjab kogu sektori
				return False
		elif lv:
			# lõikab ainult vasakut kiirt: parem ots sektori sees
			a, b = tegur(u, po, vp, pp)
			if a >= b:
				# sirgest vp-pp punkti u pool: kitsendab sektorit vasakult
				vk = po
				avp, bvp = tegur(u, vk, vp, pp)
		elif lp:
			# lõikab ainult paremat kiirt: vasak ots sektori sees
			a, b = tegur(u, vo, vp, pp)
			if a >= b:
				# sirgest vp-pp punkti u pool: kitsendab sektorit paremalt
				pk = vo
				app, bpp = tegur(u, pk, vp, pp)
		# seda, et lõik vo-po on tervenisti sektori sees, ei saa juhtuda,
		# sest me alustame ühe kiire peal oleva otspunktiga lõigust ja
		# edasi vaatame aialõike läbi piki aeda liikudes; seega, kui lõik
		# kumbagi kiirt ei lõika, on ta tervenisti sektorist väljas
	return True

naeb = []
for i in range(n):
	if nahtav(i):
		naeb.append(1 + i)

with open('aedval.txt', 'wt') as val:
	val.write('{}\n'.format(len(naeb)))
	val.write('{}\n'.format(' '.join([str(_) for _ in naeb])))
