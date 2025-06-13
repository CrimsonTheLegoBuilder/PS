#include <fstream>
#include <vector>
#include <tuple>
using namespace std;

ifstream sis("aedsis.txt");
ofstream val("aedval.txt");

const double EPS = 1e-6; // täpsus: kõik alla selle on null

// a < b asemel
bool LT(double a, double b) {
	return a < b - EPS;
}

// a <= b asemel
bool LE(double a, double b) {
	return a <= b + EPS;
}

typedef struct {
	double x, y;
} point;

int n; // aiapostide arv
point u; // Uno asukoht
vector<point> p; // aiapostide asukohad

// p0->p1 ja q0->q1 vektorkorrutis (täpsemalt selle Z-koordinaat)
double vkorr(point p0, point p1, point q0, point q1) {
	return (p1.x - p0.x) * (q1.y - q0.y) - (p1.y - p0.y) * (q1.x - q0.x);
}

// tegur, millega on vaja korrutada vektorit p0->p1, et selle lõpp
// langeks sirgele q0-q1; teguri tagastame hariliku murruna a/b
double tegur(point p0, point p1, point q0, point q1) {
	// sirge p0-p1 punktid rahuldavad võrrandit p0 + t * p0->p1
	// sirge q0-q1 punktid rahuldavad võrrandit q0 + u * q0->q1
	// sellest võrrandisüsteemist saame avaldada t
	return vkorr(p0, q0, q0, q1) / vkorr(p0, p1, q0, q1);
}

// järjestab p0 ja p1 nii, et esimene oleks punktist u vaadates
// vasakul ja teine paremal
void orient(point & p0, point & p1) {
	if (vkorr(u, p0, u, p1) > 0) {
		swap(p0, p1);
	}
}

// kas aialõik p[i]-p[i+1] on punktist u nähtav
bool nahtav(int i) {
	// mugavuseks: vp-pp on see lõik, mille nähtavust uurime
	int j0 = i, j1 = (i + 1) % n;
	point vp = p[j0], pp = p[j1];
	orient(vp, pp);
	// vektorid u->vk ja u->pk osutavad vp-pp nähtava osa otsi
	// tvp ja tpp on tegurid, millega neid vektoreid tuleks
	// korrutada, et nad ulataks täpselt sirgeni vp-pp
	point vk = vp, pk = pp;
	double tvp = 1, tpp = 1;
	// lahutame lõigust maha teiste varjus olevad sektorid
	while (j1 != i) {
		// mugavuseks: vo-po on see lõik, mille kohta nüüd
		// uurime, kas see kitsendab lõigu vp-pp nähtavat osa
		j0 = j1; j1 = (j1 + 1) % n;
		point vo = p[j0], po = p[j1];
		orient(vo, po);
		// tegurid, millega tuleks vektoreid u->vk ja u->pk
		// korrutada, et nad ulataks täpselt sirgeni vo-po
		double tvo = tegur(u, vk, vo, po), tpo = tegur(u, pk, vo, po);
		// kas lõik vo-po lõikab kiiri u->vk ja u->pk
		bool lv = LE(vkorr(u, vk, u, vo) * vkorr(u, vk, u, po), 0) && LE(0, tvo);
		bool lp = LE(vkorr(u, pk, u, vo) * vkorr(u, pk, u, po), 0) && LE(0, tpo);
		if (lv && lp) {
			// lõikab mõlemat kiirt: läbib sektori kogu laiuses
			if (LE(tvo, tvp) && LE(tpo, tpp)) {
				// sirgest vp-pp punkti u pool: varjab kogu sektori
				return false;
			}
		} else if (lv) {
			// lõikab ainult vasakut kiirt: parem ots sektori sees
			if (LT(1, tegur(u, po, vp, pp))) {
				// sirgest vp-pp punkti u pool: kitsendab sektorit vasakult
				vk = po; tvp = tegur(u, vk, vp, pp);
			}
		} else if (lp) {
			// lõikab ainult paremat kiirt: vasak ots sektori sees
			if (LT(1, tegur(u, vo, vp, pp))) {
				// sirgest vp-pp punkti u pool: kitsendab sektorit paremalt
				pk = vo; tpp = tegur(u, pk, vp, pp);
			}
		}
		// seda, et lõik vo-po on tervenisti sektori sees, ei saa juhtuda,
		// sest me alustame ühe kiire peal oleva otspunktiga lõigust ja
		// edasi vaatame aialõike läbi piki aeda liikudes; seega, kui lõik
		// kumbagi kiirt ei lõika, on ta tervenisti sektorist väljas
	}
	return true;
}

int main() {
	sis >> n;
	sis >> u.x >> u.y;
	p.resize(n);
	for (int i = 0; i < n; ++i) {
		sis >> p[i].x >> p[i].y;
	}

	vector<int> naeb;
	for (int i = 0; i < n; ++i) {
		if (nahtav(i)) {
			naeb.push_back(1 + i);
		}
	}

	val << naeb.size() << '\n';
	for (int i = 0; i < naeb.size(); ++i) {
		if (i > 0) {
			val << ' ';
		}
		val << naeb[i];
	}
	val << '\n';
}
