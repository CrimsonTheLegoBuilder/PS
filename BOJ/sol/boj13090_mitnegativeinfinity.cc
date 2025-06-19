#include <bits/stdc++.h>
using namespace std;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()
typedef pair<int, int> pii;
typedef vector<int> vi;

template <class T> int sgn(T x) { return (x > 0) - (x < 0); }
template<class T>
struct Point {
	typedef Point P;
	T x, y;
	explicit Point(T _x=0, T _y=0) : x(_x), y(_y) {}
	bool operator<(P p) const { return tie(x,y) < tie(p.x,p.y); }
	bool operator==(P p) const { return tie(x,y)==tie(p.x,p.y); }
	P operator+(P p) const { return P(x+p.x, y+p.y); }
	P operator-(P p) const { return P(x-p.x, y-p.y); }
	P operator*(T d) const { return P(x*d, y*d); }
	P operator/(T d) const { return P(x/d, y/d); }
	T dot(P p) const { return x*p.x + y*p.y; }
	T cross(P p) const { return x*p.y - y*p.x; }
	T cross(P a, P b) const { return (a-*this).cross(b-*this); }
	T dist2() const { return x*x + y*y; }
	long double dist() const { return sqrt((long double)dist2()); }
	// angle to x-axis in interval [-pi, pi]
	long double angle() const { return atan2(y, x); }
	P unit() const { return *this/dist(); } // makes dist()=1
	P perp() const { return P(-y, x); } // rotates +90 degrees
	P normal() const { return perp().unit(); }
	// returns point rotated 'a' radians ccw around the origin
	P rotate(long double a) const {
		return P(x*cos(a)-y*sin(a),x*sin(a)+y*cos(a)); }
	friend ostream& operator<<(ostream& os, P p) {
		return os << "(" << p.x << "," << p.y << ")"; }
};

template<class P> vector<P> segInter(P a, P b, P c, P d) {
	auto oa = c.cross(d, a), ob = c.cross(d, b),
	     oc = a.cross(b, c), od = a.cross(b, d);
	// Checks if intersection is single non-endpoint point.
	if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
		return {(a * ob - b * oa) / (ob - oa)};
	set<P> s;
	return {all(s)};
}

using ll = int64_t;
using ld = long double;
using Pll = Point<ll>;
using Pld = Point<ld>;

struct frac {
	ll n, d;
	frac (ll _n, ll _d) : n(_n), d(_d) {}
	friend bool operator < (const frac& a, const frac& b) {
		return __int128_t(a.n) * b.d < __int128_t(b.n) * a.d;
	}
	friend bool operator == (const frac& a, const frac& b) {
		return __int128_t(a.n) * b.d == __int128_t(b.n) * a.d;
	}
};

frac squareSegDist(Pll& s, Pll& e, Pll& p) {
	if (s==e) return frac((p-s).dist2(), 1);
	auto d = (e-s).dist2(), t = min(d,max((ll)0,(p-s).dot(e-s)));
	return frac(((p-s)*d-(e-s)*t).dist2(), d * d);
}

ld segDist(Pld& s, Pld& e, Pld& p) {
	if (s==e) return (p-s).dist();
	auto d = (e-s).dist2(), t = min(d,max(ld(0),(p-s).dot(e-s)));
	return ((p-s)*d-(e-s)*t).dist()/d;
}

ld lineDist(const Pld& a, const Pld& b, const Pld& p) {
	return (ld)(b-a).cross(p-a)/(b-a).dist();
}


Pld closestPoint(Pld s, Pld e, Pld p){
	if (s == e) return s;
	return s + (e-s) * max(ld(0), min(ld(1), (p-s).dot(e-s) / (e-s).dot(e-s)));
}

Pld toPld(Pll p){
	return Pld(p.x, p.y);
}

ld EPS = 1e-9; // for distance

template <class P>
bool sameDir(P s, P t) {
	return abs(s.unit().cross(t.unit())) <= 1e-11 && s.dot(t) > 0;
}
// checks 180 <= s..t < 360?
template <class P>
bool isReflex(P s, P t) {
	auto c = s.cross(t);
	return c ? (c < 0) : (s.dot(t) < 0);
}
// operator < (s,t) for angles in [base,base+2pi)
template <class P>
bool angleCmp(P base, P s, P t) {
	int r = isReflex(base, s) - isReflex(base, t);
	return r ? (r < 0) : (0 < s.cross(t));
}
// is x in [s,t] taken ccw? 1/0/-1 for in/border/out
template <class P>
int angleBetween(P s, P t, P x) {
	if (sameDir(x, s) || sameDir(x, t)) return 0;
	return angleCmp(s, x, t) ? 1 : -1;
}

// checks whether s to e is completely inside the polygon, s and e are guaranteed to be inside the polygon
bool segInsidePolygon(Pld s, Pld e, vector<Pld> poly){
	if((s - e).dist() <= EPS) return true;
	// check if s->e is outside some edge
	for(int i = 0; i < (int)poly.size(); i++){
		Pld pa = poly[i];
		Pld pb = poly[(i + 1) % poly.size()];
		// left of pa -> pb = inside
		{
			ld da = lineDist(s, e, pa);
			ld db = lineDist(s, e, pb);
			if((da >= 0) == (db >= 0)) continue;
			if(abs(da) <= EPS/2 || abs(db) <= EPS/2) continue;
		}
		{
			Pld p1 = s;
			Pld p2 = e;
			if((p2 - p1).cross(pb - pa) < 0) swap(p2, p1);
			//   pb
			// p1  p2
			//   pa
			ld dist = (p2 - p1).dist();
			ld loc_scaled = (p2 - pa).cross(pb - pa) / (p2 - p1).cross(pb - pa) * dist;
			if(loc_scaled >= EPS/2 && loc_scaled <= dist + EPS/2){
				return false;
			}
		}
	}
	// check if s->e cuts through some vertex
	for(int i = 0; i < (int)poly.size(); i++){
		Pld p = poly[i];
		Pld pprv = poly[(i-1+poly.size()) % poly.size()];
		Pld pnxt = poly[(i+1) % poly.size()];
		// pnxt -> pprv is the counterclockwise
		if(segDist(s, e, p) > EPS) continue;
		for(Pld v : {s, e}){
			if((v - p).dist() <= EPS) continue;
			if(angleBetween(pnxt-p, pprv-p, v-p) == -1) return false;
		}
	}
	return true;
}

int main() {
	vector<Pll> st(2);
	for(Pll& p : st) cin >> p.x >> p.y;
	vector<int> len(2);
	vector<vector<Pll> > bank(2);
	for(int z = 0; z < 2; z++){
		cin >> len[z];
		bank[z].resize(len[z]);
		for(Pll& p : bank[z]) cin >> p.x >> p.y;
		assert(bank[z].front().y == 0);
		assert(bank[z].back().y == 1000);
	}
	vector<vector<Pld> > ccw_poly(2), pts(2);
	for(int z = 0; z < 2; z++){
		for(Pll p : bank[z]){
			ccw_poly[z].push_back(Pld(p.x, p.y));
			pts[z].push_back(Pld(p.x, p.y));
		}
		pts[z].push_back(toPld(st[z]));
	}
	reverse(ccw_poly[1].begin(), ccw_poly[1].end());
	ccw_poly[0].push_back(Pld(-1, 1000));
	ccw_poly[0].push_back(Pld(-1, 0));
	ccw_poly[1].push_back(Pld(1001, 0));
	ccw_poly[1].push_back(Pld(1001, 1000));
	frac best_dist(ll(1e12), ll(1));
	vector<pair<Pld, Pld> > cands;

	for(int z = 0; z < 2; z++){
		for(int i = 0; i < len[z]; i++){
			for(int j = 0; j + 1 < len[1^z]; j++){
				frac dist = squareSegDist(bank[1^z][j], bank[1^z][j+1], bank[z][i]);
				if(dist < best_dist){
					best_dist = dist;
					cands = {};
				}
				if(dist == best_dist){
					Pld p0 = toPld(bank[z][i]);
					Pld p1 = closestPoint(toPld(bank[1^z][j]), toPld(bank[1^z][j+1]), toPld(bank[z][i]));
					if(z) swap(p0, p1);
					cands.push_back({p0, p1});
				}
			}
		}
	}
	ld bridge_len = (sqrt(best_dist.n) / sqrt(best_dist.d));
	for(int i = 0; i + 1 < len[0]; i++){
		for(int j = 0; j + 1 < len[1]; j++){
			if((bank[0][i] - bank[0][i+1]).cross(bank[1][j] - bank[1][j+1]) != 0) continue;
			bool best = false;
			if(squareSegDist(bank[1][j], bank[1][j+1], bank[0][i]) == best_dist) best = true;
			if(squareSegDist(bank[1][j], bank[1][j+1], bank[0][i+1]) == best_dist) best = true;
			if(squareSegDist(bank[0][i], bank[0][i+1], bank[1][j]) == best_dist) best = true;
			if(squareSegDist(bank[0][i], bank[0][i+1], bank[1][j+1]) == best_dist) best = true;
			if(!best) continue;
			Pld normal = (toPld(bank[0][i+1] - bank[0][i])).normal();
			Pld shift = normal * toPld(bank[1][j] - bank[0][i]).dot(normal);
			for(int a = 0; a < (int)pts[0].size(); a++){
				for(int b = 0; b < (int)pts[1].size(); b++){
					if(a == i || a == i+1 || b == j || b == j+1) continue;
					Pld pst = pts[0][a];
					Pld pen = pts[1][b] - shift;
					vector<Pld> i0 = segInter(pst, pen, toPld(bank[0][i]), toPld(bank[0][i+1]));
					vector<Pld> i1 = segInter(pst + shift, pen + shift, toPld(bank[1][j]), toPld(bank[1][j+1]));
					if(i0.empty() || i1.empty()) continue;
					Pld p0 = i0.front();
					Pld p1 = p0 + shift;
					cands.push_back({p0, p1});
				}
			}
		}
	}
	vector<vector<ld> > dist(2);
	for(int z = 0; z < 2; z++){
		int n = (int)pts[z].size();
		vector<vector<ld> > d(n, vector<ld>(n, 1e15));
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				if(segInsidePolygon(pts[z][i], pts[z][j], ccw_poly[z])){
					d[i][j] = (pts[z][i] - pts[z][j]).dist();
				}
			}
			d[i][i] = 0;
		}
		for(int k = 0; k < n; k++){
			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
				}
			}
		}
		dist[z] = d.back();
	}
	ld road_best = 1e15;
	for(auto [p0, p1] : cands){
		vector<Pld> side_cand = {p0, p1};
		ld road_dist = 0;
		for(int z = 0; z < 2; z++){
			ld road_side = 1e15;
			int n = (int)pts[z].size();
			for(int i = 0; i < n; i++){
				if(segInsidePolygon(side_cand[z], pts[z][i], ccw_poly[z])){
					road_side = min(road_side, (side_cand[z] - pts[z][i]).dist() + dist[z][i]);
				}
			}
			road_dist += road_side;
		}
		road_best = min(road_best, road_dist);
	}
	cout << fixed << setprecision(15);
	cout << bridge_len << ' ' << (bridge_len + road_best) << '\n';
}
