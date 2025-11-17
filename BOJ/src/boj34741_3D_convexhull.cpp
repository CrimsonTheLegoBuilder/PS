#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <queue>
#include <deque>
#include <random>
#include <array>
#include <tuple>
#include <complex>
typedef long long ll;
typedef double ld;
typedef std::vector<ld> Vld;
typedef std::pair<int, int> pi;
const ld INF = 1e17;
const ld TOL = 1e-9;
const ld PI = acos(-1);
const int LEN = 505;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ll sq(const ll& x) { return x * x; }
ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }

int N;
struct P3D {
	int x, y, z;
	bool operator < (const P3D& p) const { return x == p.x ? y == p.y ? z < p.z : y < p.y : x < p.x; }
	bool operator == (const P3D& p) const { return x == p.x && y == p.y && z == p.z; }
};
ll gcd(const P3D& p) {
	ll w = gcd(std::abs(p.x), std::abs(p.y));
	return gcd(std::abs(p.z), w);
}
struct Pos3D {
	ld x, y, z;
	Pos3D(ld x_ = 0, ld y_ = 0, ld z_ = 0) : x(x_), y(y_), z(z_) {}
	bool operator == (const Pos3D& p) const { return zero(x - p.x) && zero(y - p.y) && zero(z - p.z); }
	bool operator != (const Pos3D& p) const { return !zero(x - p.x) || !zero(y - p.y) || !zero(z - p.z); }
	bool operator < (const Pos3D& p) const { return zero(x - p.x) ? zero(y - p.y) ? z < p.z : y < p.y : x < p.x; }
	ld operator * (const Pos3D& p) const { return x * p.x + y * p.y + z * p.z; }
	Pos3D operator / (const Pos3D& p) const { return { y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x }; }
	Pos3D operator + (const Pos3D& p) const { return { x + p.x, y + p.y, z + p.z }; }
	Pos3D operator - (const Pos3D& p) const { return { x - p.x, y - p.y, z - p.z }; }
	Pos3D operator * (const ld& n) const { return { x * n, y * n, z * n }; }
	Pos3D operator / (const ld& n) const { return { x / n, y / n, z / n }; }
	Pos3D& operator += (const Pos3D& p) { x += p.x; y += p.y; z += p.z; return *this; }
	Pos3D& operator -= (const Pos3D& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
	Pos3D& operator *= (const ld& n) { x *= n; y *= n; z *= n; return *this; }
	Pos3D& operator /= (const ld& n) { x /= n; y /= n; z /= n; return *this; }
	ld Euc() const { return x * x + y * y + z * z; }
	ld mag() const { return sqrtl(Euc()); }
	ld lon() const { return atan2(y, x); }
	ld lat() const { return atan2(z, sqrtl(x * x + y * y)); }
	Pos3D operator - () const { return { -x, -y, -z }; }
	Pos3D unit() const { return *this / mag(); }
	Pos3D norm(const Pos3D& p) const { return (*this / p).unit(); }
	Pos3D rodrigues_rotate(const ld& th, const Pos3D& axis) const {
		ld s = sin(th), c = cos(th); Pos3D n = axis.unit();
		return n * (*this * n) * (1 - c) + (*this * c) + (n / *this) * s;
	}
	friend ld rad(const Pos3D& p1, const Pos3D& p2) { return atan2l((p1 / p2).mag(), p1 * p2); }
	friend std::istream& operator >> (std::istream& is, Pos3D& p) { is >> p.x >> p.y >> p.z; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos3D& p) { os << p.x << " " << p.y << " " << p.z; return os; }
};
typedef std::vector<Pos3D> Polyhedron;
const Pos3D O3D = { 0, 0, 0 };
const Pos3D X_axis = { 1, 0, 0 };
const Pos3D Y_axis = { 0, 1, 0 };
const Pos3D Z_axis = { 0, 0, 1 };
const Pos3D MAXP3D = { INF, INF, INF };
std::vector<Pos3D> pos;
Pos3D point(const Pos3D Xaxis, const Pos3D Yaxis, const ld& th) { return Xaxis * cos(th) + Yaxis * sin(th); }
ld angle(const Pos3D Xaxis, const Pos3D Yaxis, const Pos3D& p) { return norm(atan2(Yaxis * p, Xaxis * p)); }
Pos3D cross(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) / (d3 - d2); }
ld dot(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return (d2 - d1) * (d3 - d2); }
int ccw(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3, const Pos3D& norm) { return sign(cross(d1, d2, d3) * norm); }
ld area(const std::vector<Pos3D>& H, const Pos3D& norm) {
	if (H.size() < 3) return 0;
	ld ret = 0;
	int sz = H.size();
	for (int i = 0; i < sz; i++) {
		const Pos3D& cur = H[i], nxt = H[(i + 1) % sz];
		ret += cross(H[0], cur, nxt) * norm / norm.mag();
	}
	return std::abs(ret * .5);
}
bool on_seg_strong(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).mag()) && sign(dot(d1, d3, d2)) >= 0; }
bool on_seg_weak(const Pos3D& d1, const Pos3D& d2, const Pos3D& d3) { return zero(cross(d1, d2, d3).mag()) && sign(dot(d1, d3, d2)) > 0; }
bool collinear(const Pos3D& a, const Pos3D& b, const Pos3D& c) { return zero(((b - a) / (c - b)).Euc()); }
bool coplanar(const Pos3D& a, const Pos3D& b, const Pos3D& c, const Pos3D& p) { return zero(cross(a, b, c) * (p - a)); }
bool above(const Pos3D& a, const Pos3D& b, const Pos3D& c, const Pos3D& p) { return cross(a, b, c) * (p - a) > 0; }
int prep(std::vector<Pos3D>& p) {//refer to Koosaga'
	shuffle(p.begin(), p.end(), std::mt19937(0x14004));
	int dim = 1;
	for (int i = 1; i < p.size(); i++) {
		if (dim == 1) {
			if (p[0] != p[i])
				std::swap(p[1], p[i]), ++dim;
		}
		else if (dim == 2) {
			if (!collinear(p[0], p[1], p[i]))
				std::swap(p[2], p[i]), ++dim;
		}
		else if (dim == 3) {
			if (!coplanar(p[0], p[1], p[2], p[i]))
				std::swap(p[3], p[i]), ++dim;
		}
	}
	//assert(dim == 4);
	return dim;
}
struct Face {
	int v[3];
	Face(int a = 0, int b = 0, int c = 0) { v[0] = a; v[1] = b; v[2] = c; }
	Pos3D norm(const std::vector<Pos3D>& C) const { return cross(C[v[0]], C[v[1]], C[v[2]]); }
	bool visible(const std::vector<Pos3D>& P, int i) const { return (P[i] - P[v[0]]) * norm(P) > 0; }
	bool outer(const std::vector<Pos3D>& P, const Pos3D& q) const { return sign((q - P[v[0]]) * norm(P)) >= 0; }
};
ld dist(const std::vector<Pos3D>& C, const Face& F, const Pos3D& p) {
	Pos3D norm = cross(C[F.v[0]], C[F.v[1]], C[F.v[2]]);
	ld ret = norm * (p - C[F.v[0]]);
	return std::abs(ret / (norm.mag()));
}
ld get_min_dist(const std::vector<Pos3D>& C, const std::vector<Face>& F, const Pos3D& p) {
	ld MIN = INF;
	for (const Face& face : F) MIN = std::min(MIN, dist(C, face, p));
	return MIN;
}
struct Edge {
	int face_num, edge_num;
	Edge(int t = 0, int v = 0) : face_num(t), edge_num(v) {}
};
bool col = 0, cop = 0;
bool ff = 0;
std::vector<Face> convex_hull_3D(std::vector<Pos3D>& candi) {//incremental construction
	// 3D Convex Hull in O(n log n)
	// Very well tested. Good as long as not all points are coplanar
	// In case of collinear faces, returns arbitrary triangulation
	// Credit: Benq
	// refer to Koosaga'
	col = 0, cop = 0;
	int suf = prep(candi);
	if (suf <= 2) { col = 1; return {}; };
	if (suf == 3) { cop = 1; return {}; };
	int sz = candi.size();
	if (ff) for (Pos3D& p : candi) p = p.unit();
	std::vector<Face> faces;
	std::vector<int> active;//whether face is active - face faces outside 
	std::vector<std::vector<int>> vis(sz);//faces visible from each point
	std::vector<std::vector<int>> rvis;//points visible from each face
	std::vector<std::array<Edge, 3>> other;//other face adjacent to each edge of face
	auto ad = [&](const int& a, const int& b, const int& c) -> void {//add face
		faces.push_back({ a, b, c });
		active.push_back(1);
		rvis.emplace_back();
		other.emplace_back();
		return;
		};
	auto visible = [&](const int& a, const int& b) -> void {
		vis[b].push_back(a);
		rvis[a].push_back(b);
		return;
		};
	auto abv = [&](const int& a, const int& b) -> bool {//above
		Face tri = faces[a];
		return above(candi[tri.v[0]], candi[tri.v[1]], candi[tri.v[2]], candi[b]);
		};
	auto edge = [&](const Edge& e) -> pi {
		return { faces[e.face_num].v[e.edge_num], faces[e.face_num].v[(e.edge_num + 1) % 3] };
		};
	auto glue = [&](const Edge& a, const Edge& b) -> void {//link two faces by an edge
		pi x = edge(a); assert(edge(b) == pi(x.second, x.first));
		other[a.face_num][a.edge_num] = b;
		other[b.face_num][b.edge_num] = a;
		return;
		};//ensure face 0 is removed when i = 3

	ad(0, 1, 2), ad(0, 2, 1);
	if (abv(1, 3)) std::swap(candi[1], candi[2]);
	for (int i = 0; i < 3; i++) glue({ 0, i }, { 1, 2 - i });
	for (int i = 3; i < sz; i++) visible(abv(1, i), i);//coplanar points go in rvis[0]

	std::vector<int> label(sz, -1);
	for (int i = 3; i < sz; i++) {//incremental construction
		std::vector<int> rem;
		for (auto& v : vis[i]) if (active[v]) { active[v] = 0, rem.push_back(v); }
		if (!rem.size()) continue;//hull unchanged

		int st = -1;//start idx
		for (const int& v : rem) {
			for (int j = 0; j < 3; j++) {
				int o = other[v][j].face_num;
				if (active[o]) {//create new face!
					int idx1, idx2;
					std::tie(idx1, idx2) = edge({ v, j });
					ad(idx1, idx2, i);
					st = idx1;
					int cur = rvis.size() - 1;
					label[idx1] = cur;

					std::vector<int> tmp;
					set_union(rvis[v].begin(), rvis[v].end(), rvis[o].begin(), rvis[o].end(), back_inserter(tmp));
					//merge sorted vectors ignoring duplicates

					for (auto& x : tmp) if (abv(cur, x)) visible(cur, x);
					//if no rounding errors then guaranteed that only x > i matters

					glue({ cur, 0 }, other[v][j]);//glue old, new face
				}
			}
		}
		for (int x = st, y; ; x = y) {//glue new faces together
			int X = label[x];
			glue({ X, 1 }, { label[y = faces[X].v[1]], 2 });
			if (y == st) break;
		}
	}
	std::vector<Face> hull3D;
	for (int i = 0; i < faces.size(); i++) if (active[i]) hull3D.push_back(faces[i]);
	return hull3D;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(13);
	std::cin >> N;
	if (N < 4) { std::cout << "90\n"; return; }
	std::vector<P3D> C(N);
	for (int i = 0; i < N; i++) {
		std::cin >> C[i].x >> C[i].y >> C[i].z;
		ll g = gcd(C[i]); C[i].x /= g; C[i].y /= g; C[i].z /= g;
	}
	std::sort(C.begin(), C.end());
	C.erase(unique(C.begin(), C.end()), C.end());
	int sz = C.size();
	Polyhedron P(sz);
	for (int i = 0; i < sz; i++) {
		P[i] = Pos3D(C[i].x, C[i].y, C[i].z);
	}
	std::vector<Face> H = convex_hull_3D(P);
	if (col || cop) { std::cout << "90\n"; return; }
	bool fuck = 0;
	for (const Face& f : H) {
		if (f.outer(P, O3D)) { fuck = 1; break; }
	}
	if (fuck) { std::cout << "90\n"; return; }
	for (int i = 0; i < sz; i++) {
		P[i] = Pos3D(C[i].x, C[i].y, C[i].z);
		//P[i] = P[i].unit();
	}
	ff = 1;
	H = convex_hull_3D(P);
	ld r = 0;
	for (const Face& f : H) {
		Pos3D n = f.norm(P).unit();
		//assert(eq(rad(n, P[f.v[0]]), rad(n, P[f.v[1]])));
		//assert(eq(rad(n, P[f.v[1]]), rad(n, P[f.v[2]])));
		//assert(eq(rad(n, P[f.v[2]]), rad(n, P[f.v[0]])));
		ld t = rad(n, P[f.v[0]]);
		r = std::max(r, t);
	}
	std::cout << r * 180 / PI << "\n";
	return;
}
int main() { solve(); return 0; }//boj34741