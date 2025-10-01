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
#include <numeric>
typedef long long ll;
typedef long double ld;
//typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ll> Vll;
typedef std::vector<ld> Vld;
typedef std::vector<bool> Vbool;
const ld INF = 1e18;
const ll LINF = 1e17;
const ld TOL = 1e-9;
const ld PI = acos(-1);
const int LEN = 9e4;
inline ll sqr(const int& a) { return (ll)a * a; }
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ll sq(const ll& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo = 0, const ld& hi = 1) { return std::min(hi, std::max(lo, x)); }
inline ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }
//inline ll gcd(ll x, ll y, ll z) {
//	x = std::abs(x); y = std::abs(y); z = std::abs(z);
//	ll w = gcd(x, y);
//	return gcd(w, z);
//}

//Struct Pos { int x, y, i }
//Struct Pdd - 보로노이 다이어그램 제작에만 쓰임
//Class or Struct Delaunay trianlges :: 재홍이가 접근할 일 없도록 구성
//Class or Struct KDtree - Tree, Pos :: 재홍이가 접근할 일 없도록 구성
//입력 및 전처리는 전부 기하모듈에서 일어나야하는데 후에 하는 접근은 전부 다른 자료구조들에서 행해짐
//전처리를 전부 끝내놓고 나면 기하모듈에 접근할 일 없도록 해야함
//Class or Struct Query - int t, Pos s, Pos e, int u, int p :: 타입, 점 2개, 부모별 번호. 순례끝점 위치는 KDtree 조회 후 i 에 번호 기록

int N, M, T, Q, Q1, Q2;
Vint DT[LEN];//delaunay triangle
std::string name[LEN];

/* GEOMETRY STRUCT */
struct Pii {
	ll x, y; int i;
	Pii(ll x_ = 0, ll y_ = 0, int i_ = 0) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pii& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pii& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pii& p) const { return x == p.x ? y < p.y : x < p.x; }
	bool operator <= (const Pii& p) const { return x == p.x ? y <= p.y : x <= p.x; }
	Pii operator + (const Pii& p) const { return { x + p.x, y + p.y }; }
	Pii operator - (const Pii& p) const { return { x - p.x, y - p.y }; }
	Pii operator * (const ll& n) const { return { x * n, y * n }; }
	Pii operator / (const ll& n) const { return { x / n, y / n }; }
	ll operator * (const Pii& p) const { return { x * p.x + y * p.y }; }
	ll operator / (const Pii& p) const { return { x * p.y - y * p.x }; }
	Pii& operator += (const Pii& p) { x += p.x; y += p.y; return *this; }
	Pii& operator -= (const Pii& p) { x -= p.x; y -= p.y; return *this; }
	Pii& operator *= (const ll& n) { x *= n; y *= n; return *this; }
	Pii& operator /= (const ll& n) { x /= n; y /= n; return *this; }
	Pii operator - () const { return { -x, -y }; }
	Pii operator ~ () const { return { -y, x }; }
	Pii operator ! () const { return { y, x }; }
	ll xy() const { return x * y; }
	ll Euc() const { return x * x + y * y; }
	ll Man() const { return std::abs(x) + std::abs(y); }
	ld mag() const { return hypot(x, y); }
	friend std::istream& operator >> (std::istream& is, Pii& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pii& p) { os << p.x << " " << p.y; return os; }
};
const Pii Oii = { 0, 0 };
const Pii INF_PT = { (int)1e9, (int)1e9 };
typedef std::vector<Pii> Vpii;
ll cross(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pii& d1, const Pii& d2, const Pii& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return sign(cross(d1, d2, d3, d4)); }
/* GEOMETRY STRUCT */

/* DELAUNAY TRIANGULATION */
/*

tested in range -1e6 < x, y < 1e6;
Delaunator - https://github.com/abellgithub/delaunator-cpp/blob/master/include/delaunator.cpp

*/
struct Pos {
	ld x, y; int i;
	Pos(ld x_ = 0, ld y_ = 0, int i_ = 0) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pos& p) const { return zero(x - p.x) && zero(y - p.y); }
	bool operator != (const Pos& p) const { return !zero(x - p.x) || !zero(y - p.y); }
	bool operator < (const Pos& p) const { return zero(x - p.x) ? y < p.y : x < p.x; }
	Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
	Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
	Pos operator * (const ld& n) const { return { x * n, y * n }; }
	Pos operator / (const ld& n) const { return { x / n, y / n }; }
	ld operator * (const Pos& p) const { return { x * p.x + y * p.y }; }
	ld operator / (const Pos& p) const { return { x * p.y - y * p.x }; }
	ld operator ^ (const Pos& p) const { return { x * p.y - y * p.x }; }
	Pos operator - () const { return Pos(-x, -y); }
	Pos operator ~ () const { return { -y, x }; }
	Pos operator ! () const { return { -x, -y }; }
	Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
	Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
	Pos& operator *= (const ld& n) { x *= n; y *= n; return *this; }
	Pos& operator /= (const ld& n) { x /= n; y /= n; return *this; }
	ld xy() const { return x * y; }
	Pos rot(const ld& t) { return { x * cos(t) - y * sin(t), x * sin(t) + y * cos(t) }; }
	ld Euc() const { return x * x + y * y; }
	ld mag() const { return sqrt(Euc()); }
	Pos unit() const { return *this / mag(); }
	ld rad() const { return atan2(y, x); }
	friend ld rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
	int quad() const { return sign(y) == 1 || (sign(y) == 0 && sign(x) >= 0); }
	friend bool cmpq(const Pos& a, const Pos& b) { return (a.quad() != b.quad()) ? a.quad() < b.quad() : sign(a / b) > 0; }
	bool close(const Pos& p) const { return zero((*this - p).Euc()); }
	bool close(const Pos& rhs,
		const ld span = 1.,
		const ld tol = 1e-20) const {
		return ((*this - rhs).Euc() / span) < tol;
	}
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = { 0, 0 };
typedef std::vector<Pos> Polygon;
Pos conv(const Pii& p) { return Pos(p.x, p.y, p.i); }
Pii conv(const Pos& p) { ll x = round(p.x); ll y = round(p.y); return Pii(x, y); }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
ld cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return sign(cross(d1, d2, d3, d4)); }
bool counterclockwise(const Pos& p, const Pos& q, const Pos& r) { return ccw(p, q, r) == 1; }
bool clockwise(const Pos& p, const Pos& q, const Pos& r) { return ccw(p, q, r) == -1; }
bool collinear(const Pos& p, const Pos& q, const Pos& r) { return !ccw(p, q, r); }
bool collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return !ccw(d1, d2, d3) && !ccw(d1, d2, d4); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
ld dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) >= 0; }
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) { return !ccw(d1, d2, d3) && sign(dot(d1, d3, d2)) > 0; }
ld projection(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d1) / (d2 - d1).mag(); }
ld dist(const Pos& d1, const Pos& d2, const Pos& t) { return cross(d1, d2, t) / (d1 - d2).mag(); }
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) { ld a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2); return (p1 * a2 + p2 * a1) / (a1 + a2); }
ld circumradius(const Pos& p1, const Pos& p2, const Pos& p3) {
	Pos d = p2 - p1;
	Pos e = p3 - p1;
	const ld bl = d.Euc();
	const ld cl = e.Euc();
	const ld det = d / e;
	Pos radius((e.y * bl - d.y * cl) * 0.5 / det,
		(d.x * cl - e.x * bl) * 0.5 / det);
	if ((bl > 0.0 || bl < 0.0) && (cl > 0.0 || cl < 0.0) &&
		(det > 0.0 || det < 0.0))
		return radius.Euc();
	return (std::numeric_limits<double>::max)();
}
Pos circumcenter(const Pos& p1, const Pos& p2, const Pos& p3) {
	Pos d = p2 - p1;
	Pos e = p3 - p1;
	const ld bl = d.Euc();
	const ld cl = e.Euc();
	const ld det = d / e;
	Pos radius((e.y * bl - d.y * cl) * 0.5 / det,
		(d.x * cl - e.x * bl) * 0.5 / det);
	return p1 + radius;
}
//struct Circle {
//	Pos c;
//	ld r;
//	Circle(Pos C = Pos(0, 0), ld R = 0) : c(C), r(R) {}
//	bool operator == (const Circle& C) const { return c == C.c && std::abs(r - C.r) < TOL; }
//	bool operator != (const Circle& C) const { return !(*this == C); }
//	bool operator < (const Circle& q) const {
//		ld dist = (c - q.c).mag();
//		return r < q.r && dist + r < q.r + TOL;
//	}
//	bool operator > (const Pos& p) const { return r > (c - p).mag(); }
//	bool operator >= (const Pos& p) const { return r + TOL > (c - p).mag(); }
//	bool operator < (const Pos& p) const { return r < (c - p).mag(); }
//	Circle operator + (const Circle& C) const { return { c + C.c, r + C.r }; }
//	Circle operator - (const Circle& C) const { return { c - C.c, r - C.r }; }
//	ld H(const ld& th) const { return sin(th) * c.x + cos(th) * c.y + r; }//coord trans | check right
//	ld A() const { return 1. * r * r * PI; }
//	friend std::istream& operator >> (std::istream& is, Circle& c) { is >> c.c >> c.r; return is; }
//	friend std::ostream& operator << (std::ostream& os, const Circle& c) { os << c.c << " " << c.r; return os; }
//} INVAL = { { 0, 0 }, (std::numeric_limits<ld>::max)() };
//Circle enclose_circle(const Pos& u, const Pos& v, const Pos& w) {
//	if (!ccw(u, v, w)) return INVAL;
//	Pos m1 = (u + v) * .5, v1 = ~(v - u);
//	Pos m2 = (u + w) * .5, v2 = ~(w - u);
//	Pos c = intersection(m1, m1 + v1, m2, m2 + v2);
//	return Circle(c, (u - c).mag());
//}
//Circle circumcircle(const Pos& p1, const Pos& p2, const Pos& p3) {
//	Pos d = p2 - p1;
//	Pos e = p3 - p1;
//	const ld bl = d.Euc();
//	const ld cl = e.Euc();
//	const ld det = d / e;
//	Pos radius((e.y * bl - d.y * cl) * 0.5 / det,
//		(d.x * cl - e.x * bl) * 0.5 / det);
//	Pos c;
//	ld r;
//	if ((bl > 0.0 || bl < 0.0) && (cl > 0.0 || cl < 0.0) &&
//		(det > 0.0 || det < 0.0))
//		c = p1 + radius, r = radius.Euc();
//	else c = Pos(INF, INF), r = std::numeric_limits<ld>::max();
//	return Circle(c, r);
//}
bool in_circle(const Pos& a, const Pos& b, const Pos& c, const Pos& p) {
	const Pos d = a - p;
	const Pos e = b - p;
	const Pos f = c - p;

	const double ap = d.Euc();
	const double bp = e.Euc();
	const double cp = f.Euc();

	return d / (e * cp - f * bp) + ap * (e / f) < 0.0;
}
class BBox2 {//jinhwanlazy
	constexpr static auto INF = std::numeric_limits<double>::max();
private:
	Pos bottom_left_;
	Pos top_right_;
	Pos center_;
	double span_;

public:
	BBox2(const std::vector<Pos>& points) {
		top_right_ = Pos(-INF, -INF);
		bottom_left_ = Pos(INF, INF);
		for (const Pos& p : points) {
			bottom_left_.x = std::min(bottom_left_.x, p.x);
			bottom_left_.y = std::min(bottom_left_.y, p.y);
			top_right_.x = std::max(top_right_.x, p.x);
			top_right_.y = std::max(top_right_.y, p.y);
		}
		center_ = (bottom_left_ + top_right_) / 2;
		span_ = (bottom_left_ - top_right_).Euc();
	}

	const Pos& bottomLeft() const { return bottom_left_; }
	const Pos& topRight() const { return top_right_; }
	const Pos& center() const { return center_; }
	const double& span() const { return span_; }
};
class Delaunator {
public:
	constexpr static auto INF = std::numeric_limits<double>::max();

	std::vector<Pos> points_;
	std::vector<std::size_t> triangles_;
	std::vector<std::size_t> halfedges_;
	std::vector<std::size_t> hull_prev_;
	std::vector<std::size_t> hull_next_;
	std::vector<std::size_t> hull_tri_;
	std::size_t hull_start_;

private:
	static constexpr std::size_t INVALID_INDEX = -1;

	std::vector<std::size_t> hull_hash_;
	Pos center_;
	std::size_t hash_size_;

public:
	Delaunator(std::vector<Pos> const& points) : points_(points) {
		std::size_t n = points.size();

		BBox2 bbox(points_);
		Pos center = bbox.center();

		std::size_t i0 = INVALID_INDEX;
		std::size_t i1 = INVALID_INDEX;
		std::size_t i2 = INVALID_INDEX;

		double min_dist = INF;
		for (size_t i = 0; i < points_.size(); ++i)
		{
			const Pos& p = points_[i];
			const double d = (p - center).Euc();
			if (d < min_dist) {
				i0 = i;
				min_dist = d;
			}
		}
		Pos p0 = points_[i0];

		min_dist = (std::numeric_limits<double>::max)();
		for (std::size_t i = 0; i < n; i++) {
			if (i == i0) continue;
			const double d = (p0 - points_[i]).Euc();
			if (d < min_dist && d > 0.0) {
				i1 = i;
				min_dist = d;
			}
		}
		Pos p1 = points_[i1];

		double min_radius = INF;
		for (std::size_t i = 0; i < n; i++) {
			if (i == i0 || i == i1)
				continue;
			const double r = circumradius(p0, p1, points_[i]);
			if (r < min_radius) {
				i2 = i;
				min_radius = r;
			}
		}
		Pos p2 = points_[i2];

		if (!(min_radius < INF)) {
			throw std::runtime_error("not triangulation");
		}

		if (counterclockwise(p0, p1, p2)) {
			std::swap(i1, i2);
			std::swap(p1, p2);
		}

		center_ = circumcenter(p0, p1, p2);

		std::vector<double> dists;
		dists.reserve(points_.size());
		for (const auto& p : points_)
			dists.push_back((p - center_).Euc());

		// sort the points by distance from the seed triangle circumcenter
		std::vector<std::size_t> ids(n);
		std::iota(ids.begin(), ids.end(), 0);
		std::sort(ids.begin(), ids.end(),
			[&dists](std::size_t i, std::size_t j) { return dists[i] < dists[j]; });

		// initialize a hash table for storing edges of the advancing convex hull
		hash_size_ = static_cast<std::size_t>(std::ceil(std::sqrt(n)));
		hull_hash_.resize(hash_size_, INVALID_INDEX);

		// initialize arrays for tracking the edges of the advancing convex hull
		hull_prev_.resize(n);
		hull_next_.resize(n);
		hull_tri_.resize(n);

		hull_start_ = i0;

		size_t hull_size = 3;

		hull_next_[i0] = hull_prev_[i2] = i1;
		hull_next_[i1] = hull_prev_[i0] = i2;
		hull_next_[i2] = hull_prev_[i1] = i0;

		hull_tri_[i0] = 0;
		hull_tri_[i1] = 1;
		hull_tri_[i2] = 2;

		hull_hash_[hash_key(p0)] = i0;
		hull_hash_[hash_key(p1)] = i1;
		hull_hash_[hash_key(p2)] = i2;

		std::size_t max_triangles_ = n < 3 ? 1 : 2 * n - 5;
		triangles_.reserve(max_triangles_ * 3);
		halfedges_.reserve(max_triangles_ * 3);
		add_triangle(i0, i1, i2, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);
		Pos p_prev{ std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };

		// Go through points based on distance from the center.
		for (std::size_t k = 0; k < n; k++) {
			const std::size_t i = ids[k];
			const Pos p = points_[i];

			// skip near-duplicate points
			if (k > 0 && p == p_prev)
				continue;
			p_prev = p;

			if (p == p0 || p == p1 || p == p2) {
				continue;
			}

			// find a visible edge on the convex hull using edge hash
			std::size_t start = 0;

			size_t key = hash_key(p);
			for (size_t j = 0; j < hash_size_; j++) {
				start = hull_hash_[(key + j) % hash_size_];
				if (start != INVALID_INDEX && start != hull_next_[start])
					break;
			}

			assert(hull_prev_[start] != start);
			assert(hull_prev_[start] != INVALID_INDEX);

			start = hull_prev_[start];
			size_t e = start;
			size_t q;

			// Advance until we find a place in the hull where our current point can be added.
			while (true) {
				q = hull_next_[e];
				if (p.close(points_[e], bbox.span()) || p.close(points_[q], bbox.span())) {
					e = INVALID_INDEX;
					break;
				}
				if (counterclockwise(p, points_[e], points_[q]))
					break;
				e = q;
				if (e == start) {
					e = INVALID_INDEX;
					break;
				}
			}

			if (e == INVALID_INDEX)  // likely a near-duplicate point; skip it
				continue;

			// add the first triangle from the point
			std::size_t t = add_triangle(e, i, hull_next_[e], INVALID_INDEX, INVALID_INDEX, hull_tri_[e]);

			hull_tri_[i] = legalize(t + 2);  // Legalize the triangle we just added.
			hull_tri_[e] = t;
			hull_size++;

			// walk forward through the hull, adding more triangles_ and flipping recursively
			std::size_t next = hull_next_[e];
			while (true) {
				q = hull_next_[next];
				if (!counterclockwise(p, points_[next], points_[q]))
					break;
				t = add_triangle(next, i, q, hull_tri_[i], INVALID_INDEX, hull_tri_[next]);
				hull_tri_[i] = legalize(t + 2);
				hull_next_[next] = next;  // mark as removed
				hull_size--;
				next = q;
			}

			// walk backward from the other side, adding more triangles_ and flipping
			if (e == start) {
				while (true) {
					q = hull_prev_[e];
					if (!counterclockwise(p, points_[q], points_[e]))
						break;
					t = add_triangle(q, i, e, INVALID_INDEX, hull_tri_[e], hull_tri_[q]);
					legalize(t + 2);
					hull_tri_[q] = t;
					hull_next_[e] = e;  // mark as removed
					hull_size--;
					e = q;
				}
			}

			// update the hull indices
			hull_prev_[i] = e;
			hull_start_ = e;
			hull_prev_[next] = i;
			hull_next_[e] = i;
			hull_next_[i] = next;

			hull_hash_[hash_key(p)] = i;
			hull_hash_[hash_key(points_[e])] = e;
		}
	}

private:
	std::size_t legalize(std::size_t a) {
		std::size_t i = 0;
		std::size_t ar = 0;
		std::vector<std::size_t> edges_stack;

		// recursion eliminated with a fixed-size stack
		while (true) {
			const size_t b = halfedges_[a];

			/* if the pair of triangles_ doesn't satisfy the Delaunay condition
			 * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
			 * then do the same check/flip recursively for the new pair of triangles_
			 *
			 *           pl                    pl
			 *          /||\                  /  \
			 *       al/ || \bl            al/    \a
			 *        /  ||  \              /      \
			 *       /  a||b  \    flip    /___ar___\
			 *     p0\   ||   /p1   =>   p0\---bl---/p1
			 *        \  ||  /              \      /
			 *       ar\ || /br             b\    /br
			 *          \||/                  \  /
			 *           pr                    pr
			 */
			const size_t a0 = 3 * (a / 3);
			ar = a0 + (a + 2) % 3;

			if (b == INVALID_INDEX) {
				if (i > 0) {
					i--;
					a = edges_stack[i];
					continue;
				}
				else {
					// i = INVALID_INDEX;
					break;
				}
			}

			const size_t b0 = 3 * (b / 3);
			const size_t al = a0 + (a + 1) % 3;
			const size_t bl = b0 + (b + 2) % 3;

			const std::size_t p0 = triangles_[ar];
			const std::size_t pr = triangles_[a];
			const std::size_t pl = triangles_[al];
			const std::size_t p1 = triangles_[bl];

			const bool illegal = in_circle(points_[p0], points_[pr], points_[pl], points_[p1]);

			if (illegal) {
				triangles_[a] = p1;
				triangles_[b] = p0;

				auto hbl = halfedges_[bl];

				// Edge swapped on the other side of the hull (rare).
				// Fix the halfedge reference
				if (hbl == INVALID_INDEX) {
					std::size_t e = hull_start_;
					do {
						if (hull_tri_[e] == bl) {
							hull_tri_[e] = a;
							break;
						}
						e = hull_prev_[e];
					} while (e != hull_start_);
				}
				link(a, hbl);
				link(b, halfedges_[ar]);
				link(ar, bl);
				std::size_t br = b0 + (b + 1) % 3;

				if (i < edges_stack.size()) {
					edges_stack[i] = br;
				}
				else {
					edges_stack.push_back(br);
				}
				i++;

			}
			else {
				if (i > 0) {
					i--;
					a = edges_stack[i];
					continue;
				}
				else {
					break;
				}
			}
		}
		return ar;
	};

	// monotonically increases with real angle, but doesn't need expensive trigonometry
	static inline double pseudo_angle(const double dx, const double dy) {
		const double p = dx / (std::abs(dx) + std::abs(dy));
		return (dy > 0.0 ? 3.0 - p : 1.0 + p) / 4.0;  // [0..1)
	}

	std::size_t hash_key(double x, double y) const {
		const double dx = x - center_.x;
		const double dy = y - center_.y;
		size_t key = std::llround(std::floor(pseudo_angle(dx, dy) * static_cast<double>(hash_size_)));
		return key % hash_size_;
	};

	std::size_t hash_key(const Pos& p) const {
		const Pos v = p - center_;
		size_t key = std::llround(std::floor(pseudo_angle(v.x, v.y) * static_cast<double>(hash_size_)));
		return key % hash_size_;
	};

	std::size_t add_triangle(std::size_t i0,
		std::size_t i1,
		std::size_t i2,
		std::size_t a,
		std::size_t b,
		std::size_t c) {
		std::size_t t = triangles_.size();
		triangles_.push_back(i0);
		triangles_.push_back(i1);
		triangles_.push_back(i2);
		link(t, a);
		link(t + 1, b);
		link(t + 2, c);
		return t;
	}

	void link(std::size_t a, std::size_t b) {
		std::size_t s = halfedges_.size();
		if (a == s) {
			halfedges_.push_back(b);
		}
		else if (a < s) {
			halfedges_[a] = b;
		}
		else {
			throw std::runtime_error("Cannot link edge");
		}
		if (b != INVALID_INDEX) {
			std::size_t s2 = halfedges_.size();
			if (b == s2) {
				halfedges_.push_back(a);
			}
			else if (b < s2) {
				halfedges_[b] = a;
			}
			else {
				throw std::runtime_error("Cannot link edge");
			}
		}
	}
};
/* DELAUNAY TRIANGULATION */

/* KD-TREE */
Pii pos[LEN];
struct KDNode {
	Pii p;//mid point
	bool spl;//dx < dy ?
	int sx, ex, sy, ey;
	KDNode(Pii p_ = Pii(0, 0), bool spl_ = 0, int sx_ = 0, int ex_ = 0, int sy_ = 0, int ey_ = 0) :
		p(p_), spl(spl_), sx(sx_), ex(ex_), sy(sy_), ey(ey_) {
	}
} kdtree[LEN << 2];//segment tree
bool V[LEN << 2];//visited
bool cmpx(const Pii& p, const Pii& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pii& p, const Pii& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
bool cmpi(const Pii& p, const Pii& q) { return p.i < q.i; }
void init(int s = 0, int e = N - 1, int n = 1) {//divide & conquer
	ll MAXX = -INF, MAXY = -INF, MINX = INF, MINY = INF;
	int m = s + e >> 1;
	for (int i = s; i <= e; i++) {
		MAXX = std::max(MAXX, pos[i].x);
		MAXY = std::max(MAXY, pos[i].y);
		MINX = std::min(MINX, pos[i].x);
		MINY = std::min(MINY, pos[i].y);
	}
	kdtree[n].spl = (MAXX - MINX) < (MAXY - MINY);//slope cmp :: dx < dy
	if (kdtree[n].spl) std::sort(pos + s, pos + e + 1, cmpy);//if dy is dominant, sort by y
	else std::sort(pos + s, pos + e + 1, cmpx);//if dx is dominant, sort by x
	V[n] = 1;
	kdtree[n] = KDNode(pos[m], kdtree[n].spl, MINX, MAXX, MINY, MAXY);
	if (s <= m - 1) init(s, m - 1, n << 1);
	if (m + 1 <= e) init(m + 1, e, n << 1 | 1);
	return;
}
#define IDX x
#define DIST y
Pii search(const Pii& q, ll X = LINF, int n = 1) {//divide & conquer | refer to koosaga
	ll S = LINF; Pii D = Pii(-1, LINF);
	int i = -1;
	if (q != kdtree[n].p) {
		ll A = (q - kdtree[n].p).Euc();
		if (A < X) { i = kdtree[n].p.i, S = A; }
		else S = X;
	}
	if (kdtree[n].spl) {//if dy is dominant
		if (!cmpy(kdtree[n].p, q)) {//p.y >= t.y
			if (V[n << 1]) D = search(q, S, n << 1);//L search
			if (D.DIST < S) { i = D.IDX, S = D.DIST; }
			if (V[n << 1 | 1] && sqr(kdtree[n << 1 | 1].sy - q.y) < S) D = search(q, S, n << 1 | 1);
			//if R < MIN, execute R search
		}
		else {//p.y < t.y
			if (V[n << 1 | 1]) D = search(q, S, n << 1 | 1);//R search
			if (D.DIST < S) { i = D.IDX, S = D.DIST; }
			if (V[n << 1] && sqr(kdtree[n << 1].ey - q.y) < S) D = search(q, S, n << 1);
			//if L < MIN, execute L search
		}
	}
	else {//if dx is dominant
		if (!cmpx(kdtree[n].p, q)) {//p.x >= t.x
			if (V[n << 1]) D = search(q, S, n << 1);
			if (D.DIST < S) { i = D.IDX, S = D.DIST; }
			if (V[n << 1 | 1] && sqr(kdtree[n << 1 | 1].sx - q.x) < S) D = search(q, S, n << 1 | 1);
		}
		else {//p.x < t.x
			if (V[n << 1 | 1]) D = search(q, S, n << 1 | 1);
			if (D.DIST < S) { i = D.IDX, S = D.DIST; }
			if (V[n << 1] && sqr(kdtree[n << 1].ex - q.x) < S) D = search(q, S, n << 1);
		}
	}
	if (D.DIST < S) { i = D.IDX, S = D.DIST; }
	return Pii(i, S);
}
/* KD-TREE */

/* QUERY */
struct Query {
	int t;
	//bool t;
	Pii s, e;
	int u, p;
	Query(int t_ = 0, Pii s_ = Pii(), Pii e_ = Pii(), int u_ = 0, int p_ = 0) : 
		t(t_), s(s_), e(e_), u(u_), p(p_) {}
	friend std::istream& operator >> (std::istream& is, Query& q) {
		is >> q.t;
		if (q.t == 1) { std::cin >> q.s >> q.e; }
		else { std::cin >> q.u >> q.p; }
		return is;
	}
} qry[LEN << 1];
/* QUERY */

/* MAIN */
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(13);

	/* INPUT */
	std::cin >> N >> Q1 >> Q2; Q = Q1 + Q2;
	Vpii star(N);
	for (int i = 0; i < N; i++) {
		std::cin >> star[i] >> name[i];
		star[i].i = i;
		pos[i] = star[i];
	}
	/* INPUT */

	/* KD TREE INIT */
	memset(V, 0, sizeof V);
	init();
	/* KD TREE INIT */

	/* DELAUNAY & PRIM INIT */
	Polygon C;
	for (Pii& p : star) C.push_back(conv(p));
	Delaunator DTR(C);
	for (int i = 0; i < DTR.triangles_.size(); i += 3) {
		const int& a = DTR.points_[DTR.triangles_[i]].i;
		const int& b = DTR.points_[DTR.triangles_[i + 1]].i;
		const int& c = DTR.points_[DTR.triangles_[i + 2]].i;
		DT[a].push_back(b); DT[a].push_back(c);
		DT[b].push_back(a); DT[b].push_back(c);
		DT[c].push_back(a); DT[c].push_back(b);
	}
	for (int i = 0; i < N; i++) {
		for (const int& j : DT[i]) {
			//compare name - LCP
			//name[i] ? name[j]
			//connect node i - j
		}
	}
	//mst - PRIM
	/* DELAUNAY & PRIM INIT */

	/* QUERY INPUT & VD CELL SEARCH */
	for (int q = 0; q < Q; q++) {
		std::cin >> qry[q];
		if (qry[q].t == 1) {
			int i = search(qry[q].s).IDX;
			int j = search(qry[q].e).IDX;
			qry[q].u = i; qry[q].p = j;
		}
	}
	/* QUERY INPUT & VD CELL SEARCH */



}
int main() { solve(); return 0; }//boj18349 The Creation
/* MAIN */