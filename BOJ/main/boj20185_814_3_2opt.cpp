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
#include <set>
#include <chrono>
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::vector<int> Vint;
typedef std::vector<ll> Vll;
typedef std::vector<ld> Vld;
typedef std::vector<bool> Vbool;
const ld INF = 1e18;
const ld TOL = 1e-9;
const ld PI = acos(-1);
const int LEN = 1e4;
inline int sign(const ll& x) { return x < 0 ? -1 : !!x; }
inline int sign(const ld& x) { return x < -TOL ? -1 : x > TOL; }
inline bool zero(const ld& x) { return !sign(x); }
inline bool eq(const ld& x, const ld& y) { return zero(x - y); }
inline ld sq(const ld& x) { return x * x; }
inline ll sq(const ll& x) { return x * x; }
inline ld norm(ld th) { while (th < 0) th += 2 * PI; while (sign(th - 2 * PI) >= 0) th -= 2 * PI; return th; }
inline ld fit(const ld& x, const ld& lo = 0, const ld& hi = 1) { return std::min(hi, std::max(lo, x)); }
inline ll gcd(ll a, ll b) { while (b) { ll tmp = a % b; a = b; b = tmp; } return a; }
inline ll gcd(ll x, ll y, ll z) {
	x = std::abs(x); y = std::abs(y); z = std::abs(z);
	ll w = gcd(x, y);
	return gcd(w, z);
}

bool DEBUG = 0;
bool HILBERT_ONLY = 0;

//#define LOCAL_TEST

/*

tested in range -1e6 < x, y < 1e6;
Delaunator - https://github.com/abellgithub/delaunator-cpp/blob/master/include/delaunator.cpp
modify : jinhwanlazy
I'm : stupid

*/

const int N_ = 8000;
const int K_ = 140;
const int DX = 14;
const int DY = 10;

int N, K;
struct Pii {
	int x, y; int i;
	Pii(int x_ = 0, int y_ = 0, int i_ = -1) : x(x_), y(y_), i(i_) {}
	bool operator == (const Pii& p) const { return x == p.x && y == p.y; }
	bool operator != (const Pii& p) const { return x != p.x || y != p.y; }
	bool operator < (const Pii& p) const { return x == p.x ? y < p.y : x < p.x; }
	bool operator <= (const Pii& p) const { return x == p.x ? y <= p.y : x <= p.x; }
	Pii operator + (const Pii& p) const { return { x + p.x, y + p.y }; }
	Pii operator - (const Pii& p) const { return { x - p.x, y - p.y }; }
	Pii operator * (const int& n) const { return { x * n, y * n }; }
	Pii operator / (const int& n) const { return { x / n, y / n }; }
	ll operator * (const Pii& p) const { return { (ll)x * p.x + (ll)y * p.y }; }
	ll operator / (const Pii& p) const { return { (ll)x * p.y - (ll)y * p.x }; }
	Pii& operator += (const Pii& p) { x += p.x; y += p.y; return *this; }
	Pii& operator -= (const Pii& p) { x -= p.x; y -= p.y; return *this; }
	Pii& operator *= (const int& n) { x *= n; y *= n; return *this; }
	Pii& operator /= (const int& n) { x /= n; y /= n; return *this; }
	Pii operator - () const { return { -x, -y }; }
	Pii operator ~ () const { return { -y, x }; }
	Pii operator ! () const { return { y, x }; }
	ll xy() const { return (ll)x * y; }
	ll Euc() const { return (ll)x * x + (ll)y * y; }
	ll Man() const { return std::abs(x) + std::abs(y); }
	ld mag() const { return hypot(x, y); }
	friend std::istream& operator >> (std::istream& is, Pii& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pii& p) { os << p.x << " " << p.y; return os; }
};
const Pii Oii = { 0, 0 };
const Pii INF_PT = { (int)1e9, (int)1e9 };
typedef std::vector<Pii> Vpii;
bool cmpx(const Pii& p, const Pii& q) { return p.x == q.x ? p.y < q.y : p.x < q.x; }
bool cmpy(const Pii& p, const Pii& q) { return p.y == q.y ? p.x < q.x : p.y < q.y; }
bool cmpi(const Pii& p, const Pii& q) { return p.i < q.i; }
ll cross(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pii& d1, const Pii& d2, const Pii& d3) { return sign(cross(d1, d2, d3)); }
int ccw(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return sign(cross(d1, d2, d3, d4)); }
std::vector<Pii> graham_scan(std::vector<Pii>& C) {
	std::vector<Pii> H;
	if (C.size() < 3) {
		std::sort(C.begin(), C.end());
		return C;
	}
	std::swap(C[0], *min_element(C.begin(), C.end()));
	std::sort(C.begin() + 1, C.end(), [&](const Pii& p, const Pii& q) -> bool {
		int ret = ccw(C[0], p, q);
		if (!ret) return (C[0] - p).Euc() < (C[0] - q).Euc();
		return ret > 0;
		}
	);
	C.erase(unique(C.begin(), C.end()), C.end());
	int sz = C.size();
	for (int i = 0; i < sz; i++) {
		while (H.size() >= 2 && ccw(H[H.size() - 2], H.back(), C[i]) <= 0)
			H.pop_back();
		H.push_back(C[i]);
	}
	return H;
}
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
struct Order {
	ll o;
	int i;
	bool operator < (const Order& q) const { return o < q.o; }
};
ll hilbert_order(const Pii& p, int pow2) {
	int x = p.x, y = p.y;
	ll d = 0;
	for (int s = pow2 >> 1; s > 0; s >>= 1) {
		int rx = (x & s) > 0;
		int ry = (y & s) > 0;
		d += (ll)s * s * ((3ll * rx) ^ ry);
		if (ry == 0) {
			if (rx == 1) {
				x = pow2 - 1 - x;
				y = pow2 - 1 - y;
			}
			int temp = x; x = y; y = temp;
		}
	}
	return d;
}
Vpii clst[DX][DY];
int clst_cnt[K_], clst_cnt_2d[DX][DY];
int fst_belt_cnt[DX];
int grp[LEN];
Vint dt[LEN];
int prv[LEN], nxt[LEN];//graph
bool F[K_];
struct ClusterMeta {
	ld d;
	ll sx;
	ll sy;
} meta[DX][DY];
void update_meta(int x, int y) {
	const Vpii& path = clst[x][y];
	int sz = path.size();
	if (sz == 0) { meta[x][y] = { 0.0, 0, 0 }; return; }
	ld d = 0;
	ll sx = 0, sy = 0;
	for (int i = 0; i < sz; i++) {
		const Pii& p0 = path[i], & p1 = path[(i + 1) % sz];
		d += (p0 - p1).mag();
		sx += p0.x; sy += p0.y;
	}
	meta[x][y] = { d, sx, sy };
	return;
}
void init_all_meta() { for (int x = 0; x < DX; x++) for (int y = 0; y < DY; y++) update_meta(x, y); }
Pii centroid(int x, int y) {
	ll sx = 0, sy = 0;
	int sz = clst[x][y].size();
	if (sz == 0) return { -1, -1 };
	for (const Pii& p : clst[x][y]) { sx += p.x; sy += p.y; }
	return { (int)(sx / sz), (int)(sy / sz) };
}
Pii centroid_fast(int x, int y) {
	int sz = clst[x][y].size();
	if (sz == 0) return { -1, -1 };
	return { (int)(meta[x][y].sx / sz), (int)(meta[x][y].sy / sz) };
}
void first_clustering(Vpii& P) {
	for (int i = 0; i < K; i++) {
		if (!((i + 1) % 7)) clst_cnt[i] = clst_cnt_2d[i / DY][i % DY] = 58;
		else clst_cnt[i] = clst_cnt_2d[i / DY][i % DY] = 57;
	}
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			fst_belt_cnt[x] += clst_cnt_2d[x][y];
		}
	}
	int prv_x = 0;
	int g = 0;
	memset(grp, -1, sizeof grp);
	std::sort(P.begin(), P.end(), cmpx);
	for (int x = 0; x < DX; x++) {
		Vpii C;
		for (int i = 0; i < fst_belt_cnt[x]; i++) C.push_back(P[prv_x + i]);
		std::sort(C.begin(), C.end(), cmpy);
		int prv_y = 0;
		for (int y = 0; y < DY; y++) {
			for (int i = 0; i < clst_cnt_2d[x][y]; i++) {
				grp[C[prv_y + i].i] = g;
				clst[x][y].push_back(C[prv_y + i]);
			}
			g++;
			prv_y += clst_cnt_2d[x][y];
		}
		prv_x += fst_belt_cnt[x];
	}
	std::sort(P.begin(), P.end(), cmpi);
	return;
}
void init_hilbert_paths(Vpii& P) {
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			//if (clst[x][y].empty()) continue;
			std::vector<Order> V;
			int min_x = 1e9, min_y = 1e9;
			int max_x = -1e9, max_y = -1e9;
			int sz = clst[x][y].size();
			for (int i = 0; i < sz; i++) {
				min_x = std::min(min_x, clst[x][y][i].x);
				min_y = std::min(min_y, clst[x][y][i].y);
				max_x = std::max(max_x, clst[x][y][i].x);
				max_y = std::max(max_y, clst[x][y][i].y);
			}
			ll pow2 = 1;
			int tx = max_x - min_x;
			int ty = max_y - min_y;
			int t = std::max(tx, ty);
			while (pow2 <= t) pow2 <<= 1;
			for (int i = 0; i < clst_cnt_2d[x][y]; i++) {
				Pii p = clst[x][y][i] - Pii(min_x, min_y);
				V.push_back({ hilbert_order(p, pow2), clst[x][y][i].i });
			}
			std::sort(V.begin(), V.end());
			for (int i = 0, i0, i2; i < sz; i++) {
				i0 = (i - 1 + sz) % sz;
				i2 = (i + 1) % sz;
				prv[V[i].i] = V[i0].i;
				nxt[V[i].i] = V[i2].i;
			}
		}
	}
	if (HILBERT_ONLY) {
		ld D[K_];
		for (int x = 0; x < DX; x++) {
			for (int y = 0; y < DY; y++) {
				//if (clst[x][y].empty()) continue;
				int s = clst[x][y][0].i;
				int u = s;
				ld d = 0;
				do {
					int v = nxt[u];
					d += (P[u] - P[v]).mag();
					u = v;
				} while (u != s);
				D[x * DY + y] = d;
			}
		}
		std::cout << "HILBERT ONLY::\n";
		for (int i = 0; i < K_; i++) std::cout << D[i] << "\n";
		std::cout << "HILBERT ONLY::\n";
	}
	init_all_meta();
	return;
}
void delaunay_triagulation(const Vpii& P) {
	int sz = P.size();
	Polygon C;
	for (const Pii& p : P) C.push_back(conv(p));
	Delaunator dtr(C);
	for (int i = 0; i < dtr.triangles_.size(); i += 3) {
		const int& a = dtr.points_[dtr.triangles_[i]].i;
		const int& b = dtr.points_[dtr.triangles_[i + 1]].i;
		const int& c = dtr.points_[dtr.triangles_[i + 2]].i;
		dt[a].push_back(b); dt[a].push_back(c);
		dt[b].push_back(a); dt[b].push_back(c);
		dt[c].push_back(a); dt[c].push_back(b);
	}
	for (int i = 0; i < sz; i++) std::sort(dt[i].begin(), dt[i].end());
	return;
}
void optimize_2opt_single(int x, int y, int thr = 50) {
	Vpii& path = clst[x][y];
	int sz = path.size();
	if (sz < 3) return;
	bool imp = 1;
	while (imp && thr--) {
		imp = 0;
		for (int i = 0; i < sz - 1; i++) {
			for (int j = i + 2; j < sz; j++) {
				if (i == 0 && j == sz - 1) continue;
				Pii p1 = path[i];
				Pii p2 = path[i + 1];
				Pii p3 = path[j];
				Pii p4 = path[(j + 1) % sz];
				ld d12 = (p1 - p2).mag();
				ld d34 = (p3 - p4).mag();
				ld d13 = (p1 - p3).mag();
				ld d24 = (p2 - p4).mag();
				if (d13 + d24 < d12 + d34 - TOL) {
					std::reverse(path.begin() + i + 1, path.begin() + j + 1);
					imp = 1;
				}
			}
		}
	}
	return;
}
void optimize_2opt(int thr = 100000) {
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			Vpii& path = clst[x][y];
			int sz = path.size();
			if (sz < 3) continue;
			bool imp = 1;
			while (imp && thr--) {
				imp = 0;
				for (int i = 0; i < sz - 1; i++) {
					for (int j = i + 2; j < sz; j++) {
						if (j == sz - 1 && i == 0) continue;
						Pii p1 = path[i];
						Pii p2 = path[i + 1];
						Pii p3 = path[j];
						Pii p4 = path[(j + 1) % sz];
						ld d12 = (p1 - p2).mag();
						ld d34 = (p3 - p4).mag();
						ld d13 = (p1 - p3).mag();
						ld d24 = (p2 - p4).mag();
						if (d13 + d24 < d12 + d34 - TOL) {
							std::reverse(path.begin() + i + 1, path.begin() + j + 1);
							imp = true;
						}
					}
				}
			}
			for (int i = 0; i < sz; i++) {
				int u = path[i].i;
				int v = path[(i + 1) % sz].i;
				int w = path[(i - 1 + sz) % sz].i;
				nxt[u] = v;
				prv[u] = w;
			}
		}
	}
	init_all_meta();
	return;
}
bool try_move_point(int tx, int ty, int hx, int hy) {
	Vpii& T = clst[tx][ty];
	Vpii& H = clst[hx][hy];
	int sz = T.size();
	if (sz <= 3) return 0;
	if (H.empty()) return 0;
	Pii cen = { (int)(meta[hx][hy].sx / sz), (int)(meta[hx][hy].sy / sz) };
	int best_idx = -1;
	ll min_dist = 4e18;
	sz = T.size();
	for (int i = 0; i < sz; i++) {
		ll d = sq((ll)T[i].x - cen.x) + sq((ll)T[i].y - cen.y);
		if (d < min_dist) {
			min_dist = d;
			best_idx = i;
		}
	}
	if (best_idx == -1) return 0;
	Vpii BT = T;//backup
	Vpii BH = H;//backup
	Pii p = T[best_idx];
	T.erase(T.begin() + best_idx);
	ld best_inc = 1e18;
	int best_pos = 0;
	sz = H.size();
	for (int i = 0; i < sz; i++) {
		const Pii& u = H[i], & v = H[(i + 1) % sz];
		ld inc = (u - p).mag() + (p - v).mag() - (u - v).mag();
		if (inc < best_inc) {
			best_inc = inc;
			best_pos = i;
		}
	}
	H.insert(H.begin() + best_pos + 1, p);
	optimize_2opt_single(tx, ty);
	optimize_2opt_single(hx, hy);

	ld new_t_dist = 0; for (int i = 0; i < T.size(); i++) new_t_dist += (T[i] - T[(i + 1) % T.size()]).mag();
	ld new_h_dist = 0; for (int i = 0; i < H.size(); i++) new_h_dist += (H[i] - H[(i + 1) % H.size()]).mag();

	ld old_max = std::max(meta[tx][ty].d, meta[hx][hy].d);
	ld new_max = std::max(new_t_dist, new_h_dist);

	if (new_max < old_max - 1e-5) {
		//optimize_2opt_single(tx, ty);
		update_meta(tx, ty);
		//optimize_2opt_single(hx, hy);
		update_meta(hx, hy);
		return 1;
	}
	clst[tx][ty] = BT;
	clst[hx][hy] = BH;
	return 0;
}
void balancing_step() {
	ld mx = -1;
	int tx = -1, ty = -1;
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			if (meta[x][y].d > mx) {
				mx = meta[x][y].d;
				tx = x; ty = y;
			}
		}
	}
	Pii cen = centroid_fast(tx, ty);
	ld mn = 1e18;
	int hx = -1, hy = -1;
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			if (x == tx && y == ty) continue;
			Pii c = centroid_fast(x, y);
			ll d = (cen - c).Euc();
			if (d > 25000LL * 25000LL) continue;
			if (meta[x][y].d < mn) {
				mn = meta[x][y].d;
				hx = x; hy = y;
			}
		}
	}
	if (hx != -1) {
		if (try_move_point(tx, ty, hx, hy)) {
			optimize_2opt_single(tx, ty);
			optimize_2opt_single(hx, hy);
		}
	}
	return;
}
void reset_data() {
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			clst[x][y].clear();
		}
	}
	memset(clst_cnt, 0, sizeof(clst_cnt));
	memset(clst_cnt_2d, 0, sizeof(clst_cnt_2d));
	memset(fst_belt_cnt, 0, sizeof(fst_belt_cnt));
	memset(grp, 0, sizeof(grp));
	memset(prv, -1, sizeof(prv));
	memset(nxt, -1, sizeof(nxt));
	return;
}
ld run_solver() {
	std::cin >> N >> K;
	Vpii P(N);
	for (Pii& p : P) std::cin >> p;
	for (int i = 0; i < N; i++) P[i].i = i;
	first_clustering(P);
	init_hilbert_paths(P);
	optimize_2opt();
	auto start_time = std::chrono::steady_clock::now();
	while (true) {
		auto now = std::chrono::steady_clock::now();
		if (std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count() > 4000) break;
		balancing_step();
	}
	ld max_dist = 0;
	int min_idx = 1e9, max_idx = -1;
	std::set<int> S;
	for (int x = 0; x < DX; x++) {
		for (int y = 0; y < DY; y++) {
			if (clst[x][y].empty()) continue;
			std::cout << clst[x][y].size() << "\n";
			for (const Pii& p : clst[x][y]) std::cout << p.i << " ";
			std::cout << "\n";
			//int s = clst[x][y][0].i;
			//int u = s;
			//ld cur_dist = 0;
			//int safety = 0;
			//do {
			//	min_idx = std::min(min_idx, u);
			//	max_idx = std::max(max_idx, u);
			//	S.insert(u);
			//	int v = nxt[u];
			//	cur_dist += sqrt(pow(P[u].x - P[v].x, 2) + pow(P[u].y - P[v].y, 2));
			//	u = v;
			//	safety++;
			//} while (u != s && safety < N + 5);
			//max_dist = std::max(max_dist, cur_dist);
		}
	}
#ifdef LOCAL_TEST
	std::cout << "DEBUG::\n";
	std::cout << "min:: " << min_idx << "\n";
	std::cout << "max:: " << max_idx << "\n";
	std::cout << "sz :: " << S.size() << "\n";
	std::cout << "DEBUG::\n";
#endif
	return max_dist;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cout << std::fixed;
	std::cout.precision(10);
#ifdef LOCAL_TEST
	ld total_score = 0;
	ld max_score = 0;
	std::cout << "========= [LOCAL TEST START] =========\n";
	for (int i = 1; i <= 50; i++) {
		std::string filename = (i < 10 ? "0" : "") + std::to_string(i) + ".in";
		std::string path = "../../tests/814_3/" + filename;
		if (freopen(path.c_str(), "r", stdin) == NULL) { std::cout << "File Not Found: " << path << "\n"; continue; }
		reset_data();
		ld score = run_solver();
		std::cout << "Test " << filename << " : " << score << "\n";
		total_score += score;
		max_score = std::max(max_score, score);
	}
	std::cout << "======================================\n";
	std::cout << "2-opt:\n";
	std::cout << "Total Score: " << total_score << "\n";
	std::cout << "Average Score: " << total_score / 50.0 << "\n";
	std::cout << "Worst Case   : " << max_score << "\n";
	std::cout << "======================================\n";
#else
	reset_data();
	run_solver();
#endif
	return;
}
int main() { solve(); return 0; }
