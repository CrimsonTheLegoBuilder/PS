#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <queue>
#include <numeric>
#include <unordered_map>
typedef long long ll;
typedef double db;
typedef long double ld;
const db INF = 1e18;
const db TOL = 1e-7;
const db PI = acos(-1);
const int LEN = 1e5 + 5;
int N, M, T, Q;
ll V[LEN];
int P[LEN * 6];
inline bool zero(const db& x) { return std::abs(x) < TOL; }
inline int sign(const db& x) { return x < -TOL ? -1 : x > TOL; }
inline db sqr(const db& x) { return x * x; }
ll gcd(ll a, ll b) { return !b ? a : gcd(b, a % b); }

/*

tested in range -1e5 < x, y < 1e5;
Delaunator - https://github.com/abellgithub/delaunator-cpp/blob/master/include/delaunator.cpp
modify : jinhwanlazy
I'm : stupid

*/

struct Info {
    int u, v;
    db c;
    Info(int u = 0, int v = 0, db c = 0) : u(u), v(v), c(c) {}
    bool operator < (const Info& x) const { return c > x.c; }
};
int find(int v) {
    if (P[v] < 0) return v;
    int p = find(P[v]);
    return P[v] = p;
}
int join(int u, int v) {
    int i, j;
    i = find(u);
    j = find(v);
    if (i == j) return 0;
    if (P[i] < P[j]) {
        P[i] += P[j];
        P[j] = i;
    }
    else {
        P[j] += P[i];
        P[i] = j;
    }
    return 1;
}
struct Pii {
    int x, y, i;
    Pii(int X = 0, int Y = 0, int I = 0) : x(X), y(Y), i(I) {}
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
    Pii& operator *= (const int& scale) { x *= scale; y *= scale; return *this; }
    Pii& operator /= (const int& scale) { x /= scale; y /= scale; return *this; }
    Pii operator - () const { return { -x, -y }; }
    Pii operator ~ () const { return { -y, x }; }
    inline Pii operator ! () const { return { y, x, 0 }; }
    ll xy() const { return (ll)x * y; }
    ll Euc() const { return (ll)x * x + (ll)y * y; }
    int Man() const { return std::abs(x) + std::abs(y); }
    //ld mag() const { return hypot(x, y); }
    ld mag() const { return sqrt(Euc()); }
    friend std::istream& operator >> (std::istream& is, Pii& p) { is >> p.x >> p.y; return is; }
    friend std::ostream& operator << (std::ostream& os, const Pii& p) { os << p.x << " " << p.y; return os; }
};
const Pii Oii = { 0, 0 };
const Pii INF_PT = { (int)INF, (int)INF };
ll cross(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) / (d3 - d2); }
ll cross(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) / (d4 - d3); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3) { return (d2 - d1) * (d3 - d2); }
ll dot(const Pii& d1, const Pii& d2, const Pii& d3, const Pii& d4) { return (d2 - d1) * (d4 - d3); }
int ccw(const Pii& d1, const Pii& d2, const Pii& d3) {
    ll ret = cross(d1, d2, d3);
    return !ret ? 0 : ret > 0 ? 1 : -1;
}
namespace std {
    template<>
    struct hash<Pii> {
        std::size_t operator()(const Pii& p) const {
            std::size_t h1 = std::hash<int>()(p.x);
            std::size_t h2 = std::hash<int>()(p.y);
            return h1 ^ (h2 << 1);
        }
    };
}
struct Pos {
    db x, y, i;
    Pos(db X = 0, db Y = 0, int I = 0) : x(X), y(Y), i(I) {}
    bool operator == (const Pos& p) const { return zero(x - p.x) && zero(y - p.y); }
    bool operator != (const Pos& p) const { return !zero(x - p.x) || !zero(y - p.y); }
    bool operator < (const Pos& p) const { return zero(x - p.x) ? y < p.y : x < p.x; }
    Pos operator + (const Pos& p) const { return { x + p.x, y + p.y }; }
    Pos operator - (const Pos& p) const { return { x - p.x, y - p.y }; }
    Pos operator * (const db& scalar) const { return { x * scalar, y * scalar }; }
    Pos operator / (const db& scalar) const { return { x / scalar, y / scalar }; }
    db operator * (const Pos& p) const { return { x * p.x + y * p.y }; }
    inline db operator / (const Pos& p) const { return { x * p.y - y * p.x }; }
    db operator ^ (const Pos& p) const { return { x * p.y - y * p.x }; }
    Pos operator - () const { return Pos(-x, -y); }
    Pos operator ~ () const { return { -y, x }; }
    Pos operator ! () const { return { y, x }; }
    Pos& operator += (const Pos& p) { x += p.x; y += p.y; return *this; }
    Pos& operator -= (const Pos& p) { x -= p.x; y -= p.y; return *this; }
    Pos& operator *= (const db& scale) { x *= scale; y *= scale; return *this; }
    Pos& operator /= (const db& scale) { x /= scale; y /= scale; return *this; }
    db xy() const { return x * y; }
    Pos rot(db the) { return { x * cos(the) - y * sin(the), x * sin(the) + y * cos(the) }; }
    inline db Euc() const { return x * x + y * y; }
    inline db mag() const { return sqrt(Euc()); }
    Pos unit() const { return *this / mag(); }
    db rad() const { return atan2(y, x); }
    friend db rad(const Pos& p1, const Pos& p2) { return atan2l(p1 / p2, p1 * p2); }
    int quad() const { return sign(y) == 1 || (sign(y) == 0 && sign(x) >= 0); }
    friend bool cmpq(const Pos& a, const Pos& b) { return (a.quad() != b.quad()) ? a.quad() < b.quad() : a / b > 0; }
    bool close(const Pos& p) const { return zero((*this - p).Euc()); }
    inline bool close(const Pos& rhs,
        const db span = 1.,
        const db tol = 1e-20) const {
        return ((*this - rhs).Euc() / span) < tol;
    }
    inline friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
    friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
}; const Pos O = { 0, 0 };
typedef std::vector<Pos> Polygon;
inline db cross(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) / (d3 - d2); }
inline db cross(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) / (d4 - d3); }
inline int ccw(const Pos& d1, const Pos& d2, const Pos& d3) {
    db ret = cross(d1, d2, d3);
    return zero(ret) ? 0 : ret > 0 ? 1 : -1;
}
inline bool counterclockwise(const Pos& p, const Pos& q, const Pos& r) {
    return ccw(p, q, r) == 1;
}
inline bool clockwise(const Pos& p, const Pos& q, const Pos& r) {
    return ccw(p, q, r) == -1;
}
inline bool collinear(const Pos& p, const Pos& q, const Pos& r) {
    return !ccw(p, q, r);
}
inline bool collinear(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) {
    return !ccw(d1, d2, d3) && !ccw(d1, d2, d4);
}
inline db dot(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d2); }
inline db dot(const Pos& d1, const Pos& d2, const Pos& d3, const Pos& d4) { return (d2 - d1) * (d4 - d3); }
bool on_seg_strong(const Pos& d1, const Pos& d2, const Pos& d3) {
    db ret = dot(d1, d3, d2);
    return !ccw(d1, d2, d3) && (ret > 0 || zero(ret));
}
bool on_seg_weak(const Pos& d1, const Pos& d2, const Pos& d3) {
    db ret = dot(d1, d3, d2);
    return !ccw(d1, d2, d3) && ret > 0;
}
db projection(const Pos& d1, const Pos& d2, const Pos& d3) { return (d2 - d1) * (d3 - d1) / (d2 - d1).mag(); }
db dist(const Pos& d1, const Pos& d2, const Pos& t) {
    return cross(d1, d2, t) / (d1 - d2).mag();
}
Pos intersection(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2) {
    db a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
    return (p1 * a2 + p2 * a1) / (a1 + a2);
}
inline bool inner_check(const Pos& a, const Pos& b, const Pos& c, const Pos& t) {
    return ccw(a, b, t) > -1 && ccw(b, c, t) > -1 && ccw(c, a, t) > -1;
}
inline bool in_triangle(Pos a, Pos b, Pos c, const Pos& t) {
    if (ccw(a, b, c) < 0) std::swap(b, c);
    return inner_check(a, b, c, t);
}
inline db circumradius(const Pos& p1, const Pos& p2, const Pos& p3) {
    Pos d = p2 - p1;
    Pos e = p3 - p1;

    const db bl = d.Euc();
    const db cl = e.Euc();
    const db det = d / e;

    Pos radius((e.y * bl - d.y * cl) * 0.5 / det,
        (d.x * cl - e.x * bl) * 0.5 / det);

    if ((bl > 0.0 || bl < 0.0) && (cl > 0.0 || cl < 0.0) &&
        (det > 0.0 || det < 0.0))
        return radius.Euc();
    return (std::numeric_limits<double>::max)();
}
inline Pos circumcenter(const Pos& p1, const Pos& p2, const Pos& p3) {
    Pos d = p2 - p1;
    Pos e = p3 - p1;

    const db bl = d.Euc();
    const db cl = e.Euc();
    const db det = d / e;

    Pos radius((e.y * bl - d.y * cl) * 0.5 / det,
        (d.x * cl - e.x * bl) * 0.5 / det);
    return p1 + radius;
}
struct Circle {
    Pos c;
    db r;
    Circle(Pos C = Pos(0, 0), db R = 0) : c(C), r(R) {}
    bool operator == (const Circle& C) const { return c == C.c && std::abs(r - C.r) < TOL; }
    bool operator != (const Circle& C) const { return !(*this == C); }
    bool operator < (const Circle& q) const {
        db dist = (c - q.c).mag();
        return r < q.r && dist + r < q.r + TOL;
    }
    bool operator > (const Pos& p) const { return r > (c - p).mag(); }
    bool operator >= (const Pos& p) const { return r + TOL > (c - p).mag(); }
    bool operator < (const Pos& p) const { return r < (c - p).mag(); }
    Circle operator + (const Circle& C) const { return { c + C.c, r + C.r }; }
    Circle operator - (const Circle& C) const { return { c - C.c, r - C.r }; }
    db H(const db& th) const { return sin(th) * c.x + cos(th) * c.y + r; }//coord trans | check right
    db A() const { return 1. * r * r * PI; }
    friend std::istream& operator >> (std::istream& is, Circle& c) { is >> c.c >> c.r; return is; }
    friend std::ostream& operator << (std::ostream& os, const Circle& c) { os << c.c << " " << c.r; return os; }
} INVAL = { { 0, 0 }, (std::numeric_limits<db>::max)() };
Circle enclose_circle(const Pos& u, const Pos& v, const Pos& w) {
    if (!ccw(u, v, w)) return INVAL;
    Pos m1 = (u + v) * .5, v1 = ~(v - u);
    Pos m2 = (u + w) * .5, v2 = ~(w - u);
    Pos c = intersection(m1, m1 + v1, m2, m2 + v2);
    return Circle(c, (u - c).mag());
}
Circle circumcircle(const Pos& p1, const Pos& p2, const Pos& p3) {
    Pos d = p2 - p1;
    Pos e = p3 - p1;

    const db bl = d.Euc();
    const db cl = e.Euc();
    const db det = d / e;

    Pos radius((e.y * bl - d.y * cl) * 0.5 / det,
        (d.x * cl - e.x * bl) * 0.5 / det);

    Pos c;
    db r;
    if ((bl > 0.0 || bl < 0.0) && (cl > 0.0 || cl < 0.0) &&
        (det > 0.0 || det < 0.0))
        c = p1 + radius, r = radius.Euc();
    else c = Pos(INF, INF), r = std::numeric_limits<db>::max();

    return Circle(c, r);
}
inline bool in_circle(const Pos& a, const Pos& b, const Pos& c, const Pos& p) {
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
        Pos pPrev{ std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };

        // Go through points based on distance from the center.
        for (std::size_t k = 0; k < n; k++) {
            const std::size_t i = ids[k];
            const Pos p = points_[i];

            // skip near-duplicate points
            if (k > 0 && p == pPrev)
                continue;
            pPrev = p;

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
        std::vector<std::size_t> edgesStack;

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
                    a = edgesStack[i];
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

                if (i < edgesStack.size()) {
                    edgesStack[i] = br;
                }
                else {
                    edgesStack.push_back(br);
                }
                i++;

            }
            else {
                if (i > 0) {
                    i--;
                    a = edgesStack[i];
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
void query() {
    auto answer = [&](const db& r) -> ll {
        if (r < 2 + TOL) return 0;
        ll ans = sqr(r - 2) * PI;// +TOL;
        return ans;
        };
    memset(P, -1, sizeof P);
    std::cin >> N;
    std::vector<Pos> C(N);
    db ret = INF;
    for (int i = 0; i < N; ++i) std::cin >> C[i], C[i].i = i, ret = std::min(ret, (O - C[i]).mag());

    //soldiers are out of mines' convex hull
    if (N < 3) { std::cout << answer(ret) << "\n"; return; }

    Delaunator d(C);

    int s = -1, e = d.triangles_.size() / 3;
    std::vector<Info> info;
    std::unordered_map<Pii, int> MAP;
    for (int i = 0, j = 0; i < d.triangles_.size(); i += 3, j++) {
        int a = d.triangles_[i];
        int b = d.triangles_[i + 1];
        int c = d.triangles_[i + 2];
        if (a > b) std::swap(a, b);
        if (b > c) std::swap(b, c);
        if (a > b) std::swap(a, b);
        Pii e1 = Pii(a, b), e2 = Pii(a, c), e3 = Pii(b, c);
        auto edge = { e1, e2, e3 };
        for (const Pii& E : edge) {
            auto it = MAP.find(E);
            if (it != MAP.end()) {
                int u = it->second, v = j;
                info.push_back(Info(u, v, (C[E.x] - C[E.y]).mag() * .5));
                MAP.erase(it);
            }
            else MAP[E] = j;
        }
        if (in_triangle(C[a], C[b], C[c], Pos(0, 0))) s = j;
    }

    //soldiers are out of mines' convex hull
    if (s == -1) { std::cout << answer(ret) << "\n"; return; }

    for (auto& i : MAP) info.push_back(Info(e, i.second, (C[i.first.x] - C[i.first.y]).mag() * .5));

    std::sort(info.begin(), info.end());
    for (const Info& i : info) {
        if (find(s) == find(e)) break;
        join(i.u, i.v);
        ret = std::min(ret, i.c);
    }

    std::cout << answer(ret) << "\n";
    return;
}
void solve() {
    std::cin.tie(0)->sync_with_stdio(0);
    std::cout.tie(0);
    //freopen("../../../input_data/mine_tests/01.in", "r", stdin);
    //freopen("../../../input_data/mine_tests/out.txt", "w", stdout);
    std::cin >> T;
    while (T--) query();
    return;
}
int main() { solve(); return 0; }//BAPC 2009 C Escape from the Minefield boj5401
//refer to jinhwanlazy, JusticeHui, koosaga
