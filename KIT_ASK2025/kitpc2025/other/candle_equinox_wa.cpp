#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <cassert>
#include <vector>
#include <iomanip>
#include <chrono>
#include <set>
#include <map>
using namespace std;
typedef long long ll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
typedef pair<double, double> pdd;
#define fastio cin.tie(0)->sync_with_stdio(0); cout.tie(0);
#define all(x) x.begin(),x.end()
#define compress(x) x.erase(unique(all(x)),x.end())
#define ff first
#define ss second
#define INF 987654321
#define MAX 500010
#define SIZE 100010

const int N = 3e5 + 9;
const double inf = 1e100;
const double eps = 1e-9;
const double PI = acos(-1.0);
int sign(double x) { return (x > eps) - (x < -eps); }
int sign(int x) { return (x > 0) - (x < 0); }

struct Point {
    double x, y;
    Point() { x = 0; y = 0; }
    Point(double x, double y) : x(x), y(y) {}
    Point(const Point& p) : x(p.x), y(p.y) {}
    Point operator + (const Point& p) const { return Point(x + p.x, y + p.y); }
    Point operator - (const Point& p) const { return Point(x - p.x, y - p.y); }
    Point operator * (const double k) const { return Point(x * k, y * k); }
    friend Point operator * (const double& k, const Point& p) { return Point(k * p.x, k * p.y); }
    Point operator / (const double k) const { return Point(x / k, y / k); }
    bool operator == (const Point& p) const { return sign(x - p.x) == 0 && sign(y - p.y) == 0; }
    bool operator != (const Point& p) const { return !(*this == p); }
    bool operator < (const Point& p) const { return sign(x - p.x) != 0 ? sign(x - p.x) < 0 : sign(y - p.y) < 0; }
    bool operator > (const Point& p) const { return sign(x - p.x) != 0 ? sign(x - p.x) > 0 : sign(y - p.y) > 0; }
    double norm() { return sqrt(x * x + y * y); }
    double norm2() { return x * x + y * y; }
    Point perp() { return Point(-y, x); }
    double arg() { return atan2(y, x); }
    Point truncate(double r) { //returns vector with norm r with the same direction
        double k = norm();
        if (!sign(k)) return *this;
        r /= k;
        return Point(x * r, y * r);
    }
};
istream& operator >> (istream& in, Point& p) { return in >> p.x >> p.y; }
ostream& operator << (ostream& out, Point& p) { return out << "(" << p.x << "," << p.y << ")"; }
inline double dot(Point a, Point b) { return a.x * b.x + a.y * b.y; }
inline double dist2(Point a, Point b) { return dot(a - b, a - b); }
inline double dist(Point a, Point b) { return sqrt(dot(a - b, a - b)); }
inline double cross(Point a, Point b) { return a.x * b.y - a.y * b.x; }
inline double cross2(Point a, Point b, Point c) { return cross(b - a, c - a); }
inline int orientation(Point a, Point b, Point c) { return sign(cross(b - a, c - a)); }

bool half(Point p) {
    return sign(p.y) > 0 || (sign(p.y) == 0 && sign(p.x) < 0);
}
void polar_sort(vector<Point>& v) { // sort points in counterclockwise
    sort(v.begin(), v.end(), [](Point a, Point b) {
        if (half(a) != half(b)) return half(a) < half(b);
        if (sign(cross(a, b)) != 0) return sign(cross(a, b)) > 0;
        return sign(a.norm2() - b.norm2()) < 0;
        });
}
void polar_sort(vector<Point>& v, Point o) { // sort points in counterclockwise with respect to point o
    sort(v.begin(), v.end(), [&](Point a, Point b) {
        if (half(a - o) != half(b - o)) return half(a - o) < half(b - o);
        if (sign(cross(a - o, b - o)) != 0) return sign(cross(a - o, b - o)) > 0;
        return sign((a - o).norm2() - (b - o).norm2()) < 0;
        });
}

// returns true if  point p is on line segment ab
bool is_point_on_seg(Point a, Point b, Point p) {
    if (fabs(cross(p - b, a - b)) < eps) {
        if (p.x < min(a.x, b.x) - eps || p.x > max(a.x, b.x) + eps) return false;
        if (p.y < min(a.y, b.y) - eps || p.y > max(a.y, b.y) + eps) return false;
        return true;
    }
    return false;
}
// intersection point between segment ab and segment cd assuming unique intersection exists
bool seg_seg_intersection(Point a, Point b, Point c, Point d, Point& ans) {
    double oa = cross2(c, d, a), ob = cross2(c, d, b);
    double oc = cross2(a, b, c), od = cross2(a, b, d);
    if (sign(oa * ob) < 0 && sign(oc * od) < 0) {
        ans = (a * ob - b * oa) / (ob - oa);
        return 1;
    }
    else return 0;
}
// intersection point between segment ab and segment cd assuming unique intersection may not exists
// se.size()==0 means no intersection
// se.size()==1 means one intersection
// se.size()==2 means range intersection
set<Point> seg_seg_intersection_inside(Point a, Point b, Point c, Point d) {
    Point ans;
    if (seg_seg_intersection(a, b, c, d, ans)) return { ans };
    set<Point> se;
    if (is_point_on_seg(c, d, a)) se.insert(a);
    if (is_point_on_seg(c, d, b)) se.insert(b);
    if (is_point_on_seg(a, b, c)) se.insert(c);
    if (is_point_on_seg(a, b, d)) se.insert(d);
    return se;
}

bool is_point_on_polygon(vector<Point>& p, const Point& z) {
    int n = p.size();
    for (int i = 0; i < n; i++) {
        if (is_point_on_seg(p[i], p[(i + 1) % n], z)) return 1;
    }
    return 0;
}
// returns 1e9 if the point is on the polygon 
int winding_number(vector<Point>& p, const Point& z) { // O(n)
    if (is_point_on_polygon(p, z)) return 1e9;
    int n = p.size(), ans = 0;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        bool below = p[i].y < z.y;
        if (below != (p[j].y < z.y)) {
            auto orient = orientation(z, p[j], p[i]);
            if (sign(orient) == 0) return 0;
            if (below == (sign(orient) > 0)) ans += below ? 1 : -1;
        }
    }
    return ans;
}
// -1 if strictly inside, 0 if on the polygon, 1 if strictly outside
int is_point_in_polygon(vector<Point>& p, const Point& z) { // O(n)
    int k = winding_number(p, z);
    return k == 1e9 ? 0 : k == 0 ? 1 : -1;
}

typedef vector<Point> Poly;
#define RED 1
#define GREEN 2
#define BLUE 4
#define MAGENTA RED|BLUE
#define YELLOW RED|GREEN
#define CYAN BLUE|GREEN
#define WHITE 7

double poly_area(Poly& p) {
    ll sz = (ll)p.size();
    double ret = 0;
    for (int i = 0; i < sz; i++) {
        ret += cross(p[i], p[(i + 1) % sz]);
    }
    return ret * 0.5;
}

ll n, m, ina[50];
double area[10];
vector<Point> out, in[50], inp, ins;
map<Point, ll> mp;
ll pidx = 0;
vector<ll> g[100010];
map<pll, ll> col;

void get_intersections(Poly& p, Poly& q, Poly& r, int c) {
    ll s1 = (ll)p.size();
    ll s2 = (ll)q.size();
    ll s3 = (ll)r.size();
    for (int i = 0; i < s1; i++) {
        vector<Point> v;
        v.push_back(p[i]); v.push_back(p[(i + 1) % s1]);
        for (int j = 0; j < s2; j++) {
            Point c1 = p[i], n1 = p[(i + 1) % s1];
            Point c2 = q[j], n2 = q[(j + 1) % s2];
            set<Point> st = seg_seg_intersection_inside(c1, n1, c2, n2);
            for (auto t : st) {
                if (!mp[t]) {
                    mp[t] = ++pidx;
                    ins.push_back(t);
                }
                v.push_back(t);
            }
        }
        for (int j = 0; j < s3; j++) {
            Point c1 = p[i], n1 = p[(i + 1) % s1];
            Point c2 = r[j], n2 = r[(j + 1) % s3];
            set<Point> st = seg_seg_intersection_inside(c1, n1, c2, n2);
            for (auto t : st) {
                if (!mp[t]) {
                    mp[t] = ++pidx;
                    ins.push_back(t);
                }
                v.push_back(t);
            }
        }
        sort(all(v));
        compress(v);
        sort(all(v), [&](Point x, Point y) {
            return sign(dist(p[i], x) - dist(p[i], y)) < 0;
            });
        //cout << v.size() << " " << p.size() << "\n";
        for (int j = 1; j < (ll)v.size(); j++) {
            ll p = j - 1;
            ll id1 = mp[v[p]], id2 = mp[v[j]];
            g[id1].push_back(id2);
            g[id2].push_back(id1);
            col[{id1, id2}] |= c;
            //cout << mp[v[j]] << " " << mp[v[p]] << "\n";
        }
    }
}

Poly get_visibility_polygon(Point p) {
    //cout << p << "\n";
    vector<Point> cand = inp;
    polar_sort(cand, p);
    vector<Point> vis;
    for (auto t : cand) {
        ll f = 0;
        for (int i = 0; i < m; i++) {
            ll sz = (ll)in[i].size();
            for (int j = 0; j < sz; j++) {
                Point cur = in[i][j], nxt = in[i][(j + 1) % sz];
                if (cur == t || nxt == t) continue;
                set<Point> st = seg_seg_intersection_inside(cur, nxt, p, t);
                if (st.size()) {
                    f = 1;
                    break;
                }
            }
            if (f) break;
        }
        if (!f) vis.push_back(t);
    }

    vector<Point> ret;

    ll sz = (ll)vis.size();
    Point prev;
    for (int i = 0; i < sz; i++) {
        Point cur = vis[i];
        Point far1 = p + (cur - p) * (ll)1e4;
        vector<Point> v1, v2;
        for (int j = 0; j < n; j++) {
            Point cr = out[j], nx = out[(j + 1) % n];
            set<Point> st = seg_seg_intersection_inside(p, far1, cr, nx);
            for (auto k : st) {
                v1.push_back(k);
                v2.push_back(k);
            }
        }
        for (int t = 0; t < m; t++) {
            ll km = (ll)in[t].size();
            for (int j = 0; j < km; j++) {
                Point cr = in[t][j], nx = in[t][(j + 1) % km];
                set<Point> st = seg_seg_intersection_inside(p, far1, cr, nx);
                if ((ll)st.size() > 1) {
                    Point pr = in[t][(j - 1 + km) % km];
                    if (sign(cross2(p, far1, pr)) > 0) {
                        for (auto k : st) {
                            v1.push_back(k);
                        }
                    }
                    if (sign(cross2(p, far1, pr)) < 0) {
                        for (auto k : st) {
                            v2.push_back(k);
                        }
                    }
                }
                else if ((ll)st.size() == 1) {
                    Point inter = *(st.begin());
                    if (inter != cr && inter != nx) {
                        v1.push_back(inter);
                        v2.push_back(inter);
                    }
                    else if (inter == cr) {
                        Point pr = in[t][(j - 1 + km) % km];
                        if (sign(cross2(p, far1, pr)) > 0 || sign(cross2(p, far1, nx)) > 0) {
                            v1.push_back(inter);
                        }
                        if (sign(cross2(p, far1, pr)) < 0 || sign(cross2(p, far1, nx)) < 0) {
                            v2.push_back(inter);
                        }
                    }
                    else if (inter == nx) {
                        Point pr = in[t][(j + 2) % km];
                        if (sign(cross2(p, far1, pr)) > 0 || sign(cross2(p, far1, cr)) > 0) {
                            v1.push_back(inter);
                        }
                        if (sign(cross2(p, far1, pr)) < 0 || sign(cross2(p, far1, cr)) < 0) {
                            v2.push_back(inter);
                        }
                    }
                }
            }
        }

        sort(all(v1), [&](Point x, Point y) {
            return sign(dist2(x, p) - dist2(y, p)) < 0;
            });
        sort(all(v2), [&](Point x, Point y) {
            return sign(dist2(x, p) - dist2(y, p)) < 0;
            });
        // cout << i << " " << vis[i] << "\n";
        // cout << cur << " " << nxt << "\n";
        // cout << v1.size() << " " << v2.size() << "\n";
        // cout << v2[0] << " " << v1[0] << "\n";
        if (v2[0] == v1[0]) ret.push_back(v1[0]);
        else {
            ret.push_back(v2[0]);
            ret.push_back(v1[0]);
        }
    }

    return ret;
}

vector<vector<ll> > find_faces() {
    vector<vector<int> > used(pidx + 1);
    for (int i = 1; i <= pidx; i++) {
        used[i].resize(g[i].size() + 1);
        used[i].assign(g[i].size() + 1, 0);
        auto compare = [&](int l, int r) {
            Point pl = ins[l - 1] - ins[i - 1];
            Point pr = ins[r - 1] - ins[i - 1];
            if (half(pl) != half(pr)) return half(pl) < half(pr);
            return sign(cross(pl, pr)) > 0;
            };
        sort(all(g[i]), compare);
    }
    vector<vector<ll> > faces;
    for (int i = 1; i <= pidx; i++) {
        for (int e_id = 0; e_id < (ll)g[i].size(); e_id++) {
            if (used[i][e_id]) continue;
            vector<ll> face;
            int v = i;
            int e = e_id;
            while (!used[v][e]) {
                used[v][e] = 1;
                face.push_back(v);
                int u = g[v][e];
                int e1 = lower_bound(all(g[u]), v, [&](int l, int r) {
                    Point pl = ins[l - 1] - ins[u - 1];
                    Point pr = ins[r - 1] - ins[u - 1];
                    if (half(pl) != half(pr)) return half(pl) < half(pr);
                    return sign(cross(pl, pr)) > 0;
                    }) - g[u].begin() + 1;
                if (e1 == (ll)g[u].size()) {
                    e1 = 0;
                }
                v = u;
                e = e1;
            }
            reverse(all(face));
            faces.push_back(face);
        }
    }
    return faces;
}

void solve() {
    pidx = 0;
    mp.clear();
    ins.clear();
    for (int i = 0; i < 8; i++) {
        area[i] = 0;
    }
    for (int i = 0; i < m; i++) {
        ina[i] = 0;
    }
    Point R, G, B;
    cin >> R >> G >> B;
    Poly red, green, blue;
    if (is_point_in_polygon(out, R) == 1) red = vector<Point>();
    else {
        ll f = 0;
        for (int i = 0; i < m; i++) {
            if (is_point_in_polygon(in[i], R) == -1) {
                ina[i] |= RED;
                f = 1;
            }
        }
        if (!f) red = get_visibility_polygon(R);
        else red = vector<Point>();
    }
    if (is_point_in_polygon(out, G) == 1) green = vector<Point>();
    else {
        ll f = 0;
        for (int i = 0; i < m; i++) {
            if (is_point_in_polygon(in[i], G) == -1) {
                ina[i] |= GREEN;
                f = 1;
            }
        }
        if (!f) green = get_visibility_polygon(G);
        else green = vector<Point>();
    }
    if (is_point_in_polygon(out, B) == 1) blue = vector<Point>();
    else {
        ll f = 0;
        for (int i = 0; i < m; i++) {
            if (is_point_in_polygon(in[i], B) == -1) {
                ina[i] |= BLUE;
                f = 1;
            }
        }
        if (!f) blue = get_visibility_polygon(B);
        else blue = vector<Point>();
    }
    for (auto i : red) {
        if (!mp[i]) {
            mp[i] = ++pidx;
            ins.push_back(i);
        }
    }
    for (auto i : blue) {
        if (!mp[i]) {
            mp[i] = ++pidx;
            ins.push_back(i);
        }
    }
    for (auto i : green) {
        if (!mp[i]) {
            mp[i] = ++pidx;
            ins.push_back(i);
        }
    }
    for (int i = 0; i <= 40000; i++) {
        g[i].clear();
    }
    col.clear();
    get_intersections(red, green, blue, RED);
    get_intersections(green, blue, red, GREEN);
    get_intersections(blue, red, green, BLUE);
    for (int i = 0; i <= pidx; i++) {
        sort(all(g[i]));
        compress(g[i]);
    }

    // for(auto i:blue){
    //     cout << setprecision(10);
    //     cout << fixed << i << "\n";
    // }

    // for(auto i:col){
    //     cout << i.ff.ff << " " << i.ff.ss << " " << i.ss << "\n";
    // }

    // cout << ins.size() << "\n";
    // cout << setprecision(10);
    // cout << fixed;
    // for(auto i:ins){
    //     cout << i << " " << mp[i] << "\n";
    // }
    // for(int i=1;i<=pidx;i++){
    //     cout << "(" << i << ") ";
    //     for(auto j:g[i]){
    //         cout << j << " ";
    //     }
    //     cout << "\n";
    // }
    // cout << "\n";

    vector<vector<ll> > faces = find_faces();

    // for(auto i:faces){
    //     for(auto j:i){
    //         cout << j << " ";
    //     }
    //     cout << "\n";
    // }

    for (auto face : faces) {
        vector<Point> vp;
        int c = 0;
        for (int i = 0; i < (ll)face.size(); i++) {
            ll nx = (i + 1) % (ll)face.size();
            c |= col[{face[i], face[nx]}];
            vp.push_back(ins[face[i] - 1]);
        }
        double ar = poly_area(vp);
        // for(auto i:face){
        //     cout << i << " ";
        // }
        // cout << ar << " " << c << "\n";
        if (sign(ar) >= 0 && c != 0) {
            area[c] += abs(ar);
        }
    }
    for (int i = 0; i < m; i++) {
        if (ina[i] != 0) area[ina[i]] += poly_area(in[i]);
    }
    double tot = poly_area(out);
    area[0] = tot;
    for (int i = 1; i < 8; i++) {
        area[0] -= area[i];
    }
    cout << setprecision(10);
    cout << fixed;
    cout << area[RED] << "\n";
    cout << area[GREEN] << "\n";
    cout << area[BLUE] << "\n";
    cout << area[YELLOW] << "\n";
    cout << area[MAGENTA] << "\n";
    cout << area[CYAN] << "\n";
    cout << area[WHITE] << "\n";
    cout << area[0] << "\n";
}

int main() {
    fastio;
    //freopen("../tests/candle/in/16.in", "r", stdin);
    //freopen("../tests/candle/out_eqn/16.out", "w", stdout);
    cin >> n;
    for (int i = 0; i < n; i++) {
        Point p; cin >> p;
        out.push_back(p);
        inp.push_back(p);
    }
    cin >> m;
    for (int i = 0; i < m; i++) {
        ll k; cin >> k;
        for (int j = 0; j < k; j++) {
            Point p; cin >> p;
            in[i].push_back(p);
            inp.push_back(p);
        }
    }

    // Point p = Point(3,3);
    // Poly ret = get_visibility_polygon(p);
    // reverse(all(ret));
    // cout << poly_area(ret) << "\n";
    //ll q; cin >> q;
    ll q = 1;
    while (q--) {
        solve();
    }
}