#include <bits/stdc++.h>
using namespace std;
#if __cplusplus > 201703L
#include <ranges>
using namespace numbers;
#endif
 
template<class T>
struct graph{
	struct E{
		int from, to;
		T cost;
		friend ostream operator<<(ostream &out, const E &e){
			return out << "{" << e.from << ", " << e.to << ", " << e.cost << "}";
		}
	};
	int n;
	vector<E> edge;
	vector<vector<int>> adj;
	function<bool(int)> ignore;
	graph(int n = 0): n(n), adj(n){ }
	graph(const vector<vector<int>> &adj, bool undirected = true): n((int)adj.size()), adj(n){
		if(undirected){
			for(auto u = 0; u < n; ++ u) for(auto v: adj[u]) if(u < v) link(u, v);
		}
		else for(auto u = 0; u < n; ++ u) for(auto v: adj[u]) orient(u, v);
	}
	graph(const vector<vector<pair<int, T>>> &adj, bool undirected = true): n((int)adj.size()), adj(n){
		if(undirected){
			for(auto u = 0; u < n; ++ u) for(auto [v, w]: adj[u]) if(u < v) link(u, v, w);
		}
		else for(auto u = 0; u < n; ++ u) for(auto [v, w]: adj[u]) orient(u, v, w);
	}
	graph(int n, vector<array<int, 2>> &edge, bool undirected = true): n(n), adj(n){
		for(auto [u, v]: edge) undirected ? link(u, v) : orient(u, v);
	}
	graph(int n, vector<tuple<int, int, T>> &edge, bool undirected = true): n(n), adj(n){
		for(auto [u, v, w]: edge) undirected ? link(u, v, w) : orient(u, v, w);
	}
	int operator()(int u, int id) const{
		#ifdef LOCAL
		assert(0 <= id && id < (int)edge.size());
		assert(edge[id].from == u || edge[id].to == u);
		#endif
		return u ^ edge[id].from ^ edge[id].to;
	}
	int link(int u, int v, T w = {}){ // insert an undirected edge
		int id = (int)edge.size();
		adj[u].push_back(id), adj[v].push_back(id), edge.push_back({u, v, w});
		return id;
	}
	int orient(int u, int v, T w = {}){ // insert a directed edge
		int id = (int)edge.size();
		adj[u].push_back(id), edge.push_back({u, v, w});
		return id;
	}
	graph transposed() const{ // the transpose of the directed graph
		graph res(n);
		for(auto &e: edge) res.orient(e.to, e.from, e.cost);
		res.ignore = ignore;
		return res;
	}
	int degree(int u) const{ // the degree (outdegree if directed) of u (without the ignoration rule)
		return (int)adj[u].size();
	}
	vector<vector<int>> get_adjacency_list() const{
		vector<vector<int>> res(n);
		for(auto u = 0; u < n; ++ u) for(auto id: adj[u]){
			if(ignore && ignore(id)) continue;
			auto &e = edge[id];
			int v = u ^ e.from ^ e.to;
			res[u].push_back(v);
		}
		return res;
	}
	void set_ignoration_rule(const function<bool(int)> &f){
		ignore = f;
	}
	void reset_ignoration_rule(){
		ignore = nullptr;
	}
	friend ostream &operator<<(ostream &out, const graph &g){
		for(auto id = 0; id < (int)g.edge.size(); ++ id){
			if(g.ignore && g.ignore(id)) continue;
			auto &e = g.edge[id];
			out << "{" << e.from << ", " << e.to << ", " << e.cost << "}\n";
		}
		return out;
	}
};
 
template<class T>
struct dfs_forest{
	int n;
	vector<T> dist;
	vector<int> pv;
	vector<int> pe;
	vector<int> order;
	vector<int> pos;
	vector<int> end;
	vector<int> size;
	vector<int> root;
	vector<int> depth;
	vector<int> min_depth;
	vector<int> min_depth_origin;
	vector<int> min_depth_spanning_edge;
	vector<int> was;
	T T_id;
	dfs_forest(int n, T T_id = 0): T_id(T_id){ init(n); }
	void init(int n){
		this->n = n;
		pv.assign(n, -1);
		pe.assign(n, -1);
		order.clear();
		pos.assign(n, -1);
		end.assign(n, -1);
		size.assign(n, 0);
		root.assign(n, -1);
		depth.assign(n, -1);
		min_depth.assign(n, -1);
		min_depth_origin.assign(n, -1);
		min_depth_spanning_edge.assign(n, -1);
		dist.assign(n, T_id);
		was.assign(n, -1);
		attempt = 0;
	}
	int attempt;
	// O(n + m)
	// Requires graph
	template<class Graph, class F = plus<>>
	void dfs(const Graph &g, int u, bool clear_order = true, F UT = plus<>()){
		assert(n == g.n);
		++ attempt;
		depth[u] = 0;
		dist[u] = T_id;
		root[u] = u;
		pv[u] = pe[u] = -1;
		if(clear_order) order.clear();
		auto recurse = [&](auto self, int u)->void{
			was[u] = attempt;
			pos[u] = (int)order.size();
			order.push_back(u);
			size[u] = 1;
			min_depth[u] = depth[u];
			min_depth_origin[u] = u;
			min_depth_spanning_edge[u] = -1;
			for(auto id: g.adj[u]){
				if(id == pe[u] || g.ignore && g.ignore(id)) continue;
				auto &e = g.edge[id];
				int v = e.from ^ e.to ^ u;
				if(was[v] == attempt){
					if(min_depth[u] > depth[v]){
						min_depth[u] = depth[v];
						min_depth_spanning_edge[u] = id;
					}
					continue;
				}
				depth[v] = depth[u] + 1;
				dist[v] = UT(e.cost, dist[u]);
				pv[v] = u;
				pe[v] = id;
				root[v] = root[u];
				self(self, v);
				size[u] += size[v];
				if(min_depth[u] > min_depth[v]){
					min_depth[u] = min_depth[v];
					min_depth_origin[u] = min_depth_origin[v];
					min_depth_spanning_edge[u] = min_depth_spanning_edge[v];
				}
			}
			end[u] = (int)order.size();
		};
		recurse(recurse, u);
	}
	// O(n + m)
	template<class Graph, class F = plus<>>
	void dfs_all(const Graph &g, F UT = plus<>()){
		for(auto u = 0; u < n; ++ u) if(!~depth[u]) dfs<Graph, F>(g, u, false, UT);
	}
	// O(n + m)
	template<class F = plus<>>
	void dfs_implicitly(auto get_deg, auto get_adj, int u, bool clear_order = true, F UT = plus<>()){
		++ attempt;
		depth[u] = 0;
		dist[u] = T_id;
		root[u] = u;
		pv[u] = pe[u] = -1;
		if(clear_order) order.clear();
		auto recurse = [&](auto self, int u)->void{
			was[u] = attempt;
			pos[u] = (int)order.size();
			order.push_back(u);
			size[u] = 1;
			min_depth[u] = depth[u];
			min_depth_origin[u] = u;
			min_depth_spanning_edge[u] = -1;
			for(auto i = 0, deg = get_deg(u); i < deg; ++ i){
				auto [v, w] = get_adj(u, i);
				if(!~v) continue;
				if(was[v] == attempt){
					if(min_depth[u] > depth[v]){
						min_depth[u] = depth[v];
						min_depth_spanning_edge[u] = i;
					}
					continue;
				}
				depth[v] = depth[u] + 1;
				dist[v] = UT(w, dist[u]);
				pv[v] = u;
				pe[v] = i;
				root[v] = root[u];
				self(self, v);
				size[u] += size[v];
				if(min_depth[u] > min_depth[v]){
					min_depth[u] = min_depth[v];
					min_depth_origin[u] = min_depth_origin[v];
					min_depth_spanning_edge[u] = min_depth_spanning_edge[v];
				}
			}
			end[u] = (int)order.size();
		};
		recurse(recurse, u);
	}
	// O(n + m)
	template<class F = plus<>>
	void dfs_all_implicitly(auto get_deg, auto get_adj, F UT = plus<>()){
		for(auto u = 0; u < n; ++ u) if(!~depth[u]) dfs<F>(get_deg, get_adj, u, false, UT);
	}
	bool ancestor_of(int u, int v) const{
		return pos[u] <= pos[v] && end[v] <= end[u];
	}
};
 
template<class T>
struct point{
	T x{}, y{};
	point(){ }
	template<class U> point(const point<U> &otr): x(otr.x), y(otr.y){ }
	template<class U, class V> point(U x, V y): x(x), y(y){ }
	template<class U> point(const array<U, 2> &p): x(p[0]), y(p[1]){ }
	template<class U> operator array<U, 2>() const{
		return {x, y};
	}
	T operator*(const point &otr) const{
		return x * otr.x + y * otr.y;
	}
	T operator^(const point &otr) const{
		return x * otr.y - y * otr.x;
	}
	point operator+(const point &otr) const{
		return {x + otr.x, y + otr.y};
	}
	point &operator+=(const point &otr){
		return *this = *this + otr;
	}
	point operator-(const point &otr) const{
		return {x - otr.x, y - otr.y};
	}
	point &operator-=(const point &otr){
		return *this = *this - otr;
	}
	point operator-() const{
		return {-x, -y};
	}
#define scalarop_l(op) friend point operator op(const T &c, const point &p){ return {c op p.x, c op p.y}; }
	scalarop_l(+) scalarop_l(-) scalarop_l(*) scalarop_l(/)
#define scalarop_r(op) point operator op(const T &c) const{ return {x op c, y op c}; }
	scalarop_r(+) scalarop_r(-) scalarop_r(*) scalarop_r(/)
#define scalarapply(op) point &operator op(const T &c){ return *this = *this op c; }
	scalarapply(+=) scalarapply(-=) scalarapply(*=) scalarapply(/=)
#define compareop(op) bool operator op(const point &otr) const{ return pair<T, T>(x, y) op pair<T, T>(otr.x, otr.y); }
	compareop(>) compareop(<) compareop(>=) compareop(<=) compareop(==) compareop(!=)
#undef scalarop_l
#undef scalarop_r
#undef scalarapply
#undef compareop
	double norm() const{
		return hypot(x, y);
	}
	long double norml() const{
		return hypotl(x, y);
	}
	T squared_norm() const{
		return x * x + y * y;
	}
	// [-pi, pi]
	double angle() const{
		return atan2(y, x);
	}
	// [-pi, pi]
	long double anglel() const{
		return atan2l(y, x);
	}
	point<double> unit() const{
		return point<double>(x, y) / norm();
	}
	point<long double> unitl() const{
		return point<long double>(x, y) / norml();
	}
	point perp() const{
		return {-y, x};
	}
	point<double> normal() const{
		return perp().unit();
	}
	point<long double> normall() const{
		return perp().unitl();
	}
	point<double> rotate(double theta) const{
		return {x * cos(theta) - y * sin(theta), x * sin(theta) + y * cos(theta)};
	}
	point<long double> rotatel(double theta) const{
		return {x * cosl(theta) - y * sinl(theta), x * sinl(theta) + y * cosl(theta)};
	}
	point reflect_x() const{
		return {x, -y};
	}
	point reflect_y() const{
		return {-x, y};
	}
	point reflect(const point &o = {}) const{
		return {2 * o.x - x, 2 * o.y - y};
	}
	bool operator||(const point &otr) const{
		return !(*this ^ otr);
	}
};
template<class T> istream &operator>>(istream &in, point<T> &p){
	return in >> p.x >> p.y;
}
template<class T> ostream &operator<<(ostream &out, const point<T> &p){
	return out << "{" << p.x << ", " << p.y << "}";
}
template<class T>
double distance(const point<T> &p, const point<T> &q){
	return (p - q).norm();
}
template<class T>
long double distancel(const point<T> &p, const point<T> &q){
	return (p - q).norml();
}
template<class T>
T squared_distance(const point<T> &p, const point<T> &q){
	return (p - q).squared_norm();
}
template<class T>
T doubled_signed_area(const point<T> &p, const point<T> &q, const point<T> &r){
	return q - p ^ r - p;
}
template<class T>
T doubled_signed_area(const vector<point<T>> &a){
	if(a.empty()) return 0;
	T res = a.back() ^ a.front();
	for(auto i = 1; i < (int)a.size(); ++ i) res += a[i - 1] ^ a[i];
	return res;
}
// [-pi, pi]
template<class T>
double angle(const point<T> &p, const point<T> &q){
	return atan2(p ^ q, p * q);
}
// [-pi, pi]
template<class T>
long double anglel(const point<T> &p, const point<T> &q){
	return atan2l(p ^ q, p * q);
}
// Check if p->q->r is sorted by angle with respect to the origin
template<class T>
bool is_sorted_by_angle(const point<T> &origin, const point<T> &p, const point<T> &q, const point<T> &r){
	T x = p - origin ^ q - origin;
	T y = q - origin ^ r - origin;
	if(x >= 0 && y >= 0) return true;
	if(x < 0 && y < 0) return false;
	return (p - origin ^ r - origin) < 0;
}
// Check if a is sorted by angle with respect to the origin
template<class T>
bool is_sorted_by_angle(const point<T> &origin, const vector<point<T>> &a){
	for(auto i = 0; i < (int)a.size() - 2; ++ i) if(!is_sorted_by_angle(origin, a[i], a[i + 1], a[i + 2])) return false;
	return true;
}
template<class T>
bool counterclockwise(const point<T> &p, const point<T> &q, const point<T> &r){
	return doubled_signed_area(p, q, r) > 0;
}
template<class T>
bool clockwise(const point<T> &p, const point<T> &q, const point<T> &r){
	return doubled_signed_area(p, q, r) < 0;
}
template<class T>
bool colinear(const point<T> &p, const point<T> &q, const point<T> &r){
	return doubled_signed_area(p, q, r) == 0;
}
template<class T>
bool colinear(const vector<point<T>> &a){
	int i = 1;
	while(i < (int)a.size() && a[0] == a[i]) ++ i;
	if(i >= (int)a.size()) return true;
	for(auto j = i + 1; j < (int)a.size(); ++ j) if(!colinear(a[0], a[i], a[j])) return false;
	return true;
}
 
using pointint = point<int>;
using pointll = point<long long>;
using pointlll = point<__int128_t>;
using pointd = point<double>;
using pointld = point<long double>;
 
// Requires point
template<class P>
struct compare_by_angle{
	const P origin;
	compare_by_angle(const P &origin = P()): origin(origin){ }
	bool operator()(const P &p, const P &q) const{
		return doubled_signed_area(origin, p, q) > 0;
	}
};
template<class It, class P>
void sort_by_angle(It begin, It end, const P &origin){
	begin = partition(begin, end, [&origin](const decltype(*begin) &point){ return point == origin; });
	auto pivot = partition(begin, end, [&origin](const decltype(*begin) &point) { return point > origin; });
	compare_by_angle<P> cmp(origin);
	sort(begin, pivot, cmp), sort(pivot, end, cmp);
}
 
// Given a polygonal region bounded by outer polygon and a set of inner polygonal holes strictly contained in the outer polygon, where none of the polygons intersect, find all trapzoids in the vertical trapzoidal decomposition.
// The outer polygon is given in CCW while the inner ones are in CW.
// Returns {
//  list of trapzoids in increasing order of x-coordinates,
//  list of adjacency info
// }
// where a trapzoid given by the tuple (xl, xr, {i0, j0}, {i1, j1}) is defined by the lines X=xl, X=xr, segment(i0, j0), and segment(i1, j1) where segment(i, j) is the next edge of the j-th vertex of the i-th polygon (i=0 denotes the outer one), and
// an adjacency info given by the tuple (i, j) is that the i-th trapzoid and the j-th trapzoids are adjacent with i-th one on the left.
// T must be able to hold up to max_coordinate^4
// O(n * log(n)) where n is the total number of vertices
// Requires point
template<class T>
pair<
	vector<tuple<T, T, array<int, 2>, array<int, 2>>>,
	vector<array<int, 2>>
>
trapzoidal_decomposition(const vector<vector<point<T>>> &a){
	assert(!a.empty());
	assert(doubled_signed_area(a[0]) >= 0);
	for(auto i = 1; i < (int)a.size(); ++ i) assert(doubled_signed_area(a[i]) <= 0);
	assert(!a.empty());
	assert(doubled_signed_area(a[0]) >= 0);
	for(auto i = 1; i < (int)a.size(); ++ i) assert(doubled_signed_area(a[i]) <= 0);
	vector<array<int, 2>> order;
	for(auto i = 0; i < (int)a.size(); ++ i) for(auto j = 0; j < (int)a[i].size(); ++ j) order.push_back({i, j});
	sort(order.begin(), order.end(), [&](auto i, auto j){ return a[i[0]][i[1]] < a[j[0]][j[1]]; });
	T sweep;
	struct key_type{ // stores the line p-q
		mutable point<T> p, q;
	};
	auto cmp = [&](const key_type &a, const key_type &b)->bool{
		auto ya = a.p.x == a.q.x ? array{max(a.p.y, a.q.y), T(1)} : array{a.p.y * (a.q.x - sweep) + a.q.y * (sweep - a.p.x), a.q.x - a.p.x};
		auto yb = b.p.x == b.q.x ? array{min(b.p.y, b.q.y), T(1)} : array{b.p.y * (b.q.x - sweep) + b.q.y * (sweep - b.p.x), b.q.x - b.p.x};
		if(ya[1] < 0) ya = {-ya[0], -ya[1]};
		if(yb[1] < 0) yb = {-yb[0], -yb[1]};
		return ya[0] * yb[1] < yb[0] * ya[1];
	};
	map<key_type, int, decltype(cmp)> events(cmp);
	vector<tuple<T, T, array<int, 2>, array<int, 2>>> trapzoid;
	vector<array<int, 2>> edge;
	for(auto [i, j]: order){
		auto &b = a[i];
		sweep = b[j].x;
		int pj = (j + (int)b.size() - 1) % (int)b.size();
		int nj = (j + 1) % (int)b.size();
		if(b[j] < b[pj] && b[j] < b[nj]){
			if(doubled_signed_area(b[pj], b[j], b[nj]) > 0){ // Start
				int u = (int)trapzoid.size();
				trapzoid.push_back({sweep, {}, {i, j}, {i, pj}});
				events.insert({{b[j], b[pj]}, u});
			}
			else{ // Split
				auto it = events.lower_bound({b[j], b[j]});
				assert(it != events.end());
				get<1>(trapzoid[it->second]) = sweep;
				int u = (int)trapzoid.size();
				trapzoid.push_back({sweep, {}, get<2>(trapzoid[it->second]), {i, pj}});
				edge.push_back({it->second, u});
				int v = (int)trapzoid.size();
				trapzoid.push_back({sweep, {}, {i, j}, {get<3>(trapzoid[it->second])}});
				edge.push_back({it->second, v});
				it->second = v;
				events.insert(it, {{b[j], b[pj]}, u});
			}
		}
		else if(b[j] > b[pj] && b[j] > b[nj]){ 
			if(doubled_signed_area(b[pj], b[j], b[nj]) > 0){ // End
				auto it = events.lower_bound({b[j], b[j]});
				get<1>(trapzoid[it->second]) = sweep;
				events.erase(it);
			}
			else{ // Merge
				auto l = events.lower_bound({b[j], b[j]}), r = std::next(l);
				get<1>(trapzoid[l->second]) = get<1>(trapzoid[r->second]) = sweep;
				int u = (int)trapzoid.size();
				trapzoid.push_back({sweep, {}, get<2>(trapzoid[l->second]), get<3>(trapzoid[r->second])});
				edge.push_back({l->second, u});
				edge.push_back({r->second, u});
				r->second = u;
				events.erase(l);
			}
		}
		else{ // Regular
			auto it = events.lower_bound({b[j], b[j]});
			get<1>(trapzoid[it->second]) = sweep;
			int u = (int)trapzoid.size();
			if(b[pj] < b[j]){ // Left
				trapzoid.push_back({sweep, {}, {i, j}, get<3>(trapzoid[it->second])});
				edge.push_back({it->second, u});
				it->second = u;
			}
			else{ // Right
				trapzoid.push_back({sweep, {}, {get<2>(trapzoid[it->second])}, {i, pj}});
				edge.push_back({it->second, u});
				it->first.p = it->first.q;
				it->first.q = b[pj];
				it->second = u;
			}
		}
	}
	return {trapzoid, edge};
}
 
int main(){
	cin.tie(0)->sync_with_stdio(0);
	cin.exceptions(ios::badbit | ios::failbit);
	vector<vector<pointll>> a(2);
	const int mx = 3e5 + 1;
	a[0] = {{-1, -1}, {mx + 1, -1}, {mx + 1, mx + 1}, {-1, mx + 1}};
	int n;
	cin >> n;
	a[1].resize(n);
	for(auto j = n - 1; j >= 0; -- j){
		long long x, y;
		cin >> x >> y;
		a[1][j] = {x, y};
	}
	long long minx, miny;
	for(auto rep = 2; rep; -- rep){
		minx = min_element(a[1].begin(), a[1].end())->x;
		auto [trapzoid, edge] = trapzoidal_decomposition(a);
		int m = (int)trapzoid.size();
		graph<int> g(m);
		for(auto [u, v]: edge){
			g.orient(u, v, 1);
			g.orient(v, u, -1);
		}
		dfs_forest<int> df(m);
		df.dfs(g, 0);
		for(auto u = 1; u < m; ++ u){
			auto [xl, xr, yl, yr] = trapzoid[u];
			if(yl[0] == 1 && yr[0] == 1 && g.edge[df.pe[u]].cost == 1){
				minx -= xr - xl;
			}
		}
		for(auto t = 0; t < 2; ++ t){
			for(auto &p: a[t]){
				swap(p.x, p.y);
			}
			reverse(a[t].begin(), a[t].end());
		}
		swap(minx, miny);
	}
	vector<pointll> d(n);
	{
		for(auto i = 0; i < n; ++ i){
			d[i] = a[1][i] - a[1][(i + 1) % n];
		}
		sort_by_angle(d.begin(), d.end(), pointll{});
		auto d2 = d;
		d.clear();
		for(auto p: d2){
			if(d.empty() || d.back() ^ p){
				d.push_back(p);
			}
			else{
				d.back() += p;
			}
		}
		d.insert(d.begin(), pointll{});
		d.pop_back();
		for(auto i = 1; i < (int)d.size(); ++ i){
			d[i] += d[i - 1];
		}
	}
	{
		long long cminy = numeric_limits<long long>::max();
		for(auto &[x, y]: d){
			x += minx;
			cminy = min(cminy, y);
		}
		for(auto &[x, y]: d){
			y += miny - cminy;
		}
	}
	cout << (int)d.size() << "\n";
	for(auto [x, y]: d){
		cout << x << " " << y << "\n";
	}
	return 0;
}
 
/*
 
*/
 
////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//                                   Coded by Aeren                                   //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////
