#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include <set>
//#define TIME
#ifdef TIME
#include <iomanip>
#include <chrono>
#endif
typedef long long ll;
//typedef long double ld;
typedef double ld;
typedef std::pair<int, int> pi;
typedef std::vector<int> Vint;
typedef std::vector<ld> Vld;
const ll INF = 1e17;
const int LEN = 5e3 + 5;

//#define ROTATE

int W, H, N, M, Q;
int XL[LEN], XR[LEN];
int YL[LEN], YH[LEN];
int X1, Y1, X2, Y2;
ll A;
struct Pos {
	int x, y;
	Pos(int x_ = 0, int y_ = 0) : x(x_), y(y_) {}
	Pos operator ~ () const { return { -y, x }; }
	friend std::istream& operator >> (std::istream& is, Pos& p) { is >> p.x >> p.y; return is; }
	friend std::ostream& operator << (std::ostream& os, const Pos& p) { os << p.x << " " << p.y; return os; }
} D0, D1; const Pos O = Pos(0, 0);
typedef std::vector<Pos> Polygon;
Polygon P[LEN];
struct Seg {
	int s, e, h, i;
	Seg(int s_ = 0, int e_ = 0, int h_ = 0, int i_ = -1) : s(s_), e(e_), h(h_), i(i_) {}
	bool operator < (const Seg& o) const {
		if (h == o.h) return s == o.s ? i < o.i : s < o.s;
		return h < o.h;
	}
} B[LEN];
std::vector<Seg> V, S;//segments, stack
struct Bin {
	int i;
	bool operator < (const Bin& o) const { return B[i] < B[o.i]; }
};
typedef std::set<Bin> Tree;
Tree T;//set
struct Event {
	int x, i;
	Event(int x_ = 0, int i_ = 0) : x(x_), i(i_) {}
	bool operator < (const Event& o) const { return x < o.x; }
};
std::vector<Event> E;
Vint X;
ll sweep() {
	ll a = 0;
	std::sort(E.begin(), E.end());
	std::sort(X.begin(), X.end());
	X.erase(unique(X.begin(), X.end()), X.end());
	//std::cout << "x:: ";
	//for (int& x : X) std::cout << x << " ";
	//std::cout << "\n";
	int szx = X.size();
	int sze = E.size();
	for (int x = 0, ei = 0; x < szx - 1; x++) {
		int s = X[x], e = X[x + 1];
		while (ei < sze && E[ei].x == s) {
			//std::cout << "E[ei].i:: " << E[ei].i << "\n";
			auto it = T.find({ E[ei].i });
			if (it == T.end()) T.insert({ E[ei].i });
			else T.erase({ E[ei].i });
			//std::cout << "T::\n";
			//for (auto& t : T) std::cout << t.i << " ";
			//std::cout << "\nT::\n";
			ei++;
		}
		assert(!T.empty());
		//std::cout << "T.sz:: " << T.size() << "\n";
		int h = B[T.begin()->i].h;
		V.push_back(Seg(s, e, h));
	}
	//std::cout << "v::\n";
	//for (Seg& v : V) std::cout << "s:: " << v.s << " e:: " << v.e << " h:: " << v.h << "\n";
	//std::cout << "v::\n";
	S = { Seg(XL[0] - 1, XL[0], 0) };
	for (const Seg& v : V) {
		//std::cout << "v:: s::" << v.s << " e:: " << v.e << " h:: " << v.h << "\n";
		while (S.size() && v.h < S.back().h) {
			ll h = S.back().h;
			S.pop_back();
			ll e = S.back().e;
			a = (v.s - e) * h;
			if (A < a) {
				A = a;
				X1 = e;
				X2 = v.s;
				Y1 = 0;
				Y2 = h;
			}
			//a = std::max(a, (v.s - e) * h);
		}
		S.push_back(v);
	}
	int sz = S.size();
	ll e = XR[0];
	for (int i = 1; i < sz; i++) {
		a = S[i].h * (e - S[i - 1].e);
		if (A < a) {
			A = a;
			X1 = S[i - 1].e;
			X2 = e;
			Y1 = 0;
			Y2 = S[i].h;
		}
		//a = std::max(a, S[i].h * (e - S[i - 1].e));
	}
	E.clear(); X.clear(); V.clear(); S.clear(); T.clear();
	return A;
}
void solve() {
#ifdef TIME
	auto start = std::chrono::high_resolution_clock::now();
#endif
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	//freopen("../tests/awesome/15.in", "r", stdin);
	//freopen("../tests/awesome/15_ans.out", "w", stdout);
	std::cin >> W >> H >> N;
	D0 = Pos(0, 0); D1 = Pos(W, H);
	P[0] = { D0, D1 };
	for (int i = 1; i <= N; i++) {
		Pos a, b; std::cin >> a >> b;
		P[i] = { a, b };
	}
	E.clear(); X.clear(); V.clear(); S.clear(); T.clear();
	ll ret = 0;
#ifdef ROTATE
	for (int _ = 0; _ < 4; _++) {
#endif
		for (int i = 0; i <= N; i++) {
			XL[i] = std::min(P[i][0].x, P[i][1].x);
			XR[i] = std::max(P[i][0].x, P[i][1].x);
			YL[i] = std::min(P[i][0].y, P[i][1].y);
			YH[i] = std::max(P[i][0].y, P[i][1].y);
		}
		for (int q = 0, y, b, h; q <= N; q++) {
			b = (q == 0 ? YL[q] : YH[q]);//bottom
			//std::cout << "\n";
			//std::cout << "b:: " << b << "\n";
			//std::cout << "XL[" << q << "]:: " << XL[q] << "\n";
			//std::cout << "YL[" << q << "]:: " << YL[q] << "\n";
			//std::cout << "XR[" << q << "]:: " << XR[q] << "\n";
			//std::cout << "YH[" << q << "]:: " << YH[q] << "\n";
			for (int i = 0; i <= N; i++) {
				if (i && q == i) continue;
				if (b >= YH[i]) continue;
				y = (i == 0 ? YH[i] : YL[i]);
				h = std::max(0, y - b);//height
				B[i] = { XL[i], XR[i], h, i };
				//std::cout << "h:: " << h << " XL:: " << XL[i] << " XR:: " << XR[i] << "\n";
				E.push_back(Event(XL[i], i));
				E.push_back(Event(XR[i], i));
				X.push_back(XL[i]);
				X.push_back(XR[i]);
			}
			ll tmp = sweep();
			if (ret <= tmp) {
				ret = tmp;
				//std::cout << "Y1:: " << Y1 << " Y2:: " << Y2 << " b:: " << b  << "\n";
				//std::cout << "ret:: " << ret << "\n";
				Y1 = Y1 + b;
				Y2 = Y2 + b;
			}
			//ret = std::max(ret, tmp);
		}
#ifdef ROTATE
		for (int i = 0; i <= N; i++) for (Pos& p : P[i]) p = ~p;
	}
#endif
	std::cout << X1 << " " << Y1 << " " << X2 << " " << Y2 << "\n";
	std::cout << ret << "\n";
#ifdef TIME
	std::cout << "\n\n\n";
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Start time: " << std::chrono::duration_cast<std::chrono::microseconds>(start.time_since_epoch()).count() << " us\n";
	std::cout << "End time: " << std::chrono::duration_cast<std::chrono::microseconds>(end.time_since_epoch()).count() << " us\n";
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "Execution time: " << duration.count() << " microseconds\n";
#endif
	return;
}
int main() { solve(); return 0; }//kitpc? 1X? CSS IS AWESOME