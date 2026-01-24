#include<iostream>
#include<vector>
#include<cmath>
#include<iomanip>
#include <tuple>
using namespace std;
#define pdd pair<long double, long double>
#define pll pair<long long, long long>
int n, m, q;
vector<pll> convex[2];
long long int area[2][100000];
int change(long long int tmp) {
	if (tmp > 0) return 1;
	if (!tmp) return 0;
	return -1;
}
long long int ccw(pll a, pll b) { return a.first * b.second - a.second * b.first; }
long long int dist(pll a) { return a.first * a.first + a.second * a.second; }
long double dccw(pdd a, pdd b) { return a.first * b.second - a.second * b.first; }
pll tovec(pll a, pll b) { return { b.first - a.first, b.second - a.second }; }
pdd dtovec(pdd a, pdd b) { return { b.first - a.first, b.second - a.second }; }
int inconvex(int idx, pll point) {
	int s = convex[idx].size();
	pll vec_l = tovec(convex[idx][0], convex[idx][s - 1]);
	pll vec_r = tovec(convex[idx][0], convex[idx][1]);
	pll vec_p = tovec(convex[idx][0], point);

	if (ccw(vec_l, vec_p) > 0 || ccw(vec_p, vec_r) > 0) return -1;
	if (ccw(vec_l, vec_p) == 0) {
		if (dist(tovec(convex[idx][0], convex[idx][s - 1])) >= dist(tovec(convex[idx][0], point))) return 0;
		return -1;
	}
	if (ccw(vec_p, vec_r) == 0) {
		if (dist(tovec(convex[idx][0], convex[idx][1])) >= dist(tovec(convex[idx][0], point))) return 0;
		return -1;
	}
	int l = 1, r = s - 1;
	while (r - l > 1) {
		int mid = (l + r) / 2;
		pll vec_m = tovec(convex[idx][0], convex[idx][mid]);
		if (change(ccw(vec_m, vec_p)) == 1) l = mid;
		else r = mid;
	}

	pll vec1 = tovec(convex[idx][l], point);
	pll vec2 = tovec(point, convex[idx][(l + 1) % s]);
	return -change(ccw(vec1, vec2));
}
int f(int idx, pll point) { return change(ccw(tovec(point, convex[1][idx]), tovec(point, convex[1][(idx + 1) % m]))); }
int f2(int idx, pll a, pll b) { return change(ccw(tovec(a, b), tovec(a, convex[0][idx]))); }
pll findrange(pll point) {
	int lb = 0, rb = m - 1, lpoint, rpoint;
	int p = f(0, point);
	int q = f(m - 1, point);
	if (!p) {
		lb = 1;
		p = f(1, point);
	}
	if (!q) {
		rb = m - 2;
		q = f(m - 2, point);
	}
	if (p != q) {
		if (p == 1) {
			rpoint = lb;
			int l = lb, r = rb;
			while (l < r) {
				int mid = (l + r) / 2;
				if (f(mid, point) == 1) l = mid + 1;
				else r = mid;
			}
			lpoint = r;
		}
		else {
			if (f(m - 1, point)) lpoint = 0;
			else lpoint = m - 1;
			int l = lb, r = rb;
			while (l < r) {
				int mid = (l + r) / 2;
				if (f(mid, point) <= 0) l = mid + 1;
				else r = mid;
			}
			rpoint = r;
		}
	}
	else {
		if (p == 1) {
			int l2 = lb, r2 = rb;
			while (l2 < r2) {
				int mid = (l2 + r2) / 2;
				if (change(ccw(tovec(point, convex[1][0]), tovec(point, convex[1][mid]))) >= 0) l2 = mid + 1;
				else r2 = mid;
			}
			int l = lb, r = r2 - 1;
			while (l < r) {
				int mid = (l + r) / 2;
				if (f(mid, point) == 1) l = mid + 1;
				else r = mid;
			}
			lpoint = r;
			l = r2 - 1, r = rb;
			while (l < r) {
				int mid = (l + r) / 2;
				if (f(mid, point) <= 0) l = mid + 1;
				else r = mid;
			}
			rpoint = r;
		}
		else {
			int l2 = lb, r2 = rb;
			while (l2 < r2) {
				int mid = (l2 + r2) / 2;
				if (change(ccw(tovec(point, convex[1][0]), tovec(point, convex[1][mid]))) <= 0) l2 = mid + 1;
				else r2 = mid;
			}
			int l = lb, r = r2 - 1;
			while (l < r) {
				int mid = (l + r) / 2;
				if (f(mid, point) <= 0) l = mid + 1;
				else r = mid;
			}
			rpoint = r;
			l = r2 - 1, r = rb;
			while (l < r) {
				int mid = (l + r) / 2;
				if (f(mid, point) == 1) l = mid + 1;
				else r = mid;
			}
			lpoint = r;
		}
	}
	return { lpoint, rpoint };
}
int findidx(pll a, pll b, int k) {
	int lb = 0, rb = n - 1;
	int p = f2(0, a, b);
	int q = f2(n - 1, a, b);
	if (!p) {
		lb = 1;
		p = f2(1, a, b);
	}
	if (!q) {
		rb = n - 2;
		q = f2(n - 2, a, b);
	}
	if (p != q) {
		if (p == 1) {
			if (k) return (lb - 1 + n) % n;
			if (!f2(n - 1, a, b)) return n - 1;
			return 0;
		}
		else {
			int l = lb, r = rb;
			while (l < r) {
				int mid = (l + r) / 2;
				if (k) {
					if (f2(mid, a, b) <= 0) l = mid + 1;
					else r = mid;
				}
				else {
					if (f2(mid, a, b) < 0) l = mid + 1;
					else r = mid;
				}
			}
			if (k) return (r - 1 + n) % n;
			return r;
		}
	}
	else {
		if (p == 1) {
			int l2 = lb, r2 = rb;
			while (l2 < r2) {
				int mid = (l2 + r2) / 2;
				if (change(ccw(tovec(convex[0][0], a), tovec(convex[0][0], convex[0][mid]))) < 0 || !mid) l2 = mid + 1;
				else r2 = mid;
			}
			int l = r2, r = rb;
			while (l < r) {
				int mid = (l + r) / 2;
				if (k) {
					if (f2(mid, a, b) <= 0) l = mid + 1;
					else r = mid;
				}
				else {
					if (f2(mid, a, b) < 0) l = mid + 1;
					else r = mid;
				}
			}
			if (k) return (r - 1 + n) % n;
			return r;
		}
		else {
			int l2 = lb, r2 = rb;
			while (l2 < r2) {
				int mid = (l2 + r2) / 2;
				if (change(ccw(tovec(convex[0][0], a), tovec(convex[0][0], convex[0][mid]))) <= 0) l2 = mid + 1;
				else r2 = mid;
			}
			int l = lb, r = r2 - 1;
			if (f2(r, a, b) == 0) return r;
			if (f2(r, a, b) < 0) {
				if (k) return r;
				return (r + 1) % n;
			}
			while (l < r) {
				int mid = (l + r) / 2;
				if (k) {
					if (f2(mid, a, b) <= 0) l = mid + 1;
					else r = mid;
				}
				else {
					if (f2(mid, a, b) < 0) l = mid + 1;
					else r = mid;
				}
			}
			if (k) return (r - 1 + n) % n;
			return r;
		}
	}
}
pdd findcross(pll a, pll b, pll c, pll d) {
	long double x1 = a.first, y1 = a.second, x2 = b.first, y2 = b.second, x3 = c.first, y3 = c.second, x4 = d.first, y4 = d.second;
	if (x1 == x2) return { x1, (y4 - y3)* (x1 - x3) / (x4 - x3) + y3 };
	if (x3 == x4) return { x3, (y2 - y1)* (x3 - x1) / (x2 - x1) + y1 };
	long double resx = x1 * (x4 - x3) * (y2 - y1) - x3 * (y4 - y3) * (x2 - x1) + y3 * (x2 - x1) * (x4 - x3) - y1 * (x2 - x1) * (x4 - x3);
	resx /= ((x4 - x3) * (y2 - y1) - (y4 - y3) * (x2 - x1));
	long double resy = (y2 - y1) * (resx - x1) / (x2 - x1) + y1;
	return { resx, resy };
}
long double surface(int l, int r, int idx) {
	if (l < r) return (long double)(area[idx][r] - area[idx][l] - ccw(tovec(convex[idx][0], convex[idx][l]), tovec(convex[idx][0], convex[idx][r])));
	if (!idx) return (long double)(area[idx][n - 1] + area[idx][r] - area[idx][l] + ccw(tovec(convex[idx][0], convex[idx][r]), tovec(convex[idx][0], convex[idx][l])));
	return (long double)(area[idx][m - 1] + area[idx][r] - area[idx][l] + ccw(tovec(convex[idx][0], convex[idx][r]), tovec(convex[idx][0], convex[idx][l])));
}
int main(void) {
	ios_base::sync_with_stdio(false);
	cin.tie(NULL);
	cout.tie(NULL);
	cin >> n >> m >> q;
	long long int x, y;
	for (int i = 0; i < n; i++) {
		cin >> x >> y;
		convex[0].push_back({ x, y });
		if (i >= 2) area[0][i] = area[0][i - 1] + ccw(tovec(convex[0][0], convex[0][i - 1]), tovec(convex[0][0], convex[0][i]));
	}
	for (int i = 0; i < m; i++) {
		cin >> x >> y;
		convex[1].push_back({ x, y });
		if (i >= 2) area[1][i] = area[1][i - 1] + ccw(tovec(convex[1][0], convex[1][i - 1]), tovec(convex[1][0], convex[1][i]));
	}
	for (int i = 0; i < q; i++) {
		cin >> x >> y;
		if (inconvex(0, { x, y }) <= 0) {
			cout << "OUT\n";
			continue;
		}
		if (inconvex(1, { x, y }) >= 0) {
			cout << "IN\n";
			continue;
		}
		pll range = findrange({ x, y });
		int crossright = findidx({ x, y }, convex[1][range.second], 1);
		int crossleft = findidx({ x, y }, convex[1][range.first], 0);
		pdd lcross = findcross({ x, y }, convex[1][range.first], convex[0][crossleft], convex[0][(crossleft - 1 + n) % n]);
		pdd rcross = findcross({ x, y }, convex[1][range.second], convex[0][crossright], convex[0][(crossright + 1) % n]);
		long double res = surface(crossleft, crossright, 0) + dccw(dtovec(convex[0][crossleft], convex[0][crossright]), dtovec(convex[0][crossleft], rcross)) + dccw(dtovec(convex[0][crossleft], rcross), dtovec(convex[0][crossleft], lcross));
		res -= surface(range.first, range.second, 1) + dccw(dtovec(convex[1][range.first], convex[1][range.second]), dtovec(convex[1][range.first], rcross)) + dccw(dtovec(convex[1][range.first], rcross), dtovec(convex[1][range.first], lcross));
		cout << fixed << setprecision(7);
		cout << (area[0][n - 1] - res - area[1][m - 1]) / 2 << '\n';
	}
}