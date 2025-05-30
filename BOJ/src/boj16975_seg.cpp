#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
typedef long long ll;
const int LEN = 100'001;
const int INF = 2e9;
int T, N, M, X, I, J; ll K;
ll arr[LEN];
ll segtree[LEN << 2];

void update(int l, int r, ll diff, int s = 1, int e = N, int i = 1) {
	if (r < s || e < l) return;
	if (l <= s && e <= r) { segtree[i] += diff; return; }
	int m = s + e >> 1;
	update(l, r, diff, s, m, i << 1);
	update(l, r, diff, m + 1, e, i << 1 | 1);
	return;
}
ll search(int idx, int s = 1, int e = N, int i = 1) {
	if (e < idx || idx < s) return 0;
	if (s == e) return segtree[i] + arr[idx];
	int m = s + e >> 1;
	ll L = search(idx, s, m, i << 1);
	ll R = search(idx, m + 1, e, i << 1 | 1);
	return segtree[i] + L + R;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	memset(arr, 0, sizeof arr);
	memset(segtree, 0, sizeof segtree);
	std::cin >> N;
	for (int i = 1; i <= N; i++) std::cin >> arr[i];
	std::cin >> M;
	while (M--) {
		std::cin >> X;
		if (X == 1) {
			std::cin >> I >> J >> K;
			if (I > J) std::swap(I, J);
			update(I, J, K);
		}
		else {
			std::cin >> I;
			std::cout << search(I) << "\n";
		}
	}
	return;
}
int main() { solve(); return 0; }//boj16975

/*

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <queue>
#include <cassert>
typedef long long ll;
const int LEN = 1e5 + 5;
const ll INF = 1e17;

int N, M, K;
int a, b, c;
ll d;
ll arr[LEN], tree[LEN << 2], lazy[LEN << 2];
void lazy_propagate(int s, int e, int i) {
	if (lazy[i]) {
		tree[i] += lazy[i] * ((ll)e - s + 1);
		if (s ^ e) {
			lazy[i << 1] += lazy[i];
			lazy[i << 1 | 1] += lazy[i];
		}
		lazy[i] = 0;
	}
	return;
}
void update(int l, int r, ll diff, int s = 1, int e = N, int i = 1) {
	lazy_propagate(s, e, i);
	if (r < s || e < l) return;
	if (l <= s && e <= r) {
		tree[i] += diff * ((ll)e - s + 1);
		if (s ^ e) {
			lazy[i << 1] += diff;
			lazy[i << 1 | 1] += diff;
		}
		return;
	}
	int m = s + e >> 1;
	update(l, r, diff, s, m, i << 1);
	update(l, r, diff, m + 1, e, i << 1 | 1);
	tree[i] = tree[i << 1] + tree[i << 1 | 1];
	return;
}
ll search(int l, int r, int s = 1, int e = N, int i = 1) {
	lazy_propagate(s, e, i);
	if (r < s || e < l) return 0;
	if (l <= s && e <= r) return tree[i];
	int m = s + e >> 1;
	return search(l, r, s, m, i << 1) + search(l, r, m + 1, e, i << 1 | 1);
}
ll init(int s = 1, int e = N, int i = 1) {
	if (s == e) return tree[i] = arr[s];
	int m = s + e >> 1;
	tree[i << 1] = init(s, m, i << 1);
	tree[i << 1 | 1] = init(m + 1, e, i << 1 | 1);
	return tree[i] = tree[i << 1] + tree[i << 1 | 1];
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	std::cin >> N;
	for (int i = 1; i <= N; i++) std::cin >> arr[i];
	init();
	int cnt = M;
	while (cnt--) {
		std::cin >> a;
		if (a == 1) {
			std::cin >> b >> c >> d;
			if (c < b) std::swap(b, c);
			update(b, c, d);
		}
		else if (a == 2) {
			std::cin >> b;
			std::cout << search(b, b) << "\n";
		}
	}
	return;
}
int main() { solve(); return 0; }//16975


*/
