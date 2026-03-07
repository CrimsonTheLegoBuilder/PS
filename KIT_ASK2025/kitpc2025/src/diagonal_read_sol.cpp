#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cstring>
typedef long long ll;
const int LEN = 1001;

int N, l;
char S[LEN][LEN];
int kumoh_count(const std::string& s) {
	int sz = s.size(), cnt = 0;
	for (int j = 0; j < sz - 4; j++) if (s.substr(j, 5) == "KUMOH") cnt++;
	return cnt;
}
void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);

	for (int i = 0; i < LEN; i++)
		for (int j = 0; j < LEN; j++)
			S[i][j] = 'X';

	std::cin >> N;
	for (int i = 0; i < N; i++) {
		std::string s; std::cin >> s;
		int sz = s.size();
		l = std::max(l, sz);
		for (int j = 0; j < sz; j++) S[i][j] = s[j];
	}

	if (N < 5) { std::cout << "0\n"; return; }

	ll ret = 0;
	int i = 4, j = 1;
	while (j < l - 4) {
		int y, x;
		if (i < N) y = i, x = 0, i++;
		else y = N - 1, x = j, j++;
		std::string s;
		while (y >= 0 && x < l) {
			if (S[y][x] != 'X') s += S[y][x];
			y--; x++;
		}
		int c0 = kumoh_count(s);
		std::reverse(s.begin(), s.end());
		int c1 = kumoh_count(s);
		ret += std::max(c0, c1);
	}
	std::cout << ret << "\n";
	return;
}
int main() { solve(); return 0; }

/*
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

int main(void) {
    cin.tie(0)->sync_with_stdio(0);
    int n = 0;
    cin >> n;
    vector<vector<char>>v(1001, vector<char>(1001, '?'));
    for (int i = 0; i < n; ++i) {
        string str; cin >> str;
        for (int j = 0; j < (int)str.size(); ++j) v[i][j] = str[j];
    }
    int x = 0, y = 0, t = 0, ans = 0;
    auto cal = [&](string& str) {
        int cnt = 0, len = str.size();
        for (int i = 0; i + 4 < len; ++i) {
            string pp;
            for (int j = 0; j <= 4; ++j) pp += str[i + j];
            if (pp == "KUMOH") cnt++;
        }
        return cnt;
    };
    while (x < 1000 && y < 1000) {
        int nx = x, ny = y;
        string temp;
        while (nx >= 0 && nx < 1000 && ny >= 0 && ny < 1000) {
            if (v[nx][ny] != '?') temp += v[nx][ny];
            nx++; ny--;
        }
        int a = 0, b = 0;
        a = cal(temp);
        reverse(temp.begin(), temp.end());
        b = cal(temp);
        ans += max(a, b);
        if (!t) y++;
        else x++;
        if (y == 1000) {
            t = 1;
            x++;
            y = 999;
        }
    }
    cout << ans;
    return 0;
}
*/