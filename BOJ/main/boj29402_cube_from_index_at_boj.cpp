#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <algorithm>
#include <cassert>
#include <time.h>
#include <random>
#include <vector>
#include <unordered_map>
typedef long long ll;

/*
https://github.com/thregium/practice_baekjoon/blob/main/24000-24999/24902.c
boj의 index 님의 솔루션 코드로 공부해서 제출합니다.
Thistlethwaite's 4-phase Algorithm(https://www.jaapsch.net/puzzles/thistle.htm)
*/

#define PHASE2S 85
#define PHASE3S 1077072625834
#define PHASE4C 16434824
#define PHASE4E 205163983024656
//각 숫자의 의미가 뭘까?

int N;//query num
int S[54];//현재 큐브의 상태
int pre[6];//이전 상태를 저장하는 작은 배열
int res[256];//큐브를 돌리는 순서
int len_count[256];//큐브를 돌린 횟수 각각의 횟수 ?
ll Q[3342336];//상태 저장용 큐
int comb[16][16];//combination
/*
*        +------+
*        | 1 2 3|
*        | 8 0 4|
*        | 7 6 5|
* +------+------+------+------+
* |373839|101112|192021|282930|
* |443640|17 913|261822|352731|
* |434241|161514|252423|343332|
* +------+------+------+------+
*        |464748|
*        |534549|
*        |525150|
*        +------+
* 
* 입력의 순서를 구현에 맞게 바꿔주는 배열
*/
const int tile_od[1 << 6] = {
	             1,  2,  3,
	             8,  0,  4,
	             7,  6,  5,
	37, 38, 39, 10, 11, 12, 19, 20, 21, 28, 29, 30,
	44, 36, 40, 17,  9, 13, 26, 18, 22, 35, 27, 31,
	43, 42, 41, 16, 15, 14, 25, 24, 23, 34, 33, 32,
	            46, 47, 48,
	            53, 45, 49,
	            52, 51, 50
};
void debug_print() {
	std::cout << "DEBUG::\n";
	for (int i = 0; i < 3; i++) {
		std::cout << "      ";
		for (int j = 0; j < 3; j++) {
			std::cout << S[tile_od[i * 3 + j]] << " ";
		}
		std::cout << "\n";
	}
	for (int i = 0; i < 3; i++) {
		std::cout << "      ";
		for (int j = 0; j < 12; j++) {
			std::cout << S[tile_od[9 + i * 12 + j]] << " ";
		}
		std::cout << "\n";
	}
	for (int i = 0; i < 3; i++) {
		std::cout << "      ";
		for (int j = 0; j < 3; j++) {
			std::cout << S[tile_od[45 + i * 3 + j]] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "DEBUG::\n";
	return;
}
const char M[10] = "UFRBLD";
void print_sol(int len) {
	for (int i = 0; i < len; i++) {
		assert(res[i] > 0 && res[i] < 60);
		int f = res[i] / 10;
		int p = res[i] % 10;
		std::cout << M[f];
		assert(1 <= p && p <= 3);
		if (p == 2) std::cout << "2";
		else if (p == 3) std::cout << "'";
		std::cout << " ";
	}
	std::cout << "\n";
	return;
}
/*
*      >----------v
*      | +------+ |
*      | | 1 2 3| |
*      | | 8 0 4| |
*      *<| 7 6 5|-<
* +------+------+------+------+
* |373839|101112|192021|282930|
* |443640|17 913|261822|352731|
* |434241|161514|252423|343332|
* +------+------+------+------+
*        |464748|
*        |534549|
*        |525150|
*        +------+
*
* 각 면의 번호에 대해 시계방향으로 옆면 번호들을 매핑
*/
const int side_od[6][12] = {
	//{ 30, 29, 28, 21, 20, 19, 12, 11, 10, 39, 38, 37 },
	{ 12, 11, 10, 39, 38, 37, 30, 29, 28, 21, 20, 19 },
	{  7,  6,  5, 19, 26, 25, 48, 47, 46, 41, 40, 39 },
	{  5,  4,  3, 28, 35, 34, 50, 49, 48, 14, 13, 12 },
	{  3,  2,  1, 37, 44, 43, 52, 51, 50, 23, 22, 21 },
	{  1,  8,  7, 10, 17, 16, 46, 53, 52, 32, 31, 30 },
	{ 16, 15, 14, 25, 24, 23, 34, 33, 32, 43, 42, 41 }
};
/*
*        +------+
*        |   2  |
*        | 8   4|
*        |   6  |
* +------+------+------+------+
* |  38  |  11  |  20  |  29  |
* |44  40|17  13|26  22|35  31|
* |  42  |  15  |  24  |  33  |
* +------+------+------+------+
*        |  47  |
*        |53  49|
*        |  51  |
*        +------+
*
* 각 엣지에 면 번호 매칭
*/
const int edge_od[12][2] = {
	{  6, 11 },
	{  4, 20 },
	{  2, 29 },
	{  8, 38 },
	{ 47, 15 },
	{ 49, 24 },
	{ 51, 33 },
	{ 53, 42 },
	{ 13, 26 },
	{ 22, 35 },
	{ 31, 44 },
	{ 40, 17 }
};
/*
*        +------+
*        |      |
*        |      |
*        |      |
* +------+------+------+------+
* |      |      |      |      |
* |      |      |      |      |
* |      |      |      |      |
* +------+------+------+------+
*        |      |
*        |      |
*        |      |
*        +------+
* 이건 뭘까?
* 각 색상 쌍에 대응하는 엣지 조각의 번호를 매칭
* -1은 대응하는 조각이 없음
* 색상 쌍은 앞쪽 색 * 8 + 뒤쪽 색
*/
const int color_edge[64] = {
	-1, -1, -1, -1, -1, -1, -1, -1,
	-1, -1,  0,  1,  2,  3, -1, -1,
	-1,  0, -1,  0, -1, 11,  4, -1,
	-1,  1,  8, -1,  9, -1,  5, -1,
	-1,  2, -1,  9, -1, 10,  6, -1,
	-1,  3, 11, -1, 10, -1,  7, -1,
	-1, -1,  4,  5,  6,  7, -1, -1,
	-1, -1, -1, -1, -1, -1, -1, -1
};
/*
*        +------+
*        | 1   3|
*        |      |
*        | 7   5|
* +------+------+------+------+
* |37  39|10  12|19  21|28  30|
* |      |      |      |      |
* |43  41|16  14|25  23|34  32|
* +------+------+------+------+
*        |46  48|
*        |      |
*        |52  50|
*        +------+
* 
* 코너 조각에 면 번호 매칭
*/
const int corner_od[8][3] = {
	{ 39, 10, 7 },
	{ 41, 46, 16 },
	{ 43, 32, 52 },
	{ 37,  1, 30 },
	{ 19,  5, 12 },
	{ 25, 14, 48 },
	{ 23, 50, 34 },
	{ 21, 28,  3 },
};
/*
*      >----------v
*      | +------+ |
*      | |   2  | |
*      | | 3   1| |
*      *<|   0  |-<
* +------+------+------+------+
* |   3  |   0  |   1  |   2  |
* |10  11|11   8| 8   9| 9  10|
* |   7  |   4  |   5  |   6  |
* +------+------+------+------+
*        |   4  |
*        | 7   5|
*        |   6  |
*        +------+
*
* phase1, 2에서 돌아가는 엣지 조각의 순서
*/
const int oper_edge[6][4] = {
	//{ 2,  1, 0,  3 },
	{ 0,  3, 2,  1 },
	{ 0,  8, 4, 11 },
	{ 1,  9, 5,  8 },
	{ 2, 10, 6,  9 },
	{ 3, 11, 7, 10 },
	{ 4,  5, 6,  7 }
};


/* PHASE 1 */
/*
* 페이즈 1
* 각 모서리 블록의 방향을 바꾼다
* 위아래 면을 짝수번 돌려서 원래 있어야 하는 상태로 옮길 수 있는 경우로 만든다
* ZZ공식이랑 비슷한 느낌
* 사용하는 동작: U, D, F(1, 2, 3), R(1, 2, 3), B(1, 2, 3), L(1, 2, 3)
* 경우의 수 2,048
*/
int phase1_last[4096], phase1_oper[4096], phase1_dist[4096];
/*
*      >----------v
*      | +------+ |
*      | |   2  | |
*      | | 3   1| |
*      *<|   0  |-<
* +------+------+------+------+
* |   3  |   0  |   1  |   2  |
* |10  11|11   8| 8   9| 9  10|
* |   7  |   4  |   5  |   6  |
* +------+------+------+------+
*        |   4  |
*        | 7   5|
*        |   6  |
*        +------+
*
* phase1에서 돌아가는 엣지 조각의 순서
*/
const int phase1_edge_oper[6][4] = {
	//{ 2,  1, 0,  3 },
	{ 0,  3, 2,  1 },
	{ 0,  8, 4, 11 },
	{ 1,  9, 5,  8 },
	{ 2, 10, 6,  9 },
	{ 3, 11, 7, 10 },
	{ 4,  5, 6,  7 }
};
/* PHASE 1 */


/* PHASE 2 */
int phase2_corner[6561][6][3], phase2_edge[495][6][3];
char phase2_corner_vis[6561], phase2_edge_vis[495];
int phase2_last[3342336];
char phase2_oper[3342336];
/*
*      >----------v
*      | +------+ |
*      | | 3   7| |
*      | |      | |
*      *<| 0   4|-<
* +------+------+------+------+
* | 3   0| 0   4| 4   7| 7   3|
* |      |      |      |      |
* | 2   1| 1   5| 5   6| 6   2|
* +------+------+------+------+
*        | 1   5|
*        |      |
*        | 2   6|
*        +------+
*
* phase2에서 돌아가는 코너 조각의 순서
* 이건 뭔 순서인지 모르겠다
*/
const int phase2_corner_oper[6][4] = {
	{ 0, 3, 7, 4 },
	{ 0, 4, 5, 1 },
	{ 4, 7, 6, 5 },
	{ 2, 6, 7, 3 },
	{ 0, 1, 2, 3 },
	{ 1, 5, 6, 2 }
};
/*
*      >----------v
*      | +------+ |
*      | |   2  | |
*      | | 3   1| |
*      *<|   0  |-<
* +------+------+------+------+
* |   3  |   0  |   1  |   2  |
* |10  11|11   8| 8   9| 9  10|
* |   7  |   4  |   5  |   6  |
* +------+------+------+------+
*        |   4  |
*        | 7   5|
*        |   6  |
*        +------+
*
* phase2에서 돌아가는 엣지 조각의 순서
* phase1과 동일하기 때문에 하나만 있어도 될 듯 함
*/
const int phase2_edge_oper[6][4] = {
	//{ 2,  1, 0,  3 },
	{ 0,  3, 2,  1 },
	{ 0,  8, 4, 11 },
	{ 1,  9, 5,  8 },
	{ 2, 10, 6,  9 },
	{ 3, 11, 7, 10 },
	{ 4,  5, 6,  7 }
};
/* PHASE 2 */


/* PHASE 3 */
int phase3_corner_perm_cnt = 0, phase3_corner_perm_chk[8];
int phase3_corner[40320][6][3], phase3_edge[70][6][3];
int phase3_corner_perm[40320];
char phase3_corner_vis[40320], phase3_edge_vis[70];
int phase3_last[2822400];
char phase3_oper[2822400], phase3_is_end[2822400];
/* PHASE 3 */


/* PHASE 4 */
int phase4_edge_perm_cnt = 0, phase4_edge_perm_chk[12];
int phase4_corner[96][6][3], phase4_edge[13824][6][3];
int phase4_corner_perm[96];
ll phase4_edge_perm[13824];
char phase4_corner_vis[96], phase4_edge_vis[13824];
int phase4_last[1327104];
char phase4_oper[1327104];
/*
*        +------+
*        |   2  |
*        | 3   1|
*        |   0  |
* +------+------+------+------+
* |   3  |   0  |   1  |   2  |
* |10  11|11   8| 8   9| 9  10|
* |   7  |   4  |   5  |   6  |
* +------+------+------+------+
*        |   4  |
*        | 7   5|
*        |   6  |
*        +------+
* 이건 뭘까?
* phase4일 때 엣지 조각의 가능한 상태
*/
const int oper4_edge[12][12] = {
	{ 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
	{ 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 },
	{ 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
	{ 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 },
	{ 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
	{ 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 },
	{ 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
	{ 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1 }
};
/* PHASE 4 */


void move_cw(const int& f /* face */) {
	//옆면 회전
	pre[0] = S[side_od[f][9]];
	pre[1] = S[side_od[f][10]];
	pre[2] = S[side_od[f][11]];
	for (int i = 11; i >= 3; i--) S[side_od[f][i]] = S[side_od[f][i - 3]];
	for (int i = 0; i < 3; i++) S[side_od[f][i]] = pre[i];
	//윗면 회전
	pre[1] = S[f * 9 + 7];
	pre[2] = S[f * 9 + 8];
	for (int i = 8; i >= 3; i--) S[f * 9 + i] = S[f * 9 + (i - 2)];
	S[f * 9 + 1] = pre[1];
	S[f * 9 + 2] = pre[2];
	return;
}
void move_ccw(const int& f /* face */) {
	//옆면 회전
	pre[0] = S[side_od[f][0]];
	pre[1] = S[side_od[f][1]];
	pre[2] = S[side_od[f][2]];
	for (int i = 0; i < 9; i++) S[side_od[f][i]] = S[side_od[f][i + 3]];
	for (int i = 9; i < 12; i++) S[side_od[f][i]] = pre[i - 9];
	//윗면 회전
	pre[1] = S[f * 9 + 1];
	pre[2] = S[f * 9 + 2];
	for (int i = 1; i <= 6; i++) S[f * 9 + i] = S[f * 9 + (i + 2)];
	S[f * 9 + 7] = pre[1];
	S[f * 9 + 8] = pre[2];
	return;
}
void move_180(const int& f /* face */) {
	//옆면 회전
	for (int i = 0; i < 6; i++) std::swap(S[side_od[f][i]], S[side_od[f][i + 6]]);
	//윗면 회전
	for (int i = 1; i <= 4; i++) std::swap(S[f * 9 + i], S[f * 9 + (i + 4)]);
	return;
}
void move(const int& o) {
	assert(0 < o && o < 60);
	int f = o / 10, p = o % 10;
	assert(1 <= p && p <= 3);
	if (p == 1) move_cw(f);
	if (p == 2) move_180(f);
	if (p == 3) move_ccw(f);
	return;
}
void solve(const int& len) {
	for (int i = 0; i < len; i++) move(res[i]);
	return;
}

void random25() {
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 9; j++)
			S[i * 9 + j] = i + 1;
	for (int i = 0; i < 25; i++) move(rand() % 6 * 10 + rand() % 3 + 1);
	return;
}
bool check() {
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 9; j++)
			if (S[i * 9 + j] != i + 1) return 0;
	return 1;
}


/* PHASE 1 */
//엣지 조각의 상태를 비트마스킹하여 정수형으로 반환.
//1인 엣지는 오리엔테이션이 반대
int check_edge_parity() {
	int ret = 0;
	for (int i = 0; i < 8; i++) {//위 아래 8개
		if (S[edge_od[i][1]] == 1 || S[edge_od[i][1]] == 6) {//윗면이나 아랫면의 오리엔테이션이 맞지 않다면
			ret ^= ((S[edge_od[i][0]] ^ i ^ 1) & 1) << i;//옆면의 위치가 정위치 혹은 반대 위치인지 확인
		}
		else ret ^= ((S[edge_od[i][1]] ^ i) & 1) << i;
	}
	for (int i = 8; i < 12; i++) {//옆 4개
		if (S[edge_od[i][1]] == 1 || S[edge_od[i][1]] == 6) {//
			ret ^= ((S[edge_od[i][0]] ^ i) & 1) << i;
		}
		else ret ^= ((S[edge_od[i][1]] ^ i ^ 1) & 1) << i;
	}
	return ret;
}
/*
* phase1 전처리
* bfs로 가능한 형태 전부 탐색하는 것이 목표
*/
void prec_phase1() {//bfs?
	int x, y, bit;//x: 현재 상태, y: 변할 상태, bit: 선택한 비트
	int qf = 0, qr = 0;
	memset(phase1_last, -1, sizeof phase1_last);
	phase1_last[0] = 0;
	Q[qf++] = 0;
	while (qf > qr) {
		x = Q[qr++];
		for (int i = 0; i < 6; i++) {
			y = x;
			if (i == 0) y ^= 15;//윗면, 00001111
			else if (i == 5) y ^= 240;//아랫면, 11110000
			for (int j = 1; j <= 3; j++) {
				bit = ((y >> oper_edge[i][3]) & 1);
				for (int k = 3; k >= 1; k--) {
					y -= (y & (1 << oper_edge[i][k]));
					y += (((y >> oper_edge[i][k - 1]) & 1) << oper_edge[i][k]);
				}
				y -= (y & (1 << oper_edge[i][0]));
				y += (bit << oper_edge[i][0]);
				if ((i == 0 || i == 5) && j == 2) continue;
				if (phase1_last[y] < 0) {
					Q[qf++] = y;
					phase1_last[y] = x;
					phase1_oper[y] = i * 10 + (4 - j);
					phase1_oper[y] = phase1_dist[x] + 1;
				}
			}
		}
	}
	return;
}
/* PHASE 1 */


/* PHASE 2 */
//이건 또 뭘까
int fb_edge_to_perm(int* fb_edge) {
	int ret = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = (i ? fb_edge[i - 1] + 1 : 0); j < fb_edge[i]; j++) {
			ret += comb[12 - j - 1][3 - i];
		}
	}
	return ret;
}
int get_phase_num_2() {
	int ret = 0, fb_edge[4], tedge, edgecnt = 0;
	for (int i = 7; i >= 0; i--) {
		ret *= 3;
		for (int j = 0; j < 3; j++) {
			if (S[corner_od[i][j]] == 3 || S[corner_od[i][j]] == 5) {
				ret += j;
				break;
			}
		}
	}
	ret *= 495;
	for (int i = 0; i < 12; i++) {
		tedge = color_edge[S[edge_od[i][0]]] * 8 + S[edge_od[i][1]];
		assert(tedge >= 0);
		if (tedge <= 6 && !(tedge & 1)) {
			fb_edge[edgecnt++] = i;
		}
	}
	assert(edgecnt == 4);
	ret += fb_edge_to_perm(fb_edge);
	return ret;
}
int pack_phase_num_2(ll lnum) {
	int fb_edge[4];
	int ecnt = 0;
	int corner;
	int corner_ret = 0;
	for (int i = 0; i < 12; i++) {
		if ((lnum >> i) & 1) fb_edge[ecnt++] = i;
	}
	assert(ecnt == 4);//뭐가 4개인걸까?
	corner = lnum >> 16;
	for (int i = 7; i >= 0; i--) {
		corner_ret *= 3;
		corner_ret += ((corner >> (1 << 1)) & 3);
	}
	return corner_ret * 495 + fb_edge_to_perm(fb_edge);
}
ll phase2_corner_move(ll corner_num, int f) {
	ll store, cornertrit;//이 변수들은 또 뭘까?
	if (f == 1 || f == 3) {
		for (int i = 0; i < 4; i++) {
			cornertrit = ((corner_num >> (phase2_corner_oper[f][i] << 1)) & 3);
			corner_num -= cornertrit << (phase2_corner_oper[f][i] << 1);
			if (i & 1) cornertrit--;//홀수면 1 뺀다
			else cornertrit++;
			if (cornertrit >= 3) cornertrit -= 3;
			else if (cornertrit < 0) cornertrit += 3;
			corner_num += cornertrit << (phase2_corner_oper[f][i] << 1);
		}
	}
	store = ((corner_num >> (phase2_corner_oper[f][3] << 1)) & 3);//뭔가 두 자리를 남겨놔야했다.
	for (int i = 3; i >= 1; i--) {
		corner_num -= (corner_num & (3ll << (phase2_corner_oper[f][i] << 1)));
		corner_num += (((corner_num >> (phase2_corner_oper[f][i - 1] << 1)) & 3) << (phase2_corner_oper[f][i] << 1));
	}
	corner_num -= (corner_num & (3ll << (phase2_corner_oper[f][0] << 1)));
	corner_num += (store << (phase2_corner_oper[f][0] << 1));
	return corner_num;
}
ll phase2_edge_move(ll edge_num, int f) {
	ll store;
	store = ((edge_num >> phase2_edge_oper[f][3]) & 1);
	for (int i = 3; i >= 1; i--) {
		edge_num -= (edge_num & (1ll << phase2_edge_oper[f][i]));
		edge_num += (((edge_num >> phase2_edge_oper[f][i - 1]) & 1) << phase2_edge_oper[f][i]);
	}
	edge_num -= (edge_num & (1ll << phase2_edge_oper[f][0]));
	edge_num += (store << (1ll << phase2_edge_oper[f][0]));
	return edge_num;
}
//페이즈2의 코너 조각부터 전처리한다.
void prec_phase2_corner() {
	int qf = 0, qr = 0;
	ll x = 0, xpk, y, ypk;//x: 현재 상태, y: 바뀔 상태
	memset(phase2_corner, -1, sizeof phase2_corner);
	Q[qf++] = x;
	while (qf > qr) {//bfs
		x = Q[qr++];
		xpk = pack_phase_num_2((x << 16) + 15);
		assert(xpk % 495 == 0);
		xpk /= 495;
		for (int i = 0; i < 6; i++) {
			y = x;
			for (int j = 0; j < 3; j++) {
				y = phase2_corner_move(y, i);
				//if (j != i && (i == 0 || i == 5)) continue;//위아래는 짝수번만 돌릴 수 있음. 첨부 논문 참고
				ypk = pack_phase_num_2((y << 16) + 15);
				assert(ypk % 495 == 0);
				ypk /= 495;
				assert(0 <= ypk && ypk < 6561);
				phase2_corner[xpk][i][j] == ypk;
				if (!phase2_corner_vis[ypk]) {
					Q[qf++] = y;
					phase2_corner_vis[ypk] = 1;
				}
			}
		}
	}
}
//페이즈2의 엣지 조각을 전처리한다.
void prec_phase2_edge() {
	int qf = 0, qr = 0;
	ll x = PHASE2S, xpk, y, ypk;//
	memset(phase2_edge, -1, sizeof phase2_edge);
	Q[qf++] = x;
	while (qf > qr) {//bfs
		x = Q[qr++];
		xpk = pack_phase_num_2(x);
		for (int i = 0; i < 6; i++) {
			y = x;
			for (int j = 0; j < 3; j++) {
				y = phase2_edge_move(y, i);
				//if (j != i && (i == 0 || i == 5)) continue;//위아래는 짝수번만 돌릴 수 있음. 첨부 논문 참고
				ypk = pack_phase_num_2(y);
				assert(0 <= ypk && ypk < 495);
				phase2_edge[xpk][i][j] == ypk;
				if (!phase2_edge_vis[ypk]) {
					Q[qf++] = y;
					phase2_edge_vis[ypk] = 1;
				}
			}
		}
	}
}
//페이즈2 전처리
void prec_phase2() {
	ll x, y;//x: 현재 상태, y: 바뀔 상태
	ll xpk, ypk;//각각을 패킹한 값
	ll corner, edge;//코너, 엣지 블록의 값
	ll qf = 0, qr = 0;
	prec_phase2_corner();
	prec_phase2_edge();
	memset(phase2_oper, -1, sizeof phase2_oper);
	Q[qf++] = pack_phase_num_2(PHASE2S);
	x = pack_phase_num_2(PHASE2S);
	phase2_oper[x] = 0;
	while (qf > qr) {
		x = Q[qr++];
		xpk = (x >> 16) * 495 + (x & 65535);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 3; j++) {
				if (j != 1 && (i == 0 || i == 5)) continue;//위아래는 짝수번만
				y = x;
				corner = (y >> 16);
				edge = (y & 65535);
				corner = phase2_corner[corner][i][j];
				edge = phase2_edge[edge][i][j];
				y = (corner << 16) + edge;
				ypk = (y >> 16) * 495 + (y & 65535);
				assert(ypk >= 0 && ypk < 3247695);
				if (phase2_oper[ypk] < 0) {
					phase2_oper[ypk] = i * 10 + (3 - j);
					phase2_last[ypk] = xpk;
					Q[qf++] = y;
				}
			}
		}
	}
	return;
}
/* PHASE 2 */


/* PHASE 3 */
//?
void prec_phase3_corner_perm(int cp, int cnt) {
	if (cnt < 0) {
		phase3_corner_perm[phase3_corner_perm_cnt++] = cp;//?
		return;
	}
	for (int i = 0; i < 8; i++) {
		if (phase3_corner_perm_chk[i]) continue;
		phase3_corner_perm_chk[i] = 1;
		prec_phase3_corner_perm(cp + (i << (cnt * 3)), cnt - 1);//재귀
		phase3_corner_perm_chk[i] = 0;
	}
	return;
}
//코너 조각의 해시값을 구한다
int get_corner_num_3(int cnum) {
	int c[3], h = 0;
	for (int i = 0; i < 3; i++) c[i] = S[corner_od[cnum][i]];
	std::sort(c, c + 3);
	h = c[0] * 100 + c[1] * 10 + c[2];
	if (h == 125) return 0;
	else if (h == 256) return 1;
	else if (h == 456) return 2;
	else if (h == 145) return 3;
	else if (h == 123) return 4;
	else if (h == 236) return 5;
	else if (h == 346) return 6;
	else if (h == 134) return 7;
	assert(0);
	return -1;
}
//8진법 형태의 코너 조각 순열을 이분 탐색으로 순열 번호로 바꿈
int corner_to_perm_3(int c) {
	int s = 0, e = 40319, m;
	while (s < e) {
		m = s + e >> 1;
		if (phase3_corner_perm[m] >= c) e = m;
		else s = m + 1;
	}
	assert(phase3_corner_perm[s] == c);
	return s;
}
//좌우 엣지 위치의 순서를 바꾼다. 상하 모서리는 제외하고 생각한다.
int lr_edge_to_perm(int* lr_edge) {
	int ret = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = (i ? lr_edge[i - 1] + 1 : 0); j < lr_edge[i]; j++) {
			ret += comb[8 - j - 1][3 - i];
		}
	}
	return ret;
}
int get_phase_num_3() {
	int ret = 0, corn_num = 0, lr_edge[4], tedge, edge_minus = 0, edge_cnt = 0;
	for (int i = 0; i < 8; i++) {
		corn_num = (get_corner_num_3(i) << (i * 3));
	}
	ret = corner_to_perm_3(corn_num);
	ret *= 70;
	for (int i = 0; i < 12; i++) {
		if (i <= 6 && (i & 1)) {
			edge_minus++;
			continue;
		}
		tedge = color_edge[S[edge_od[i][0]] * 8 + S[edge_od[i][1]]];
		assert(tedge >= 0);
		if (tedge <= 7 && (tedge & 1)) {
			lr_edge[edge_cnt++] = i - edge_minus;
		}
	}
	assert(edge_cnt == 4);
	ret += lr_edge_to_perm(lr_edge);
	return ret;
}
int pack_phase_num_3(ll lnum) {
	int lr_edge[4];
	ll edge_minus = 0, ecnt = 0, corner, corner_res = 0;
	for (int i = 0; i < 12; i++) {
		if (i <= 6 && (~i & 1)) {
			edge_minus++;
			continue;
		}
		if ((lnum >> i) & 1) lr_edge[ecnt++] = i - edge_minus;
	}
	assert(ecnt == 4);
	corner = lnum >> 16;
	corner_res = corner_to_perm_3(corner);
	return corner_res * 70 + lr_edge_to_perm(lr_edge);
}
//페이즈 3, 4에서 사용 가능한 코너 조각 순열의 이동
ll phase3_corner_move(ll corner_num, int f) {
	ll store;
	store = ((corner_num >> (phase2_corner_oper[f][3] * 3)) & 7);//111
	for (int i = 3; i >= 1; i--) {
		corner_num -= (corner_num & (7ll << (phase2_corner_oper[f][i] * 3)));
		corner_num += (((corner_num >> (phase2_corner_oper[f][i - 1] * 3)) & 7)
			<< (phase2_corner_oper[f][i] * 3));
	}
	corner_num -= (corner_num & (7ll << (phase2_corner_oper[f][0] * 3)));
	corner_num += (store << (phase2_corner_oper[f][0] * 3));
	return corner_num;
}
//페이즈3의 코너 전처리
void prec_phase3_corner() {
	ll qf = 0, qr = 0, x = (PHASE3S >> 16), xpk, y, ypk;
	memset(phase3_corner, -1, sizeof phase3_corner);
	Q[qf++] = x;
	while (qf > qr) {
		x = Q[qr++];
		xpk = corner_to_perm_3(x);
		for (int i = 0; i < 6; i++) {
			y = x;
			for (int j = 0; j < 3; j++) {
				y = phase3_corner_move(y, i);
				if (j != 1 && (i != 2 && i != 4)) continue;//양 옆을 제외하면 짝수번만 돌릴 수 있음
				ypk = corner_to_perm_3(y);
				assert(ypk >= 0 && ypk < 40320);
				phase3_corner[xpk][i][j] = ypk;
				if (!phase3_corner_vis[ypk]) {
					Q[qf++] = y;
					phase3_corner_vis[ypk] = 1;
				}
			}
		}
	}
	return;
}
//페이즈3의 엣지 전처리
void prec_phase3_edge() {
	ll qf = 0, qr = 0, x = PHASE3S, xpk, y, ypk;
	memset(phase3_edge, -1, sizeof phase3_edge);
	Q[qf++] = x;
	while (qf > qr) {
		x = Q[qr++];
		xpk = pack_phase_num_3(x) - 40319 * 70;
		for (int i = 0; i < 6; i++) {
			y = x;
			for (int j = 0; j < 3; j++) {
				y = phase2_edge_move(y, i);
				if (j != 1 && (i != 2 && i != 4)) continue;//양 옆을 제외하면 짝수번만 돌릴 수 있음
				ypk = pack_phase_num_3(y) - 40319 * 70;
				assert(ypk >= 0 && ypk < 70);
				phase3_edge[xpk][i][j] = ypk;
				if (!phase3_edge_vis[ypk]) {
					Q[qf++] = y;
					phase3_edge_vis[ypk] = 1;
				}
			}
		}
	}
	return;
}
void prec_phase3() {
	ll x, y, xpk, ypk, corner, edge, qf = 0, qr = 0;
	prec_phase3_corner_perm(0, 7);
	prec_phase3_corner();
	prec_phase3_edge();
	
	memset(phase3_oper, -1, sizeof phase3_oper);
	Q[qf++] = ((pack_phase_num_3(PHASE3S) / 70ll) << 16) + pack_phase_num_3(PHASE3S) % 70;
	x = pack_phase_num_3(PHASE3S);
	phase3_oper[x] = 0;
	while (qf > qr) {//도착지의 상태 96가지를 bfs로 찾는다.
		x = Q[qr++];
		xpk = (x >> 16) * 70 + (x & 65535);
		phase3_is_end[xpk] = 1;
		phase4_corner_perm[qr - 1] = phase3_corner_perm[x >> 16];//페이즈4의 가능한 경우들을 미리 찾아둔다.
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 3; j++) {
				y = x;
				corner = (y >> 16);
				edge = (y * 65535);
				corner = phase3_corner[corner][i][j];
				edge = phase3_edge[edge][i][j];
				y = (corner << 16) + edge;
				if (j != 1) continue;//짝수번만 돌린다
				ypk = (y >> 16) * 70 + (y & 65535);
				assert(ypk >= 0 && ypk < 2822400);
				if (phase3_oper[ypk] < 0) {
					phase3_oper[ypk] = 0;
					phase3_last[ypk] = 0;
					Q[qf++] = y;
				}
			}
		}
	}
	std::sort(phase4_corner_perm, phase4_corner_perm + qr);
	qr = 0;
	while (qf < qr) {
		x = Q[qr++];
		xpk = (x >> 16) * 70 + (x & 65535);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 3; j++) {
				if (j != 1 && (i != 2 && i != 4)) continue;//양 옆을 제외하면 짝수번만 돌릴 수 있음
				y = x;
				corner = (y >> 16);
				edge = (y & 65535);
				corner = phase3_corner[corner][i][j];
				edge = phase3_edge[edge][i][j];
				y = (corner << 16) + edge;

				ypk = (y >> 16) * 70 + (y & 65535);
				assert(ypk >= 0 && ypk < 2822400);
				if (phase3_oper[ypk] < 0) {
					phase3_oper[ypk] = i * 10 + (3 - j);
					phase3_last[ypk] = xpk;
					Q[qf++] = y;
				}
			}
		}
	}
	return;
}
/* PHASE 3 */


/* PHASE 4 */
void prec_phase4_edge_perm(ll ep, int cnt) {
	if (cnt < 0) {
		phase4_edge_perm[phase4_edge_perm_cnt++] = ep;
		return;
	}
	for (ll i = 0; i < 12; i++) {
		if (phase4_edge_perm_chk[i]) continue;
		if (!oper4_edge[i][cnt]) continue;
		phase4_edge_perm_chk[i] = 1;
		prec_phase4_edge_perm(ep + (i << (cnt << 2)), cnt - 1);
		phase4_edge_perm_chk[i] = 0;
	}
	return;
}
//8진법 형태의 코너 순열을 이분 탐색으로 찾는 함수
int corner_to_perm_4(int c) {
	int s = 0, e = 95, m;
	while (s < e) {
		m = s + e >> 1;
		if (phase4_corner_perm[m] >= c) e = m;
		else s = m + 1;
	}
	assert(phase4_corner_perm[s] == c);
	return s;
}
//8진법 형태의 엣지 순열을 이분 탐색으로 찾는 함수
int edge_to_perm_4(ll ed) {
	int s = 0, e = 13823, m;
	while (s < e) {
		m = s + e >> 1;
		if (phase4_edge_perm[m] >= ed) e = m;
		else s = m + 1;
	}
	assert(phase4_edge_perm[s] == ed);
	return s;
}
int get_phase_num_4() {
	ll ret = 0, corn_num = 0, edge_num = 0, tedge;
	for (int i = 0; i < 8; i++) {
		corn_num += ((ll)get_corner_num_3(i) << (i * 3ll));
	}
	ret = corner_to_perm_4(corn_num);
	ret *= 13824;
	for (ll i = 0; i < 12; i++) {
		tedge = color_edge[S[edge_od[i][0]] * 8 + S[edge_od[i][1]]];
		assert(tedge >= 0);
		edge_num += (tedge << (i << 2));
	}
	ret += edge_to_perm_4(edge_num);
	return ret;
}
ll phase4_edge_move(ll edge_num, int f) {
	ll store = ((edge_num >> (phase2_edge_oper[f][3] << 2)) & 15ll);
	for (int i = 3; i >= 1; i--) {
		edge_num -= (edge_num & (15ll << (phase2_edge_oper[f][i] << 2)));
		edge_num += (((edge_num >> (phase2_edge_oper[f][i - 1] << 2)) & 15ll)
			<< (phase2_edge_oper[f][i] << 2));
	}
	edge_num -= (edge_num & (15ll << (phase2_edge_oper[f][0] << 2)));
	edge_num += (store << (phase2_edge_oper[f][0] << 2));
	return edge_num;
}
//페이즈4의 코너 전처리
void prec_phase4_corner() {
	ll qf = 0, qr = 0, x = PHASE4C, xpk, y, ypk;
	memset(phase4_edge, -1, sizeof phase4_edge);
	Q[qf++] = x;
	while (qf > qr) {
		x = Q[qr++];
		xpk = corner_to_perm_4(x);
		for (int i = 0; i < 6; i++) {
			y = x;
			for (int j = 0; j < 3; j++) {
				y = phase3_corner_move(y, i);
				if (j != 1) continue;//짝수번만 돌릴 수 있음
				ypk = corner_to_perm_4(y);
				assert(ypk >= 0 && ypk < 96);
				phase4_corner[xpk][i][j] = ypk;
				if (!phase4_corner_vis[ypk]) {
					Q[qf++] = y;
					phase4_corner_vis[ypk] = 1;
				}
			}
		}
	}
	return;
}
void prec_phase4_edge() {
	ll qf = 0, qr = 0, x = PHASE4E, xpk, y, ypk;
	memset(phase4_edge, -1, sizeof phase4_edge);
	Q[qf++] = x;
	while (qf > qr) {
		x = Q[qr++];
		xpk = edge_to_perm_4(x);
		for (int i = 0; i < 6; i++) {
			y = x;
			for (int j = 0; j < 2; j++) {
				y = phase4_edge_move(y, i);
				if (j != 1) continue;//짝수번만 돌릴 수 있음
				ypk = edge_to_perm_4(y);
				assert(ypk >= 0 && ypk < 13824);
				phase4_edge[xpk][i][j] = ypk;
				if (!phase4_edge_vis[ypk]) {
					Q[qf++] = y;
					phase4_edge_vis[ypk] = 1;
				}
			}
		}
	}
	return;
}
void prec_phase4() {
	ll qf = 0, qr = 0, x, xpk, y, ypk, corner, edge;
	prec_phase4_edge_perm(0, 11);
	prec_phase4_corner();
	prec_phase4_edge();
	memset(phase4_oper, -1, sizeof phase4_oper);
	x = (corner_to_perm_4(PHASE4C) << 16) + edge_to_perm_4(PHASE4E);
	Q[qf++] = x;
	phase4_oper[(corner_to_perm_4(PHASE4C) * 13824) + edge_to_perm_4(PHASE4E)] = 0;
	while (qf > qr) {
		x = Q[qr++];
		xpk = (x >> 16) * 13824 + (x & 65535);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 2; j++) {
				if (j != 1) continue;//짝수번만 돌릴 수 있음
				y = x;
				corner = (y >> 16);
				edge = (y & 65535);
				corner = phase4_corner[corner][i][j];
				edge = phase4_edge[edge][i][j];
				assert(corner >= 0 && corner < 96 && edge >= 0 && edge < 13824);
				y = (corner << 16) + edge;
				ypk = (y >> 16) * 13824 + (y & 65535);
				assert(ypk >= 0 && ypk < 1327104);
				if (phase4_oper[ypk] < 0) {
					phase4_oper[ypk] = i + 10 + (3 - j);
					phase4_last[ypk] = xpk;
					Q[qf++] = y;
				}
			}
		}
	}
	return;
}
/* PHASE 4 */


void solve() {
	std::cin.tie(0)->sync_with_stdio(0);
	std::cout.tie(0);
	int edge_parity, phase2_num, phase3_num, phase4_num, reslen, last, tmp;
	ll move_cnt = 0;//전체 이동 횟수 합
	srand(time(NULL));
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) {
			if (j == 0) comb[i][j] = 1;
			else if (i == 0) comb[i][j] = 0;
			else comb[i][j] = comb[i - 1][j] + comb[i - 1][j - 1];
		}
	}
	prec_phase1();
	prec_phase2();
	prec_phase3();
	prec_phase4();
	std::cin >> N;
	while (N--) {
		for (int i = 0; i < 54; i++) std::cin >> S[tile_od[i]];
		reslen = 0;
		/*
		* 페이즈 1
		* 각 모서리 블록의 방향을 바꾼다
		* 위아래 면을 짝수번 돌려서 원래 있어야 하는 상태로 옮길 수 있는 경우로 만든다
		* ZZ공식이랑 비슷한 느낌
		* 사용하는 동작: U, D, F, F', F2, R, R', R2, L, L',L2 , B, B', B2
		* 경우의 수 2,048
		*/
		edge_parity = check_edge_parity();
		while (edge_parity > 0) {
			last = edge_parity;
			res[reslen++] = phase1_oper[edge_parity];
			edge_parity = phase1_last[edge_parity];
			move(res[reslen - 1]);
			tmp = check_edge_parity();
			if (tmp != edge_parity) {
				std::cout << "FUCK:: 1::\n";
			}
		}
		edge_parity = check_edge_parity();
		assert(!edge_parity);//모든 엣지 조각을 원하는 방향으로 이동

		/*
		* 페이즈 2
		* 각 꼭짓점의 방향을 바꾸는 동시에
		* FU, FD, BU, BD 엣지의 위치를 옮김
		* 사용하는 동작: U2, D2, F, F', F2, R, R', R2, L, L',L2 , B, B', B2
		* 경우의 수: 1,082,565가지 (corner 2,187, edge 495)
		*/
		phase2_num = get_phase_num_2();
		while (phase2_num != 54) {
			last = phase2_num;
			res[reslen++] = phase2_oper[phase2_num];
			phase2_num = phase2_last[phase2_num];
			move(res[reslen - 1]);
			tmp = get_phase_num_2();
			if (tmp != phase2_num) {
				std::cout << "FUCK:: 2::\n";
			}
		}
		phase2_num = get_phase_num_2();
		assert(phase2_num == 54);//모든 코거 조각이 원하는 방향, 가운데 면 엣지 조각을 가운데로 모았음

		/*
		* 페이즈 3
		* 모든 칸의 색상을 원래 색 또는 그 반대 위치 색으로 바꿈
		* 사용하는 동작: U2, D2, F2, B2, R, R', R2, L, L', L2
		* 경우의 수: 29,400 * 96 (corner 420 * 96, edge 70)
		*/
		phase3_num = get_phase_num_3();
		while (!phase3_is_end[phase3_num]) {
			last = phase3_num;
			res[reslen++] = phase3_oper[phase3_num];
			phase3_num = phase3_last[phase3_num];
			move(res[reslen - 1]);
			tmp = get_phase_num_3();
			if (tmp != phase3_num) {
				std::cout << "FUCK:: 3::\n";
			}
		}
		phase3_num = get_phase_num_3();
		assert(phase3_is_end[phase3_num]);//모든 면의 색이 원하는 위치 혹은 반대 면으로 이동
		
		/*
		* 페이즈 4
		* 반바퀴 회전만을 사용해 큐브를 맞춤
		* 사용하는 동작: 모든 면 2회 회전
		* 경우의 수: 663552 (corner 96 * edge 6912)
		*/
		phase4_num = get_phase_num_4();
		while (phase4_num != 1327103) {
			last = phase4_num;
			res[reslen++] = phase4_oper[phase4_num];
			phase4_num = phase4_last[phase4_num];
			move(res[reslen - 1]);
			tmp = get_phase_num_4();
			if (tmp != phase4_num) {
				std::cout << "FUCK:: 4::\n";
			}
		}
		assert(check());//complete

		print_sol(reslen);
		//move_cnt += reslen;
		//len_count[reslen]++;
	}
	return;

}
int main() { solve(); return 0; }//boj24902 Fewest Moves Challenge