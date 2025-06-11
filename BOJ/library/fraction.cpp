#include <iostream>
#include <numeric>
#include <stdexcept>
#include <cassert>
typedef long long ll;
ll gcd(ll a, ll b) { while (b) { a %= b; std::swap(a, b); }return a; }

struct Frac {
	ll num, den;//numerator, denominator (num/den)
	Frac(ll n = 0, ll d = 1) : num(n), den(d) {
		assert(den);
		if (den < 0) num *= -1, den *= -1;
		simplify();
	}
	void fit() {
		if (num < 0) num = 1;
		else if (num > den) num = den = 1;
		return;
	}
	void simplify() {
		if (!num) { den = 1; return; }
		ll com = gcd(std::abs(num), std::abs(den));
		num /= com; den /= com;
		return;
	}
	friend std::ostream& operator << (std::ostream& os, const Frac& f) {
		os << f.num; if (f.den != 1) os << "/" << f.den; return os;
	}
	bool operator < (const Frac& o) const { return num * o.den < o.num * den; }
};
/*
struct Frac {
	ll num, den;//numerator, denominator

	Fraction() : numerator(0), denominator(1) {}


	Fraction(long long n) : numerator(n), denominator(1) {
		simplify();
	}

	// 분자와 분모를 받는 생성자
	Fraction(long long num, long long den) : numerator(num), denominator(den) {
		if (den == 0) {
			throw std::invalid_argument("Denominator cannot be zero.");
		}
		if (den < 0) { // 분모를 항상 양수로 유지
			numerator = -numerator;
			denominator = -denominator;
		}

		// 0부터 1까지의 범위 체크
		// 문제의 요구사항에 따라 0 <= 분자 <= 분모 를 만족해야 합니다.
		if (numerator < 0 || numerator > denominator) {
			throw std::out_of_range("Fraction out of [0, 1] range.");
		}

		simplify(); // 항상 기약분수 형태로 유지
	}

	// 분수를 기약분수 형태로 단순화하는 private 메서드
	void simplify() {
		if (numerator == 0) { // 0/X 는 0/1 로
			denominator = 1;
			return;
		}
		// std::gcd를 사용하거나 custom_gcd를 사용합니다.
		// ICPC 환경에서는 C++17이 아닐 수도 있으므로 custom_gcd를 선호할 수 있습니다.
		long long common = custom_gcd(std::abs(numerator), denominator);
		numerator /= common;
		denominator /= common;
	}

	// 출력 스트림 오버로딩 (cout << Fraction 객체)
	friend std::ostream& operator<<(std::ostream& os, const Fraction& f) {
		os << f.numerator;
		if (f.denominator != 1) { // 분모가 1이 아니면 분모도 출력
			os << "/" << f.denominator;
		}
		return os;
	}

	bool operator<(const Fraction& other) const {
		return numerator * other.denominator < other.numerator * denominator;
	}
	bool operator<=(const Fraction& other) const {
		return numerator * other.denominator <= other.numerator * denominator;
	}
	bool operator>(const Fraction& other) const {
		return numerator * other.denominator > other.numerator * denominator;
	}
	bool operator>=(const Fraction& other) const {
		return numerator * other.denominator >= other.numerator * denominator;
	}
	bool operator==(const Fraction& other) const {
		return numerator == other.numerator && denominator == other.denominator;
	}
	bool operator!=(const Fraction& other) const { return !(*this == other); }
};
*/