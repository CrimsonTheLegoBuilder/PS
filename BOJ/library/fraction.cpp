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

	// ���ڿ� �и� �޴� ������
	Fraction(long long num, long long den) : numerator(num), denominator(den) {
		if (den == 0) {
			throw std::invalid_argument("Denominator cannot be zero.");
		}
		if (den < 0) { // �и� �׻� ����� ����
			numerator = -numerator;
			denominator = -denominator;
		}

		// 0���� 1������ ���� üũ
		// ������ �䱸���׿� ���� 0 <= ���� <= �и� �� �����ؾ� �մϴ�.
		if (numerator < 0 || numerator > denominator) {
			throw std::out_of_range("Fraction out of [0, 1] range.");
		}

		simplify(); // �׻� ���м� ���·� ����
	}

	// �м��� ���м� ���·� �ܼ�ȭ�ϴ� private �޼���
	void simplify() {
		if (numerator == 0) { // 0/X �� 0/1 ��
			denominator = 1;
			return;
		}
		// std::gcd�� ����ϰų� custom_gcd�� ����մϴ�.
		// ICPC ȯ�濡���� C++17�� �ƴ� ���� �����Ƿ� custom_gcd�� ��ȣ�� �� �ֽ��ϴ�.
		long long common = custom_gcd(std::abs(numerator), denominator);
		numerator /= common;
		denominator /= common;
	}

	// ��� ��Ʈ�� �����ε� (cout << Fraction ��ü)
	friend std::ostream& operator<<(std::ostream& os, const Fraction& f) {
		os << f.numerator;
		if (f.denominator != 1) { // �и� 1�� �ƴϸ� �и� ���
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