#include <iostream> // ����¿�
#include <numeric>  // std::gcd (C++17 ����, ������ ���� ����)
#include <stdexcept> // ���� ó����

// C++17 �̸� �����Ϸ��� ���� std::gcd�� ���� ���� (���� ����)
// �� �Լ��� ǥ�� ���̺귯���� <numeric> ����� �ִ� std::gcd�� �����ϰ� �����մϴ�.
long long custom_gcd(long long a, long long b) {
    while (b) {
        a %= b;
        std::swap(a, b);
    }
    return a;
}

// 0���� 1������ �м��� �����ϴ� ����ü
struct Fraction {
    long long numerator;   // ���� (numerator)
    long long denominator; // �и� (denominator)

    // �⺻ ������
    Fraction() : numerator(0), denominator(1) {}

    // ���� n�� �м� n/1�� ǥ���ϴ� ������
    Fraction(long long n) : numerator(n), denominator(1) {
        // 0���� 1������ ���� üũ�� �ܺο��� �� �� �ֵ��� ��.
        // ���� ���, Fraction f(5) �� 5/1 �� ������, �� ����ü ��ü�� ������ �������� ����.
        // ������ ���� ���ǿ��� "0���� 1����"��� �����Ƿ�, �׷��� ���� ���̶�� �����մϴ�.
        // ������ ������ �����ϰ� �ʹٸ� �Ʒ� �ּ� ó���� �κ��� �߰��� �� �ֽ��ϴ�.
        // if (numerator < 0 || numerator > denominator) {
        //     throw std::out_of_range("Fraction out of [0, 1] range.");
        // }
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

    // --- ��� �� ������ �����ε� ---
    // �� �м��� ���� �� ũ�ν� ���� ����մϴ�.
    // a/b vs c/d  =>  ad vs bc
    // �����÷ο� ������ ���� long long�� ����ؾ� �մϴ�.

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
        // �̹� simplify()�� ���� ���м��� �����ǹǷ�,
        // ����, �и� ������ ���� �м��Դϴ�.
        return numerator == other.numerator && denominator == other.denominator;
    }

    bool operator!=(const Fraction& other) const {
        return !(*this == other);
    }
};

// --- ��� ���� ---
int main() {
    // 0���� 1 ������ �м� ����
    Fraction f1(1, 2); // 1/2
    Fraction f2(3, 4); // 3/4
    Fraction f3(1, 3); // 1/3
    Fraction f4(2, 4); // 1/2 (�ڵ����� �ܼ�ȭ��)
    Fraction f5(0, 5); // 0/1
    Fraction f6(1);    // 1/1
    Fraction f7(0);    // 0/1

    std::cout << "f1: " << f1 << std::endl;
    std::cout << "f2: " << f2 << std::endl;
    std::cout << "f3: " << f3 << std::endl;
    std::cout << "f4: " << f4 << " (simplified from 2/4)" << std::endl;
    std::cout << "f5: " << f5 << " (simplified from 0/5)" << std::endl;
    std::cout << "f6: " << f6 << std::endl;
    std::cout << "f7: " << f7 << std::endl;

    // ��� ��
    std::cout << "\n--- Comparisons ---" << std::endl;
    std::cout << "f1 (" << f1 << ") < f2 (" << f2 << "): " << (f1 < f2 ? "true" : "false") << std::endl; // true (0.5 < 0.75)
    std::cout << "f1 (" << f1 << ") == f4 (" << f4 << "): " << (f1 == f4 ? "true" : "false") << std::endl; // true (0.5 == 0.5)
    std::cout << "f2 (" << f2 << ") > f3 (" << f3 << "): " << (f2 > f3 ? "true" : "false") << std::endl; // true (0.75 > 0.33)
    std::cout << "f5 (" << f5 << ") <= f1 (" << f1 << "): " << (f5 <= f1 ? "true" : "false") << std::endl; // true (0 <= 0.5)
    std::cout << "f6 (" << f6 << ") >= f2 (" << f2 << "): " << (f6 >= f2 ? "true" : "false") << std::endl; // true (1 >= 0.75)
    std::cout << "f3 (" << f3 << ") != f1 (" << f1 << "): " << (f3 != f1 ? "true" : "false") << std::endl; // true

    // ���� ����� ��� (���� ó�� Ȯ��)
    std::cout << "\n--- Range Check ---" << std::endl;
    try {
        Fraction f_invalid1(5, 3); // 5/3 ( > 1)
        std::cout << "Invalid fraction: " << f_invalid1 << std::endl;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
    }

    try {
        Fraction f_invalid2(-1, 2); // -1/2 ( < 0)
        std::cout << "Invalid fraction: " << f_invalid2 << std::endl;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
    }

    try {
        Fraction f_invalid3(1, 0); // Denominator is zero
        std::cout << "Invalid fraction: " << f_invalid3 << std::endl;
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}