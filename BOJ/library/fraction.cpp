#include <iostream> // 입출력용
#include <numeric>  // std::gcd (C++17 부터, 없으면 직접 구현)
#include <stdexcept> // 예외 처리용

// C++17 미만 컴파일러를 위해 std::gcd를 직접 구현 (선택 사항)
// 이 함수는 표준 라이브러리의 <numeric> 헤더에 있는 std::gcd와 동일하게 동작합니다.
long long custom_gcd(long long a, long long b) {
    while (b) {
        a %= b;
        std::swap(a, b);
    }
    return a;
}

// 0부터 1까지의 분수를 관리하는 구조체
struct Fraction {
    long long numerator;   // 분자 (numerator)
    long long denominator; // 분모 (denominator)

    // 기본 생성자
    Fraction() : numerator(0), denominator(1) {}

    // 정수 n을 분수 n/1로 표현하는 생성자
    Fraction(long long n) : numerator(n), denominator(1) {
        // 0부터 1까지의 범위 체크는 외부에서 할 수 있도록 함.
        // 예를 들어, Fraction f(5) 는 5/1 이 되지만, 이 구조체 자체는 범위를 강제하지 않음.
        // 하지만 문제 조건에서 "0부터 1까지"라고 했으므로, 그렇게 사용될 것이라고 가정합니다.
        // 생성자 내에서 강제하고 싶다면 아래 주석 처리된 부분을 추가할 수 있습니다.
        // if (numerator < 0 || numerator > denominator) {
        //     throw std::out_of_range("Fraction out of [0, 1] range.");
        // }
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

    // --- 대소 비교 연산자 오버로딩 ---
    // 두 분수를 비교할 때 크로스 곱을 사용합니다.
    // a/b vs c/d  =>  ad vs bc
    // 오버플로우 방지를 위해 long long을 사용해야 합니다.

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
        // 이미 simplify()를 통해 기약분수로 관리되므로,
        // 분자, 분모가 같으면 같은 분수입니다.
        return numerator == other.numerator && denominator == other.denominator;
    }

    bool operator!=(const Fraction& other) const {
        return !(*this == other);
    }
};

// --- 사용 예시 ---
int main() {
    // 0부터 1 사이의 분수 생성
    Fraction f1(1, 2); // 1/2
    Fraction f2(3, 4); // 3/4
    Fraction f3(1, 3); // 1/3
    Fraction f4(2, 4); // 1/2 (자동으로 단순화됨)
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

    // 대소 비교
    std::cout << "\n--- Comparisons ---" << std::endl;
    std::cout << "f1 (" << f1 << ") < f2 (" << f2 << "): " << (f1 < f2 ? "true" : "false") << std::endl; // true (0.5 < 0.75)
    std::cout << "f1 (" << f1 << ") == f4 (" << f4 << "): " << (f1 == f4 ? "true" : "false") << std::endl; // true (0.5 == 0.5)
    std::cout << "f2 (" << f2 << ") > f3 (" << f3 << "): " << (f2 > f3 ? "true" : "false") << std::endl; // true (0.75 > 0.33)
    std::cout << "f5 (" << f5 << ") <= f1 (" << f1 << "): " << (f5 <= f1 ? "true" : "false") << std::endl; // true (0 <= 0.5)
    std::cout << "f6 (" << f6 << ") >= f2 (" << f2 << "): " << (f6 >= f2 ? "true" : "false") << std::endl; // true (1 >= 0.75)
    std::cout << "f3 (" << f3 << ") != f1 (" << f1 << "): " << (f3 != f1 ? "true" : "false") << std::endl; // true

    // 범위 벗어나는 경우 (예외 처리 확인)
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