#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <vector>
typedef long double ld;
const ld TOL = 1e-4;

std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

bool doubleCompare(const ld& o, const ld& a, const ld& tol = TOL) {
    if (std::isnan(o) || std::isnan(a)) return false;
    if (std::isinf(o) || std::isinf(a)) return false;
    if (o == a) return true;
    ld diff = std::fabs(o - a);
    if (diff <= tol) return 1;
    ld minv = std::min(o * (1.0 - tol), o * (1.0 + tol));
    ld maxv = std::max(o * (1.0 - tol), o * (1.0 + tol));
    return minv <= a && a <= maxv;
}

bool OK[1000];

#define TEST 95

std::string Dir = "../tests/candle_renew/";
int main() {
    for (int i = 1; i <= TEST; i++) {
        std::string Qout = Dir + "out/" + (i < 10 ? "0" : "") + std::to_string(i) + ".out";
        std::string Qans = Dir + "out_hpi/" + (i < 10 ? "0" : "") + std::to_string(i) + ".out";
        std::ifstream out(Qout);
        std::ifstream ans(Qans);

        if (!ans || !out) {
            std::cerr << "Error: cannot open file(s)\n";
            return 1;
        }

        std::vector<ld> A, O;
        ld val;
        while (ans >> val) A.push_back(val);
        while (out >> val) O.push_back(val);

        if (A.size() != O.size()) {
            std::cerr << "test " << i << " Mismatch: expected " << A.size()
                << " values, got " << O.size() << "\n";
            return 1;
        }

        int sz = A.size();
        bool f = 1;
        for (int v = 0; v < sz; v++) {
            if (!doubleCompare(A[v], O[v])) {
                f = 0;
                std::cout << "v:: " << v << "\n";
                break;
            }
        }
        if (!f) continue;

        OK[i] = 1;
    }
    for (int i = 1; i <= TEST; i++) {
        if (OK[i]) std::cout << "test case #" << i << " Accepted.\n";
        else std::cout << "test case #" << i << " Fxxk.\n";
    }
    return 0;
}//kitpc? 14? candle & candle & candle & shadow
/*
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>
typedef long double ld;
const ld TOL = 1e-7;

std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

bool doubleCompare(const ld& o, const ld& a, const ld& tol = TOL) {
    if (std::isnan(o) || std::isnan(a)) return false;
    if (std::isinf(o) || std::isinf(a)) return false;
    if (o == a) return true;
    ld diff = std::fabs(o - a);
    if (diff <= tol) return 1;
    ld minv = std::min(o * (1.0 - tol), o * (1.0 + tol));
    ld maxv = std::max(o * (1.0 - tol), o * (1.0 + tol));
    return minv <= a && a <= maxv;
}

bool OK[1000];

#define TEST 51

std::string Dir = "../tests/candle/";
int main() {
	for (int i = 1; i <= TEST; i++) {
		std::string Qout = Dir + "ans0/" + (i < 10 ? "0" : "") + std::to_string(i) + ".out";
		std::string Qans = Dir + "out10/" + (i < 10 ? "0" : "") + std::to_string(i) + ".out";
		std::ifstream out(Qout);
		std::ifstream ans(Qans);
		if (!ans || !out) {
			std::cerr << "No file::\n";
			return 1;
		}
        int case_num = 0;
        std::string a_line, o_line;

        std::cout << "test #" << i << "\n";
        while (std::getline(ans, a_line)) {
            ++case_num;

            // Case header
            if (!std::getline(out, o_line)) {
                std::cerr << "Case #" << case_num << " None\n";
                return 1;
            }

            a_line = trim(a_line);
            o_line = trim(o_line);
            if (a_line != o_line) {
                std::cerr << "Case #" << case_num << " Header error - expected: '" << a_line << "', found: '" << o_line << "'\n";
                return 1;
            }

            for (int i = 0; i < 8; ++i) {
                if (!std::getline(ans, a_line) || !std::getline(out, o_line)) {
                    std::cerr << "Case #" << case_num << ", line " << i + 1 << " None\n";
                    return 1;
                }

                a_line = trim(a_line);
                o_line = trim(o_line);

                size_t acol = a_line.find(':');
                size_t ocol = o_line.find(':');

                if (acol == std::string::npos || ocol == std::string::npos) {
                    std::cerr << "Case #" << case_num << ", line " << i + 1 << " error - ':' None\n";
                    return 1;
                }

                std::string alabel = trim(a_line.substr(0, acol));
                std::string olabel = trim(o_line.substr(0, ocol));

                if (alabel != olabel) {
                    std::cerr << "Case #" << case_num << ", label wrong - expected: '" << alabel << "', found: '" << olabel << "'\n";
                    return 1;
                }

                ld aval = std::stod(trim(a_line.substr(acol + 1)));
                ld oval = std::stod(trim(o_line.substr(ocol + 1)));

                ld jmax = std::max(aval * (1 + TOL + 1e-15), aval * (1 - TOL - 1e-15));
                ld jmin = std::min(aval * (1 + TOL + 1e-15), aval * (1 - TOL - 1e-15));

                if (!doubleCompare(aval, oval)) {
                    std::cerr << "Case #" << case_num << ", " << alabel
                        << " Not equal - expected: " << aval
                        << ", found: " << oval
                        << ", abs_error: " << std::fabs(aval - oval)
                        << ", tolerance: [" << jmin << ", " << jmax << "]\n";
                    return 1;
                }
            }
        }

        if (std::getline(out, o_line)) {
            std::cerr << "OUT FUCK: -> '" << o_line << "'\n";
            return 1;
        }

        OK[i] = 1;
	}
    for (int i = 1; i <= TEST; i++) {
        if (OK[i]) std::cout << "test case #" << i << " Accepted.\n";
        else std::cout << "test case #" << i << " Fxxk.\n";
    }
	return 0;
}//kitpc? 14? candle & candle & candle & shadow
*/