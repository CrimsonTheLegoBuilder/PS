#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <chrono>
#include <iomanip>
const int NUM_FILES = 50;
const int NUM_POINTS = 8000;
const int NUM_GROUPS = 140;
const int MAX_COORD = 814000;
int main() {
	std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<int> dist(0, MAX_COORD);
	for (int i = 1; i <= NUM_FILES; i++) {
		std::string base_path = "../../tests/814_3/";
		std::string filename = (i < 10 ? "0" : "") + std::to_string(i) + ".in";
		std::string full_path = base_path + filename;
		std::ofstream out(full_path);
		if (!out.is_open()) {
			std::cerr << "Error: " << full_path << "\n";
			return 1;
		}
		out << NUM_POINTS << " " << NUM_GROUPS << "\n";
		for (int j = 0; j < NUM_POINTS; j++) {
			out << dist(rng) << " " << dist(rng) << "\n";
		}
		std::cout << "Generated: " << full_path << "\n";
		out.close();
	}
	std::cout << "Create complete!!\n";
	return 0;
}