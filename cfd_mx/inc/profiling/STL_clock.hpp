#pragma once

#include <string>			// for C++ string
#include <vector>			// for vector container
#include <chrono>
class StopWatch {

	std::chrono::high_resolution_clock::time_point t_start, t_stop;

	double total=0;

public:
	void start() {
		t_start = std::chrono::high_resolution_clock::now();
	}

	void stop() {
		t_stop = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> d = t_stop - t_start;
		total += d.count();
	}

	double elapsedTime() {
		return total;
	}

	void init_start() { 
		total = 0.0; 
		t_start = std::chrono::high_resolution_clock::now();
	}

	static double resolution() {
		auto tmp = std::chrono::high_resolution_clock::period();
		return (double)tmp.num / tmp.den;
	}
};
