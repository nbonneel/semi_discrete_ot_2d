// Small semi-discrete optimal transport class
// used to compute the uniformity of a point set by 
// computing the optimal transport cost to a uniform
// distribution on the square
// 
// by Nicolas BONNEEL, CNRS. July 2023.
// 
// Copyright (c) 2023 CNRS Université de Lyon All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met :

// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#define _CRT_SECURE_NO_WARNINGS 1


#include <sstream>
#include <filesystem>
#include <fstream>

#include "transport.h"


// testing on toy problems: random sampling and some rank1 point sets
void testOT2D() {

	// for random points: 
	// 1k points in 0.05 s
	// 4k points in 2 s
	// 16k points in 16 s
	// 65k points in 3 min
	// 262k points: 19 min

	int N = 64 * 64;
	OptimalTransport2D ot;
	ot.N = N;
	ot.V.samples.resize(N);

	std::default_random_engine engine;
	std::uniform_real_distribution<double> uniform(0, 1);

	Vector v(uniform(engine), uniform(engine));
	for (int i = 0; i < N; i++) {

		ot.V.samples[i] = Vector(uniform(engine), uniform(engine));
		//ot.V.samples[i] = Vector(v[0]*i-floor(v[0]*i), v[1] * i - floor(v[1] * i)); // crappy rank 1
	}

	/*for (int i = 0; i < 64; i++) {
		for (int j = 0; j < 64; j++) {
			ot.V.samples[i * 64 + j] = Vector(i / 64.0 + 1. / 128., j / 64.0 + 1. / 128.); // regular grid
		}
	}*/

	const auto start = std::chrono::steady_clock::now();

	double squaredOTdist = ot.optimize();

	const auto end = std::chrono::steady_clock::now();
	const std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	std::cout << squaredOTdist << std::endl;

	ot.V.build();
	ot.V.save_svg(ot.V.cells, "OTSmall.svg");

}



int main(int argc, char* argv[]) {


	testOT2D();


	return 0;
}