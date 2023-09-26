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

using namespace transport;

// testing on toy problems: random sampling and some rank1 point sets
void test_optimal_transport_2d() {

	// for random points: 
	// 1k points in 0.03 s  (1 thread)
	// 4k points in 0.6 s (1 thread)
	// 16k points in 2.4 s
	// 65k points in 6.6 s
	// 262k points: 1 min 20 s

	// 1M points : 863 s

	int N = 64 * 64*16;  //number of points
	OptimalTransport2D ot;


	omp_set_num_threads(16); // play with the number of threads on your machine. Many threads will hurt for low sample count

	std::default_random_engine engine;
	std::uniform_real_distribution<double> uniform(0, 1);


	ot.V.vertices.resize(N);
	for (int i = 0; i < N; i++) {
		ot.V.vertices[i] = Vector(uniform(engine), uniform(engine));
	}

	const auto start = std::chrono::steady_clock::now();

	double squaredOTdist = ot.optimize(100); //100 Newton iterations max.

	const auto end = std::chrono::steady_clock::now();
	const std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	std::cout << "OT squared distance: "<<squaredOTdist << std::endl;
	std::cout << "worst relative Voronoi cell area error: " << ot.V.worst_cell_area_relative_diff_percent() << "%"<<std::endl;
	std::cout << "average relative Voronoi cell area error: " << ot.V.avg_cell_area_relative_diff_percent() << "%"<<std::endl;


	ot.V.save_svg("OT_result.svg", false, true, true);

}


// computes the weighted (or also unweighted) delaunay and its power diagram quite efficiently
void test_power_diagram() {

	// for random points with random weights of the form 1+uniform(0,1) (delaunay + dual): 
	// 16k points: 0.007 s
	// 64k points: 0.029 s
	// 256k points: 0.12 s
	// 1M points: 0.55 s
	// 4M points: 2.7 s
	// 8M points: 6.2 s

	// for unweighted points (slower because more polygons to construct : all vertices have one cell):
	// 256k points: 1s
	// 1M points: 4 s

	std::default_random_engine engine;
	std::uniform_real_distribution<double> uniform(0, 1);

	engine.seed(1385);

	Bowyer2D vd;
	for (int i = 0; i < 1024 * 16; i++) {
		vd.vertices.push_back(Vector(uniform(engine), uniform(engine)));
		vd.weights.push_back(1+uniform(engine));
	}



	const auto start = std::chrono::steady_clock::now();
	vd.compute_delaunay();
	vd.compute_dual();

	const auto end = std::chrono::steady_clock::now();
	const std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


	vd.save_svg("voronoi_delaunay.svg", true, true, true, false);
}

void test_OT_lloyd() {
	std::default_random_engine engine;
	std::uniform_real_distribution<double> uniform(0, 1);

	engine.seed(1385);

	OptimalTransport2D ot;

	for (int i = 0; i < 1024; i++) {
		ot.V.vertices.push_back(Vector(uniform(engine), uniform(engine)));
	}

	const auto start = std::chrono::steady_clock::now();
	ot.ot_lloyd();
	const auto end = std::chrono::steady_clock::now();
	const std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


	ot.V.save_svg("lloyd.svg", false, true, true);
}




int main(int argc, char* argv[]) {


	test_optimal_transport_2d();

	//test_power_diagram();
    //test_OT_lloyd();


	return 0;
}