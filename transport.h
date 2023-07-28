#include <vector>
#include <random>
#include <omp.h>

#include <chrono>

#include "picoflann.h"



// a 3d Vector class ; 2d points are embedded in 3d for kNN search with power distance
class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) { coords[0] = x; coords[1] = y; coords[2] = z; }
	double operator[](int i) const {
		return coords[i];
	}
	double& operator[](int i) {
		return coords[i];
	}
	Vector& operator+=(const Vector& V) {
		coords[0] += V[0];
		coords[1] += V[1];
		coords[2] += V[2];
		return *this;
	}
	Vector& operator/=(double s) {
		coords[0] /= s;
		coords[1] /= s;
		coords[2] /= s;
		return *this;
	}
	double getNorm2_3d() const {
		return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
	}
	double getNorm2_2d() const {
		return coords[0] * coords[0] + coords[1] * coords[1];
	}
	double getNorm_2d() const {
		return sqrt(coords[0] * coords[0] + coords[1] * coords[1]);
	}
	void normalize_3d() {
		double n = sqrt(getNorm2_3d());
		coords[0] /= n;
		coords[1] /= n;
		coords[2] /= n;
	}

	double coords[3];
};


bool operator==(const Vector& a, const Vector& b);
Vector operator*(double a, const Vector& B);
Vector operator*(const Vector& B, double a);
Vector operator*(const Vector& A, const Vector& B);
Vector operator/(const Vector& B, double a);
Vector operator+(const Vector& A, const Vector& B);
Vector operator-(const Vector& A, const Vector& B);
Vector operator-(const Vector& A);
double dot_3d(const Vector& A, const Vector& B);
double dot_2d(const Vector& A, const Vector& B);

template<typename T> T sqr(T x) {
	return x * x;
}



class Facet {
public:
	Facet() {}

	std::vector<Vector> vertices;
	std::vector<int> neighbors;

	// https://www.seas.upenn.edu/~sys502/extra_materials/Polygon%20Area%20and%20Centroid.pdf

	double area() const {
		// https://wwwf.imperial.ac.uk/~rn/centroid.pdf
		double a = 0;
		size_t Nv = vertices.size();
		for (size_t i = 0; i < Nv; i++) {
			a += vertices[i][0] * vertices[(i + 1) % Nv][1] - vertices[(i + 1) % Nv][0] * vertices[i][1];
		}
		return abs(a / 2);
	}

	Vector centroid() const {
		double A = 0;
		Vector C(0, 0, 0);
		size_t Nv = vertices.size();
		if (Nv == 0) return Vector(0., 0., 0.);
		if (Nv == 1) return vertices[0];
		if (Nv == 2) return (vertices[0] + vertices[1]) * 0.5;
		Vector avg(0., 0., 0.);
		for (size_t i = 0; i < Nv; i++) {
			double cc = vertices[i][0] * vertices[(i + 1) % Nv][1] - vertices[(i + 1) % Nv][0] * vertices[i][1];
			A += cc;
			avg += vertices[i];
			for (int k = 0; k < 2; k++) {
				C[k] += cc * (vertices[i][k] + vertices[(i + 1) % Nv][k]);
			}
		}
		A *= 0.5;
		if (A == 0) return avg / (double)Nv;
		Vector c = C / (6 * A);
		return c;
	}

	double integrateSquaredDistance2D(const Vector& Pi) { //https://people.sc.fsu.edu/~jburkardt/cpp_src/polygon_integrals/polygon_integrals.html
		double s = 0;
		size_t Nv = vertices.size();
		double npi2 = Pi.getNorm2_2d();
		size_t im = Nv - 1;
		for (size_t i = 0; i < Nv; i++) {
			double sm = vertices[im].getNorm2_2d() + vertices[i].getNorm2_2d() + dot_2d(vertices[im], vertices[i]) - 4 * dot_2d(Pi, vertices[im] + vertices[i]) + 6 * npi2;
			s += (vertices[im][0] * vertices[i][1] - vertices[i][0] * vertices[im][1]) * sm;
			im = i;
		}
		double v1 = abs(s) / 12.;
		return v1;
	}



	void save_svg(std::string filename, std::string fillcol = "none", bool header = true, bool footer = true) const {
		FILE* f = fopen(filename.c_str(), "a+");
		if (header)
			fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
		fprintf(f, "<g>\n");
		fprintf(f, "<polygon points = \"");
		for (int j = 0; j < vertices.size(); j++) {
			fprintf(f, "%u, %u ", (int)(vertices[j][0] * 1000), (int)(1000 - vertices[j][1] * 1000));
		}
		fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
		fprintf(f, "</g>\n");
		if (footer)
			fprintf(f, "</svg>\n");
		fclose(f);
	}
};




struct PicoFlann_Array3d_Adapter {
	inline   double operator( )(const Vector& elem, int dim)const { return elem.coords[dim]; }
};
struct PicoFlann_Array3d_Container {
	const Vector* _array;
	size_t _size;
	PicoFlann_Array3d_Container(double* array, size_t Size) :_array((Vector*)array), _size(Size) {}
	inline size_t size()const { return _size; }
	inline const Vector& at(int idx)const { return _array[idx]; }
};

class VoronoiDiagram2D {
public:
	VoronoiDiagram2D() {	
		square.vertices.resize(4);
		square.neighbors.resize(4, -1);
		square.vertices[0] = Vector(0., 0., 0.);
		square.vertices[3] = Vector(1., 0., 0.);
		square.vertices[2] = Vector(1., 1., 0.);
		square.vertices[1] = Vector(0., 1., 0.);
	};

	mutable Vector tmpMemVector[128][512];
	mutable int tmpMemInt[128][512];

	Facet square;

	void clipCellByBissectorPlane(std::vector<Facet>& cells, int si, int sj) const {

		int tid = omp_get_thread_num();
		size_t Nv = cells[si].vertices.size();


		//https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
		double	weightsi = weights[si];
		double	weightsj = weights[sj];
		const Vector diff = samples[si] - samples[sj];
		double off = (weightsi - weightsj) / diff.getNorm2_2d();
		Vector M = (samples[si] + samples[sj] - off * diff) * 0.5;

		int numNewVertices = 0;
		Vector* newVertices = &tmpMemVector[tid][0];
		int* newNeighbors = &tmpMemInt[tid][0];
		const Vector* pv2 = &cells[si].vertices[0];
		const Vector* pv1 = &cells[si].vertices[Nv - 1];


		for (size_t j = 0; j < Nv; j++, ++pv2) {

			const Vector& v1 = *pv1;
			const Vector& v2 = *pv2;

			if (dot_2d(v2 - M, diff) > -1E-12) { // v2 inside
				double d1 = dot_2d(M - v1, diff);
				if (d1 > 1E-12) { // v1 outside
					const Vector v1v2 = v2 - v1;
					double t = d1 / dot_2d(v1v2, diff);
					if (t > 1E-12 && t < 1 - 1E-12) {
						newVertices[numNewVertices][0] = v1[0] + t * v1v2[0];
						newVertices[numNewVertices][1] = v1[1] + t * v1v2[1];
						newNeighbors[numNewVertices] = sj;
						numNewVertices++;
					}
				}
				newVertices[numNewVertices] = v2;
				newNeighbors[numNewVertices] = cells[si].neighbors[j];
				numNewVertices++;
			} else {
				double d1 = dot_2d(M - v1, diff);
				if (d1 < -1E-12) {
					const Vector v1v2 = v2 - v1;
					double d2 = dot_2d(v1v2, diff);
					double t = d1 / d2;
					if (t > 1E-12 && t < 1 - 1E-12) {
						newVertices[numNewVertices][0] = v1[0] + t * v1v2[0];
						newVertices[numNewVertices][1] = v1[1] + t * v1v2[1];
						newNeighbors[numNewVertices] = cells[si].neighbors[j];
						numNewVertices++;
					}
				}
			}
			pv1 = pv2;
		}

		cells[si].vertices.resize(numNewVertices);
		cells[si].neighbors.resize(numNewVertices);
		memcpy(&cells[si].vertices[0], newVertices, numNewVertices * sizeof(Vector));
		memcpy(&cells[si].neighbors[0], newNeighbors, numNewVertices * sizeof(int));
	}

	void clear() {
		cells.clear();
	}


	std::vector<int> num_neigh;
	picoflann::KdTreeIndex<3, PicoFlann_Array3d_Adapter> kdtree;

	void build() {

		cells.resize(samples.size());

		double maxweight = -1E9;
		for (int i = 0; i < weights.size(); i++) {
			maxweight = std::max(maxweight, weights[i]);
		}
		for (int i = 0; i < samples.size(); i++) {
			samples[i].coords[2] = sqrt(maxweight - weights[std::min(i, (int)weights.size() - 1)]); // embedding to get a power diagram out of a voronoi
		}

		if (samples.size() < 1500) { // low sample count, do not bother with kNN
#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < samples.size(); i++) {
				cells[i] = square;
				for (int j = 0; j < samples.size(); j++) {
					if (i != j)
						clipCellByBissectorPlane(cells, i, j);
				}
			}
		} else {

			kdtree.build(PicoFlann_Array3d_Container(&samples[0].coords[0], samples.size()));
			PicoFlann_Array3d_Container p3container(&samples[0].coords[0], samples.size());

			if (num_neigh.size() != samples.size()) {
				num_neigh.resize(samples.size(), 64); // default number of neighbor per cell
			}

#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < samples.size(); i++) {
				double query_pt[3] = { samples[i].coords[0],  samples[i].coords[1],  samples[i].coords[2] };


				std::vector<std::pair<uint32_t, double> > resultSet;
				bool tooshort = true;

				while (tooshort) {
					cells[i] = square;

					resultSet = kdtree.searchKnn(p3container, samples[i], num_neigh[i]);

					size_t num_results = resultSet.size();

					for (size_t j = 1; j < num_results; j++) {

						clipCellByBissectorPlane(cells, i, resultSet[j].first);
						double furthest = 0;
						for (size_t k = 0; k < cells[i].vertices.size(); k++) {
							furthest = std::max(furthest, (cells[i].vertices[k] - samples[i]).getNorm2_3d());
						}

						if (resultSet[j].second > 4 * furthest) { // security radius met
							tooshort = false;
							break;
						}
					}

					if (num_results == samples.size()) break;

					if (tooshort) {
						num_results *= 2;
						num_results = std::min(num_results, samples.size());
						num_neigh[i] = (int)num_results;
					} else {
						num_neigh[i] = (int) num_results + 10;
						num_results = std::min(num_results, samples.size());
					}

				}

			}
		}

		for (int i = 0; i < samples.size(); i++) {
			samples[i].coords[2] = 0;
		}
	}


	// save an svg file showing the power diagram contained in cells
	void save_svg(const std::vector<Facet>& cells, std::string filename) {

		FILE* f = fopen(filename.c_str(), "w+");
		fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");

		fprintf(f, "<g>\n");
		for (int i = 0; i < cells.size(); i++) {
			fprintf(f, "<polygon points = \"");
			for (int j = 0; j < cells[i].vertices.size(); j++) {
				fprintf(f, "%3.3f, %3.3f ", (cells[i].vertices[j][0] * 1000), (1000 - cells[i].vertices[j][1] * 1000));
			}
			fprintf(f, "\"\nfill = \"#005BBB\" stroke = \"black\"/>\n");
		}
		// if you want to show samples, uncomment
		/*for (int i = 0; i < cells.size(); i++) {
			Vector K = cells[i].centroid();
			Vector L = samples[i];

			//fprintf(f, "<circle cx=\"%3.3f\" cy=\"%3.3f\" r=\"1\" stroke=\"black\" stroke-width=\"1\" fill = \"black\" />", (K[0] * 1000), (1000 - K[1] * 1000)); // cell centroids instead
			fprintf(f, "<circle cx=\"%3.3f\" cy=\"%3.3f\" r=\"3\" stroke=\"deeppink\" stroke-width=\"1\" fill = \"deeppink\" />", (L[0] * 1000), (1000 - L[1] * 1000));
		}*/

		fprintf(f, "</g>\n");
		fprintf(f, "</svg>\n");

		fclose(f);
	}

	std::vector<Vector> samples;
	std::vector<double> weights;
	std::vector<Facet> cells;
	size_t N;
};




class OptimalTransport2D {
public:
	OptimalTransport2D() {
		V.weights.resize(10000, 0);
	};


	// finally, more efficient to pre-store the Hessian than computing on the fly
	void precompute_hessian(std::vector<std::vector<std::pair<int, double> > >& HessianValues) {

#pragma omp parallel for
		for (int i = 0; i < N; i++) {

			double result = 0;
			double sum_contribs = 0;
			size_t num_nn = V.cells[i].neighbors.size();
			const int* nn = &V.cells[i].neighbors[0];
			const Vector* vtx = &V.cells[i].vertices[0];
			const Vector* prev_vtx = &V.cells[i].vertices[num_nn - 1];
			Vector curSample = V.samples[i];

			HessianValues[i].clear();
			HessianValues[i].reserve(V.cells[i].vertices.size());

			double sumcontrib = 0;
			for (int j = 0; j < num_nn; j++, ++nn, ++vtx) {
				double edgeLength = (*vtx - *prev_vtx).getNorm_2d();
				prev_vtx = vtx;
				if (*nn < 0) continue;
				double contrib = -edgeLength / (2. * (curSample - V.samples[*nn]).getNorm_2d());
				sumcontrib += contrib;
				HessianValues[i].push_back(std::pair<int, double>(*nn, contrib));
			}
			diagHessian[i] = -sumcontrib;

		}

	}

	double hessian_row_mult(int row_si, double* rhs) {
		double result = 0;
		size_t num_nn = HessianValues[row_si].size();
		const std::pair<int, double>* H = &HessianValues[row_si][0];
		for (size_t i = 0; i < num_nn; i++, ++H) {
			double contrib = H->second;
			result += contrib * rhs[H->first];
		}
		result += diagHessian[row_si] * rhs[row_si];
		return result;
	}

	void hessian_mult(double* rhs, double* result) {
#pragma omp parallel for 
		for (int i = 0; i < N; i++) {
			result[i] = hessian_row_mult(i, rhs);
		}
	}
	void residual(double* b, double* rhs, double* result) {
#pragma omp parallel for 
		for (int i = 0; i < N; i++) {
			result[i] = b[i] - hessian_row_mult(i, rhs);
		}
	}

	std::vector<std::vector<std::pair<int, double> > > HessianValues;
	std::vector<double> diagHessian; // also used for preconditionning
	std::vector<double> r, p, z, Ap;

	void conjugate_gradient_solve(double* rhs, double* result) { // result initialized ( = x0 )
		r.resize(N);
		p.resize(N);
		z.resize(N);
		Ap.resize(N);
		HessianValues.resize(N);
		diagHessian.resize(N);

		precompute_hessian(HessianValues);

		residual(rhs, result /* = x0 */, &r[0]);
		for (int i = 0; i < N; i++) {
			z[i] = r[i] / diagHessian[i]; // Jacobi preconditionning
		}
		p = z;

		for (int k = 0; k < std::min(N, N / 4 + 2); k++) { //had to limit the number of iterations for edgy cases where it starts to blow up

			hessian_mult(&p[0], &Ap[0]);
			double rz = 0, pAp = 0;
			for (int i = 0; i < N; i++) {
				rz += r[i] * z[i];
				pAp += p[i] * Ap[i];
			}
			double alpha = rz / pAp;
			double rzb = 0;
			double rr = 0;
			for (int i = 0; i < N; i++) {
				result[i] += alpha * p[i];
				r[i] -= alpha * Ap[i];
				z[i] = r[i] / diagHessian[i];
				rzb += r[i] * z[i];
				rr += r[i] * r[i];
			}
			if (rr < 1E-11) break;
			double beta = rzb / rz;
			for (int i = 0; i < N; i++) {
				p[i] = z[i] + beta * p[i];
			}
		}

	}

	// returns the worst relative difference in area, in percentage (1 = 1% area error w.r.t a uniform 1./N). Also used as stopping criterion
	double newton_optimize(int max_Newton_iter, double worst_area_relative_threshold_percent) {

		V.N = N;
		V.weights.resize(N);
		std::fill(V.weights.begin(), V.weights.end(), 0.0);

		double lambda = 1;
		double worstarea = 0;
		std::vector<double> grad(N), invGrad(N, 0.);
		for (int iter = 0; iter < max_Newton_iter; iter++) {
		start:
			V.clear();
			V.build();

			worstarea = 0;
			for (int i = 0; i < N; i++) {
				double a = V.cells[i].area();
				if (a <= 0) { // if an empty cell is detected, we roll back half the step size and reduce the step size by half. See Gallouet&Merigot
					//std::cout << "new lambda:"<<lambda*0.5 << std::endl;
					lambda *= 0.5;
					if (lambda < 1E-9) {
						std::cout << "ERROR during optimization ; resulting value is unreliable !" << std::endl;
						return -1;
					}
					for (int j = 0; j < N; j++) {
						V.weights[j] -= lambda * invGrad[j];
					}
					goto start; // boo me for that goto as much as you want
					break;
				}
				grad[i] = 1.0 / N - a;
				if (abs(grad[i]) > worstarea) worstarea = abs(grad[i]);
			}
			if (worstarea * N * 100 < worst_area_relative_threshold_percent) break;

			conjugate_gradient_solve(&grad[0], &invGrad[0]);

			for (int i = 0; i < N; i++) {
				V.weights[i] += lambda * invGrad[i];
			}
		}


		//std::cout << worstarea * N * 100 << "%" << std::endl;
		return worstarea * N * 100;
	}


	double optimize(int max_Newton_iter = 100, double worst_area_relative_threshold_percent = 0.5, double* resulting_worst_area_error = NULL) {
		V.N = N;
		V.weights.resize(N);
		std::fill(V.weights.begin(), V.weights.end(), 0.0); // need to start with voronoi: it guarantees no cells are empty

		double worst_error = newton_optimize(max_Newton_iter, worst_area_relative_threshold_percent);

		if (resulting_worst_area_error)
			*resulting_worst_area_error = worst_error;

		V.clear();
		V.build();
		double sI = 0;
		for (size_t i = 0; i < N; i++) {
			sI += V.cells[i].integrateSquaredDistance2D(V.samples[i]);
		}
		return sI;
	}

	VoronoiDiagram2D V;
	size_t N;
	std::vector<double> masses;
};
