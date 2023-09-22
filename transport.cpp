#include "transport.h"

namespace transport {

	Vector operator*(double a, const Vector& B) {
		return Vector(a * B[0], a * B[1]);
	}
	Vector operator*(const Vector& B, double a) {
		return Vector(a * B[0], a * B[1]);
	}
	Vector operator*(const Vector& A, const Vector& B) {
		return Vector(A[0] * B[0], A[1] * B[1]);
	}
	Vector operator/(const Vector& B, double a) {
		return Vector(B[0] / a, B[1] / a);
	}
	Vector operator/(const Vector& B, const Vector& A) {
		return Vector(B[0] / A[0], B[1] / A[1]);
	}
	Vector operator+(const Vector& A, const Vector& B) {
		return Vector(A[0] + B[0], A[1] + B[1]);
	}
	Vector operator-(const Vector& A, const Vector& B) {
		return Vector(A[0] - B[0], A[1] - B[1]);
	}
	Vector operator-(const Vector& A) {
		return Vector(-A[0], -A[1]);
	}
	double dot(const Vector& A, const Vector& B) {
		return A[0] * B[0] + A[1] * B[1];
	}




	// not used anymore since using Shewchuk's way to checking inCircle. Could be useful anyway.
	Vector circumcenter(const Vector& A, const Vector& B, const Vector& C) {
		Vector u = B - A, v = C - A, M = (A + B) * 0.5, N = (A + C) * 0.5;
		Vector up(u[1], -u[0]), vp(v[1], -v[0]);
		double det = u[0] * v[1] - u[1] * v[0];
		return (dot(u, M) * vp - dot(v, N) * up) / det;
	}


	double det22(const Vector& A, const Vector& B) {
		return A[0] * B[1] - A[1] * B[0];
	}


	// https://www.cs.cmu.edu/~quake/robust.html
	bool isInCircle(const Vector& P, const Vector& A, const Vector& B, const Vector& C) {
		Vector PA = A - P;
		Vector PB = B - P;
		Vector PC = C - P;
		double det1 = det22(PB, PC);
		double det2 = det22(PA, PC);
		double det3 = det22(PA, PB);
		return PA.getNorm2() * det1 - PB.getNorm2() * det2 + PC.getNorm2() * det3 > 0;

	}


	// https://www.cs.cmu.edu/~quake/robust.html but with weights and without the "with a little more effort" transform. height is the parabolic lifting with weight: x^2+y^2+w^2
	bool isInCircleLifted(const Vector& P, const Vector& A, const Vector& B, const Vector& C, double hP, double hA, double hB, double hC) {
		Vector PA = A - P;
		Vector PB = B - P;
		Vector PC = C - P;
		double det1 = det22(PB, PC);
		double det2 = det22(PA, PC);
		double det3 = det22(PA, PB);
		return (hA - hP) * det1 - (hB - hP) * det2 + (hC - hP) * det3 > 0;
	}


	// Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates, JR Shewchuk DCG 1997
	// Use expansion arithmetic to produce robust orientation predicates (not used for inCircle as it seems already much more robust)
	// basic idea: if some operation is within precision limits, convert it to an approximate value and a remainder such that the sum of both is still the exact value.


	// "normal" orient 2d with 2x2 determinant. If within precision limits, just return 0
	inline int orient_2d_filter(const Vector& p0, const Vector& p1, const Vector p2) {

		double a = (p1[0] - p0[0]) * (p2[1] - p0[1]);
		double b = (p1[1] - p0[1]) * (p2[0] - p0[0]);
		double delta = a - b;

		double eps = (4 * std::numeric_limits<double>::epsilon() * (abs(a) + abs(b)));
		if (delta > eps)
			return 1;
		if (delta < -eps)
			return -1;
		return 0;
	}

	// arithmetic presented in Shewchuk's paper
	void fast_two_sum(double a, double b, double& x, double& y) {
		x = a + b;
		double bvirtual = x - a;
		y = b - bvirtual;
	}
	void fast_two_diff(double a, double b, double& x, double& y) {
		x = a - b;
		double bvirtual = a - x;
		y = bvirtual - b;
	}
	void two_sum(double a, double b, double& x, double& y) {
		x = a + b;
		double bvirtual = x - a;
		double avirtual = x - bvirtual;
		double broundoff = b - bvirtual;
		double aroundoff = a - avirtual;
		y = aroundoff + broundoff;
	}
	void two_diff(double a, double b, double& x, double& y) {
		x = a - b;
		double bvirtual = a - x;
		double avirtual = x + bvirtual;
		double broundoff = bvirtual - b;
		double aroundoff = a - avirtual;
		y = aroundoff + broundoff;
	}

	const double dekkerSplitter = (1 << 27) + 1;

	void split(double a, double& ahi, double& alo) {
		double abig = dekkerSplitter - a;
		ahi = dekkerSplitter - abig;
		alo = a - ahi;
	}
	void two_product(double a, double b, double& x, double& y) {
		x = a * b;
		double ahi, alo, bhi, blo;
		split(a, ahi, alo);
		split(b, bhi, blo);
		double err1 = x - (ahi * bhi);
		double err2 = err1 - (alo * bhi);
		double err3 = err2 - (ahi * blo);
		y = (alo * blo) - err3;
	}


	// Fig 21 of Shewchuk's paper
	// Computes the determinant exactly
	int orient_exact(double ax, double ay, double bx, double by, double cx, double cy) {

		double x1, y1, x2, y2, x3, y3, x4, y4;
		two_diff(ax, cx, x1, y1);
		two_diff(by, cy, x2, y2);
		two_diff(ay, cy, x3, y3);
		two_diff(bx, cx, x4, y4);

		double x5a, y5a, x5b, y5b, x5c, y5c, x5d, y5d;
		two_product(x1, x2, x5a, y5a);
		two_product(y1, x2, x5b, y5b);
		two_product(x1, y2, x5c, y5c);
		two_product(y1, y2, x5d, y5d);

		double x6a, y6a, x6b, y6b, x6c, y6c, x6d, y6d;
		two_product(x3, x4, x6a, y6a);
		two_product(y3, x4, x6b, y6b);
		two_product(x3, y4, x6c, y6c);
		two_product(y3, y4, x6d, y6d);

		/*double A = x5a - x6a;
		double A2 = x5b - x6b;
		double A3 = x5c - x6c;*/

		Expansion<4> B1, B2, B3, B4;
		B1.expansion_diff(x5a, y5a, x6a, y6a);
		B2.expansion_diff(x5b, y5b, x6b, y6b);
		B3.expansion_diff(x5c, y5c, x6c, y6c);
		B4.expansion_diff(x5d, y5d, x6d, y6d);

		Expansion<8> S1, S2;
		S1.fromExp(B1);
		S1.sum(B2);
		S2.fromExp(B3);
		S2.sum(B4);

		Expansion<16> D;
		D.fromExp(S1);
		D.sum(S2);

		return D.sign();
	}


	// the actual 2D orientation code : first try the approximate way, and if it fails, the exact (but much slower) predicate
	int orient(const Vector& A, const Vector& B, const Vector C) {

		int result = orient_2d_filter(A, B, C);

		if (result == 0) {

			// also tried below a perturbation theory formula, but this results in additional degenerate triangles on the convex hull in some cases
			/*const Vector U(B - A), V(C - A);
			if (V[1] > 0) return 1;
			if (V[1] < 0) return -1;
			double x = U[0] - U[1] - V[0];
			if (x > 0) return 1;
			if (x < 0) return -1;
			return 1;*/

			return orient_exact(A[0], A[1], B[0], B[1], C[0], C[1]);

		}
		return result;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////// This part of the code below contains helpers for ordering vertices along a Morton curve                                             //////////////////////// 
	////////  It has been generated with ChatGPT 4 (no other part of the code is AI, apart from the actual comparison code Bowyer2D::cmp_zorder) //////////////////////// 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Convert a normalized double (in range [0, 1)) to its 64-bit integer representation
	int64_t doubleToRawInt(double value) {
		return static_cast<int64_t>(value * (static_cast<double>(INT64_MAX) - 1.0));
	}

	bool less_msb(int64_t x, int64_t y) {
		return x < y&& x < (x^ y);
	}

	//////////////////////////// end of ChatGPT Code /////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////

};