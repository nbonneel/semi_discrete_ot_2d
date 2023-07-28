#include "transport.h"


bool operator==(const Vector& a, const Vector& b) {
	if (abs(a[0] - b[0]) > 1E-12) return false;
	if (abs(a[1] - b[1]) > 1E-12) return false;
	if (abs(a[2] - b[2]) > 1E-12) return false;
	return true;
}
Vector operator*(double a, const Vector& B) {
	return Vector(a * B[0], a * B[1], a * B[2]);
}
Vector operator*(const Vector& B, double a) {
	return Vector(a * B[0], a * B[1], a * B[2]);
}
Vector operator*(const Vector& A, const Vector& B) {
	return Vector(A[0] * B[0], A[1] * B[1], A[2] * B[2]);
}
Vector operator/(const Vector& B, double a) {
	return Vector(B[0] / a, B[1] / a, B[2] / a);
}
Vector operator+(const Vector& A, const Vector& B) {
	return Vector(A[0] + B[0], A[1] + B[1], A[2] + B[2]);
}
Vector operator-(const Vector& A, const Vector& B) {
	return Vector(A[0] - B[0], A[1] - B[1], A[2] - B[2]);
}
Vector operator-(const Vector& A) {
	return Vector(-A[0], -A[1], -A[2]);
}
double dot_3d(const Vector& A, const Vector& B) {
	return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}
double dot_2d(const Vector& A, const Vector& B) {
	return A[0] * B[0] + A[1] * B[1];
}
