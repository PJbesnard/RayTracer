#include "Vector.h"

Vector operator+(const Vector& a, const Vector &b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(double a, const Vector &b) {
	return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector &b, double a) {
	return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator/(const Vector &a, double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}

Vector operator/(const Vector& a, const Vector &b) {
	return Vector(a[0] / b[0], a[1] / b[1], a[2] / b[2]);
}

double dot(const Vector&a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
