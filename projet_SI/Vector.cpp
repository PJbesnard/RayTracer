#include "Vector.h"

Vector operator+(const Vector& a, const Vector &b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}

Vector operator*(double a, const Vector &b) {
	return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector &b, double a) {
	return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector &a, const Vector &b) {
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
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

Vector cross(const Vector&a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
	
Vector random_cos(const Vector& N){
	double r1 = uniform(generator);
	double r2 = uniform(generator);
	Vector direction_random_local(cos(2 * M_PI * r1) * sqrt(1 - r2), sin(2 * M_PI * r1) * sqrt(1 - r2), sqrt(r2));
	Vector random_vec(uniform(generator), uniform(generator), uniform(generator));
	Vector tang1 = cross(N, random_vec);
	tang1.normalize();
	Vector tang2 = cross(tang1, N);
	return direction_random_local[2] * N + direction_random_local[0] * tang1 + direction_random_local[1] * tang2;	
}

Vector random_phong(const Vector& R, double phong_exposant){
	double r1 = uniform(generator);
	double r2 = uniform(generator);
	double facteur = sqrt(1 - std::pow(r2, 2. / (phong_exposant + 1)));
	Vector direction_random_local(cos(2 * M_PI * r1) * facteur, sin(2 * M_PI * r1) * facteur, std::pow(r2, 1./(phong_exposant + 1)));
	Vector random_vec(uniform(generator) - 0.5, uniform(generator) - 0.5, uniform(generator) - 0.5);
	Vector tang1 = cross(R, random_vec);
	tang1.normalize();
	Vector tang2 = cross(tang1, R);
	return direction_random_local[2] * R + direction_random_local[0] * tang1 + direction_random_local[1] * tang2;	

}
