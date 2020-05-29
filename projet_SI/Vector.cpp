#include "Vector.h"

Vector operator+(const Vector& x, const Vector &y) {
	return Vector(x[0] + y[0], x[1] + y[1], x[2] + y[2]);
}

Vector operator-(const Vector& x, const Vector& y) {
	return Vector(x[0] - y[0], x[1] - y[1], x[2] - y[2]);
}

Vector operator-(const Vector& x) {
	return Vector(-x[0], -x[1], -x[2]);
}

Vector operator*(double x, const Vector &y) {
	return Vector(x * y[0], x * y[1], x * y[2]);
}

Vector operator*(const Vector &y, double x) {
	return Vector(x * y[0], x * y[1], x * y[2]);
}

Vector operator*(const Vector &x, const Vector &y) {
	return Vector(x[0] * y[0], x[1] * y[1], x[2] * y[2]);
}

Vector operator/(const Vector &x, double y) {
	return Vector(x[0] / y, x[1] / y, x[2] / y);
}

Vector operator/(const Vector& x, const Vector &y) {
	return Vector(x[0] / y[0], x[1] / y[1], x[2] / y[2]);
}

double dot(const Vector&x, const Vector& y) {
	return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

Vector cross(const Vector& x, const Vector& y) {
	return Vector(x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]);
}
	
Vector random_cos(const Vector& N){
	double r1 = generator_rand(generator);
	double r2 = generator_rand(generator);
	Vector rand_loc_dir(cos(2 * M_PI * r1) * sqrt(1 - r2), sin(2 * M_PI * r1) * sqrt(1 - r2), sqrt(r2));
	Vector random_vec(generator_rand(generator), generator_rand(generator), generator_rand(generator));
	Vector tang1 = cross(N, random_vec);
	tang1.normalize();
	Vector tang2 = cross(tang1, N);
	return rand_loc_dir[2] * N + rand_loc_dir[0] * tang1 + rand_loc_dir[1] * tang2;	
}

Vector random_phong(const Vector& R, double phong_exposant){
	double r1 = generator_rand(generator);
	double r2 = generator_rand(generator);
	double fac = sqrt(1 - std::pow(r2, 2. / (phong_exposant + 1)));
	Vector rand_loc_dir(cos(2 * M_PI * r1) * fac, sin(2 * M_PI * r1) * fac, std::pow(r2, 1./(phong_exposant + 1)));
	Vector random_vec(generator_rand(generator) - 0.5, generator_rand(generator) - 0.5, generator_rand(generator) - 0.5);
	Vector tang1 = cross(R, random_vec);
	tang1.normalize();
	Vector tang2 = cross(tang1, R);
	return rand_loc_dir[2] * R + rand_loc_dir[0] * tang1 + rand_loc_dir[1] * tang2;	

}