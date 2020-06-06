#include <random>
#include <iostream>
#include <math.h>
#include <vector>

static std::uniform_real_distribution<double> generator_rand(0, 1);
static std::default_random_engine generator;

class Vector;

Vector operator+(const Vector& x, const Vector &y);
Vector operator-(const Vector& x, const Vector &y);
Vector operator-(const Vector& x);
Vector operator*(double x, const Vector &y);
Vector operator*(const Vector &y, double x);
Vector operator*(const Vector &x, const Vector &y);
Vector operator/(const Vector &x, double y);
Vector operator/(const Vector& x, const Vector &y);
double dot(const Vector& x, const Vector& y);
Vector random_phong(const Vector& R, double phong_expo);
Vector random_cos(const Vector& N);
Vector cross(const Vector& x, const Vector& y);


class Vector {
public:

	Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	const double& operator[](int i) const { 
		return coord[i]; 
	}

	Vector reflect(const Vector& N) const{
		Vector result = *this - 2 * dot(*this, N) * N;
		return result;
	}

	/* renvoir norme vecteur */
	double getNorm() {
		return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
	}

	void normalize() {
		double norm = sqrt(getNorm());
		coord[0] /= norm;
		coord[1] /= norm;
		coord[2] /= norm;
		
	}

	/* renvoi vecteur normalisé */
	Vector getNormalizedCopy() {
		Vector result(*this);
		result.normalize();
		return result;
	}



private:
	double coord[3];
};

class Ray {
public:
	Ray() {};
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
    Vector origin, direction;
};

