#include <math.h>
#include <vector>

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

	/*double& operator[](int i) const { 
		return coord[i]; 
	}*/

	/* renvoir norme vecteur */
	double getNorm2() {
		return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
	}

	void normalize() {
		double norm = sqrt(getNorm2());
		coord[0] /= norm;
		coord[1] /= norm;
		coord[2] /= norm;
		
	}

	/* renvoi vecteur normalisé */
	Vector getNormalized() {
		Vector result(*this);
		result.normalize();
		return result;
	}

private:
	double coord[3];
};

Vector operator+(const Vector& a, const Vector &b);

Vector operator-(const Vector& a, const Vector &b);

Vector operator*(double a, const Vector &b);

Vector operator*(const Vector &b, double a);

Vector operator/(const Vector &a, double b);

Vector operator/(const Vector& a, const Vector &b);

double dot(const Vector& a, const Vector& b);

class Ray {
public:
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
    Vector origin, direction;
};

class Sphere {
public:
	Sphere(const Vector &origin, double rayon, const Vector& couleur) : O(origin), R(rayon), albedo(couleur){};
	
	bool intersection(const Ray& d, Vector& P, Vector& N, double& t) {
		double a = 1;
		double b = 2 * dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm2() - R * R;

		double delta = b * b - 4 * a * c;
		if (delta < 0) {
			// Pas d'interesection
			return false;
		}
		double t1 = (-b - sqrt(delta)) / (2 * a);
		double t2 = (-b + sqrt(delta)) / (2 * a);

		if (t2 < 0) {
			// intersection derriere camera
			return false;
		}
		
		if (t1 > 0) {
			t = t1;
		}
		else {
			t = t2;
		}

		P = d.origin + t * d.direction; // point d'intersection entre notre rayon lancé et la sphere

		N = (P - O).getNormalized(); 

		return true;
	}

	Vector O;
	double R;
	Vector albedo;
};

class Scene {
public:
	Scene() {};
	void addSphere(const Sphere& s) {
		spheres.push_back(s);
	}

	bool intersection(const Ray& d, Vector& P, Vector& N, int& id) {

		double min = 1E99;
		bool inter = false;

		for(int i = 0; i < spheres.size(); i++){
			Vector lP;
			Vector lN;
			double t;
			bool got_inter = spheres[i].intersection(d, lP, lN, t);
			if (got_inter){
				inter = true;
				if (min > t){
					min = t;
					P = lP;
					N = lN;
					id = i;
				}
			}
		}
		return inter;
	}

	std::vector<Sphere> spheres;
};