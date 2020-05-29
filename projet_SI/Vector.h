#include <math.h>
#include <vector>
#include <random>
#include <iostream>

static std::default_random_engine generator;
static std::uniform_real_distribution<double> uniform(0, 1);

class Vector;

Vector operator+(const Vector& a, const Vector &b);

Vector operator-(const Vector& a, const Vector &b);

Vector operator-(const Vector& a);

Vector operator*(double a, const Vector &b);

Vector operator*(const Vector &b, double a);

Vector operator*(const Vector &a, const Vector &b);

Vector operator/(const Vector &a, double b);

Vector operator/(const Vector& a, const Vector &b);

double dot(const Vector& a, const Vector& b);

Vector cross(const Vector& a, const Vector& b);

Vector random_cos(const Vector& N);

Vector random_phong(const Vector& R, double phong_exposant);


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

	Vector reflect(const Vector& N) const{
		Vector result = *this - 2 * dot(*this, N) * N;
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

class Object {
	public:
		Object(const Vector& couleur, bool is_mirror, bool is_transparent, double ks, double phong_exposant): albedo(couleur), is_mirror(is_mirror), is_transparent(is_transparent), phong_exposant(phong_exposant), ks(ks) {};
		virtual bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const = 0;
		Vector albedo;
		bool is_mirror;
		bool is_transparent;
		double phong_exposant;
		double ks; // à changer pour plus brillant, facteur mutliplicatif de l'albedo pour que la spec soit de la meme couleur que la couleur
};


class Triangle: public Object {
	public:
		Triangle(const Vector& A, const Vector &B, const Vector& C, const Vector& couleur, bool is_mirror = false, bool is_transparent = false, double ks = 0, double phong_exposant = 1000): Object(couleur, is_mirror, is_transparent, ks, phong_exposant), A(A), B(B), C(C) {};

		bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const {
			N = cross(B -A, C-A).getNormalized();
			t = dot(C-d.origin, N) / dot(d.direction, N);
			if(t < 0) return false;

			P = d.origin + t*d.direction;
			Vector u = B-A;
			Vector v = C-A;
			Vector w = P-A;
			double m11 = u.getNorm2();
			double m12 = dot(u, v);
			double m22 = v.getNorm2();
			double detm = m11*m22 - m12*m12;

			double b11 = dot(w, u);
			double b21 = dot(w, v);
			double detb = b11*m22 - b21*m12;
			double beta = detb / detm; // coordonnée barycentrique w.r.t à B

			double g12 = b11;
			double g22 = b21;
			double detg = m11*g22 - m12*g12;
			double gamma = detg / detm; // coordonnée barycentrique w.r.t à C

			double alpha = 1 - beta - gamma;
			if(alpha<0 || alpha > 1) return false;
			if(beta<0|| beta > 1) return false;
			if(gamma<0|| gamma > 1) return false;
			if(alpha+beta+gamma > 1) return false;

			return true;
		}
		Vector A, B, C;
};

class Rectangle: public Object {
	public:
		Rectangle(const Vector& A, const Vector &B, const Vector& C, const Vector& D, const Vector& couleur, bool is_mirror = false, bool is_transparent = false, double ks = 0, double phong_exposant = 1000): Object(couleur, is_mirror, is_transparent, ks, phong_exposant), A(A), B(B), C(C), D(D) {};

		bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const {
			N = cross(B -A, D-A).getNormalized();
			t = dot(D-d.origin, N) / dot(d.direction, N);
			if(t < 0) return false;

			P = d.origin + t*d.direction;
			Vector v1 = (B - A).getNormalized();
			Vector v2 = (C - B).getNormalized();
			Vector v3 = (D - C).getNormalized();
			Vector v4 = (A - D).getNormalized();

			Vector v5 = (P - A).getNormalized();
			Vector v6 = (P - B).getNormalized();
			Vector v7 = (P - C).getNormalized();
			Vector v8 = (P - D).getNormalized();
			if(dot(v1, v5) < 0 || dot(v1, v5) > 1) return false;
			if(dot(v2, v6) < 0 || dot(v2, v6) > 1) return false;
			if(dot(v3, v7) < 0 || dot(v3, v7) > 1) return false;
			if(dot(v4, v8) < 0 || dot(v4, v8) > 1) return false;

			return true;
		}

		Vector A, B, C, D;
};

class Cylindre: public Object {
	public: 
		Cylindre(const Vector &origin, double rayon, double hauteur, const Vector& couleur, bool is_mirror = false, bool is_transparent = false, double ks = 0, double phong_exposant = 1000): Object(couleur, is_mirror, is_transparent, ks, phong_exposant), O(origin), R(rayon), H(hauteur){};

		bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const {
			
			double a = (d.direction[0] * d.direction[0]) + (d.direction[2] * d.direction[2]);
			double b = 2 * (d.direction[0]*(d.origin[0]-O[0]) + d.direction[2]*(d.origin[2]-O[2]));
			double c = (d.origin[0] - O[0]) * (d.origin[0] - O[0]) + (d.origin[2] - O[2]) * (d.origin[2] - O[2]) - (R*R);
			
			double delta = b*b - 4*a*c;
			if(delta < 0) return false;
			
			
			double t1 = (-b - std::sqrt(delta))/(2*a);
			double t2 = (-b + std::sqrt(delta))/(2*a);
			
			if (t1 <= 0) t = t2;
			else t = t1;
			if (t < 0) return false;
			
			double r = d.origin[1] + t*d.direction[1];
			P = d.origin + t * d.direction;
			N = Vector (P[0]-O[0],P[1]-O[1],P[2]-O[2]).getNormalized();
			
			if ((r >= O[1]) && (r <= O[1] + H))return true;
			else return false;
		}

		Vector O;
		double R;
		double H;
};

class Sphere : public Object {
public:
	Sphere(const Vector &origin, double rayon, const Vector& couleur, bool is_mirror, bool is_transparent, double ks = 0, double phong_exposant = 1000) : O(origin), R(rayon), Object(couleur, is_mirror, is_transparent, ks, phong_exposant) {};
	
	bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const{
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
};

class Disk : public Object{
public:
	Disk(const Vector& center, const Vector& normal, const Vector& couleur, double radius, bool is_mirror = false, bool is_transparent = false, double ks = 0, double phong_exposant = 1000): c(center),n(normal), r(radius), Object(couleur, is_mirror, is_transparent, ks, phong_exposant) {
	}
	bool intersection(const Ray& d, Vector& P, Vector &N, double &t) const {
		N = n;
		t = dot(c - d.origin, N) / dot(d.direction, N);
		if (t < 0 || t != t) return false;  //isnan

		P = d.origin + t*d.direction;
		double r2 = (P - c).getNorm2();
		return (r2 <= r*r);
	}
	const Vector& c;
	const Vector& n;
	double r;
};

class Scene {
public:
	Scene() {};

	void addSphere(const Sphere* s) {
		objects.push_back(s);
	}
	void addTriangle(const Triangle* s) {
		objects.push_back(s);
	}

	void addCylindre(const Cylindre* s) {
		objects.push_back(s);
	}
	void addRectangle(const Rectangle* s) {
		objects.push_back(s);
	}

	void addDisk(const Disk* s) {
		objects.push_back(s);
	}

	bool intersection(const Ray& d, Vector& P, Vector& N, int& id, double& min) const{

		min = 1E99;
		bool inter = false;

		for(int i = 0; i < objects.size(); i++){

			Vector lP;
			Vector lN;
			double t;
			bool got_inter = objects[i] -> intersection(d, lP, lN, t);
			
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

	// Tous les objets qui composent une scene
	//std::vector<Sphere> spheres;
	std::vector<const Object*> objects;
	Sphere* lumiere;
	double intensite_lumiere;
};