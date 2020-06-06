#include "Vector.h"

class Shape {
	public:
		Shape(const Vector& color, bool is_mirror, bool is_transparent, double spec, double phong_expo): color(color), is_mirror(is_mirror), is_transparent(is_transparent), phong_expo(phong_expo), spec(spec) {};
		virtual bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const = 0;
		Vector color;
		bool is_mirror;
		bool is_transparent;
		double phong_expo;
		double spec; // à changer pour plus brillant, facteur mutliplicatif de l'albedo pour que la spec soit de la meme couleur que la couleur
};


class Triangle: public Shape {
	public:
		Triangle(const Vector& A, const Vector &B, const Vector& C, const Vector& color, bool is_mirror = false, bool is_transparent = false, double spec = 0, double phong_expo = 1000): Shape(color, is_mirror, is_transparent, spec, phong_expo), A(A), B(B), C(C) {};

		bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const {
			N = cross(B - A, C - A).getNormalizedCopy();
			t = dot(C - d.origin, N) / dot(d.direction, N);
			if(t < 0) return false;

			P = d.origin + t * d.direction;
			Vector u = B - A;
			Vector v = C - A;
			Vector w = P - A;
			double m11 = u.getNorm();
			double m12 = dot(u, v);
			double m22 = v.getNorm();
			double detm = m11 * m22 - m12 * m12;

			double b11 = dot(w, u);
			double b21 = dot(w, v);
			double detb = b11 * m22 - b21 * m12;
			double beta = detb / detm; // coordonnée barycentrique w.r.t à B

			double g12 = b11;
			double g22 = b21;
			double detg = m11 * g22 - m12 * g12;
			double gamma = detg / detm; // coordonnée barycentrique w.r.t à C

			double alpha = 1 - beta - gamma;
			if(alpha < 0 || alpha > 1) return false;
			if(beta < 0 || beta > 1) return false;
			if(gamma < 0 || gamma > 1) return false;
			if(alpha + beta + gamma > 1) return false;

			return true;
		}
		Vector A, B, C;
};

class Rectangle: public Shape {
	public:
		Rectangle(const Vector& A, const Vector &B, const Vector& C, const Vector& D, const Vector& color, bool is_mirror = false, bool is_transparent = false, double spec = 0, double phong_expo = 1000): Shape(color, is_mirror, is_transparent, spec, phong_expo), A(A), B(B), C(C), D(D) {};

		bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const {
			N = cross(B -A, D-A).getNormalizedCopy();
			t = dot(D-d.origin, N) / dot(d.direction, N);
			if(t < 0) return false;

			P = d.origin + t*d.direction;
			Vector v1 = (B - A).getNormalizedCopy();
			Vector v2 = (C - B).getNormalizedCopy();
			Vector v3 = (D - C).getNormalizedCopy();
			Vector v4 = (A - D).getNormalizedCopy();

			Vector v5 = (P - A).getNormalizedCopy();
			Vector v6 = (P - B).getNormalizedCopy();
			Vector v7 = (P - C).getNormalizedCopy();
			Vector v8 = (P - D).getNormalizedCopy();
			if(dot(v1, v5) < 0 || dot(v1, v5) > 1) return false;
			if(dot(v2, v6) < 0 || dot(v2, v6) > 1) return false;
			if(dot(v3, v7) < 0 || dot(v3, v7) > 1) return false;
			if(dot(v4, v8) < 0 || dot(v4, v8) > 1) return false;

			return true;
		}

		Vector A, B, C, D;
};

class Cylinder: public Shape {
	public: 
		Cylinder(const Vector &origin, double rayon, double hauteur, const Vector& color, bool is_mirror = false, bool is_transparent = false, double spec = 0, double phong_expo = 1000): Shape(color, is_mirror, is_transparent, spec, phong_expo), O(origin), R(rayon), H(hauteur){};

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
			N = Vector (P[0]-O[0],P[1]-O[1],P[2]-O[2]).getNormalizedCopy();
			
			if ((r >= O[1]) && (r <= O[1] + H))return true;
			else return false;
		}

		Vector O;
		double R;
		double H;
};

class Sphere : public Shape {
public:
	Sphere(const Vector &origin, double ray, const Vector& color, bool is_mirror, bool is_transparent, double spec = 0, double phong_expo = 1000) : O(origin), R(ray), Shape(color, is_mirror, is_transparent, spec, phong_expo) {};
	
	bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const{
		double a = 1;
		double b = 2 * dot(d.direction, d.origin - O);
		double c = (d.origin - O).getNorm() - R * R;
		double delta = b * b - 4 * a * c;
		if (delta < 0) return false;
		double t1 = (-b - sqrt(delta)) / (2 * a);
		double t2 = (-b + sqrt(delta)) / (2 * a);

		// intersection derriere camera
		if (t2 < 0) return false;
		if (t1 > 0) t = t1;
		else t = t2;
		P = d.origin + t * d.direction; // point d'intersection entre notre rayon lancé et la sphere
		N = (P - O).getNormalizedCopy(); 
		return true;
	}
	Vector O;
	double R;
};

