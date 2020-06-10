/*!
 * \file Shape.h
 * \brief Define all objects that can be placed in a scene
 * \author Pierre-Jean Besnard & Louis Billaut
 * \version 1.0
 */
#include "Vector.h"

/*! \class Shape
 * \brief Allows to create shape and define a color, or if this shape will be a mirror, transparent, his brightness and his phong exponant
 */
class Shape {
	public:

		/*!
		 * \brief Constructor
		 * Constructor of the Shape class
		 * 
		 * \param color : the color of the shape, representing by a vector, RGB
		 * \param is_mirror : true if the shape is a mirror, false either
		 * \param is_transparent : true if the shape is transparent, false either
		 * \param spec : the brightness of the shape, representing by a double
		 * \param phong_expo : the phong exponent of the shape, representing by a double
		 */
		Shape(const Vector& color, bool is_mirror, bool is_transparent, double spec, double phong_expo): color(color), is_mirror(is_mirror), is_transparent(is_transparent), phong_expo(phong_expo), spec(spec) {};
		
		/*!
		 * \brief Calculate the intersection of a ray with the shape if it exist
		 * Calculate the P intersection point, N the normal intersection and t the intensity 
		 * 
		 * \param d : the ray within the triangle intersect
		 * \param P : the intersection point wich will be calculated
		 * \param N : the normal wich will be calculated
		 * \param t : the intensity wich will be calculated
		 * 
		 * \return true if the light intersect the shape, false either
		 */ 
		virtual bool intersection(const Ray& d, Vector& P, Vector& N, double& t) const = 0;
		Vector color; /*!< color of the vector, RGB */ 
		bool is_mirror;/*!< define if the shape is a mirror or not */ 
		bool is_transparent;/*!< define if the shape is transparent or not */ 
		double phong_expo;/*!< the phong exponent of the shape */ 
		double spec; /*!< the brightness of the shape */ 
};

/*! \class Triangle
 * \brief Allows to create Triangle shape and to calculate intersections with him
 */
class Triangle: public Shape {
	public:
		/*!
		 * \brief Constructor
		 * Constructor of the triangle class
		 * 
		 * \param A : the A point of the triangle
		 * \param B : the B point of the triangle
		 * \param C : the C point of the triangle
		 * \param color : the color of the triangle, representing by a vector, RGB
		 * \param is_mirror : true if the triangle is a mirror, false either
		 * \param is_transparent : true if the triangle is transparent, false either
		 * \param spec : the brightness of the triangle, representing by a double
		 * \param phong_expo : the phong exponent of the triangle, representing by a double
		 */
		Triangle(const Vector& A, const Vector &B, const Vector& C, const Vector& color, bool is_mirror = false, bool is_transparent = false, double spec = 0, double phong_expo = 1000): Shape(color, is_mirror, is_transparent, spec, phong_expo), A(A), B(B), C(C) {};

		/*!
		 * \brief Calculate the intersection of a ray with the triangle if it exist
		 * Calculate the P intersection point, N the normal intersection and t the intensity with the d ray
		 * 
		 * \param d : the ray within the triangle intersect
		 * \param P : the intersection point wich will be calculated
		 * \param N : the normal wich will be calculated
		 * \param t : the intensity wich will be calculated
		 * 
		 * \return true if the light intersect the triangle, false either
		 */ 
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
		Vector A, B, C;/*!< the A, B, and C points of the triangle */ 
};

/*! \class Rectangle
 * \brief Allows to create Rectangle shape and to calculate intersections with him
 */
class Rectangle: public Shape {
	public:

		/*!
		 * \brief Constructor
		 * Constructor of the rectangle class
		 * 
		 * \param A : the A point of the rectangle
		 * \param B : the B point of the rectangle
		 * \param C : the C point of the rectangle
		 * \param D : the D point of the rectangle
		 * \param color : the color of the rectangle, representing by a vector, RGB
		 * \param is_mirror : true if the rectangle is a mirror, false either
		 * \param is_transparent : true if the rectangle is transparent, false either
		 * \param spec : the brightness of the rectangle, representing by a double
		 * \param phong_expo : the phong exponent of the rectangle, representing by a double
		 */
		Rectangle(const Vector& A, const Vector &B, const Vector& C, const Vector& D, const Vector& color, bool is_mirror = false, bool is_transparent = false, double spec = 0, double phong_expo = 1000): Shape(color, is_mirror, is_transparent, spec, phong_expo), A(A), B(B), C(C), D(D) {};

		/*!
		 * \brief Calculate the intersection of a ray with the rectangle if it exist
		 * Calculate the P intersection point, N the normal intersection and t the intensity with the d ray
		 * 
		 * \param d : the ray within the rectangle intersect
		 * \param P : the intersection point wich will be calculated
		 * \param N : the normal wich will be calculated
		 * \param t : the intensity wich will be calculated
		 * 
		 * \return true if the light intersect the rectangle, false either
		 */ 
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

		Vector A, B, C, D;/*!< the A, B, C and D points of the rectangle */ 
};

/*! \class Cylinder
 * \brief Allows to create Cylinder shape and to calculate intersections with him
 */
class Cylinder: public Shape {
	public: 

		/*!
		 * \brief Constructor
		 * Constructor of the cylinder class
		 * 
		 * \param origin : the origin point of the cylinder
		 * \param rayon : the ray of the cylinder
		 * \param hauteur : the C height of the cylinder
		 * \param color : the color of the cylinder, representing by a vector, RGB
		 * \param is_mirror : true if the cylinder is a mirror, false either
		 * \param is_transparent : true if the cylinder is transparent, false either
		 * \param spec : the brightness of the cylinder, representing by a double
		 * \param phong_expo : the phong exponent of the cylinder, representing by a double
		 */
		Cylinder(const Vector &origin, double rayon, double hauteur, const Vector& color, bool is_mirror = false, bool is_transparent = false, double spec = 0, double phong_expo = 1000): Shape(color, is_mirror, is_transparent, spec, phong_expo), O(origin), R(rayon), H(hauteur){};

		/*!
		 * \brief Calculate the intersection of a ray with the cylinder if it exist
		 * Calculate the P intersection point, N the normal intersection and t the intensity with the d ray
		 * 
		 * \param d : the ray within the cylinder intersect
		 * \param P : the intersection point wich will be calculated
		 * \param N : the normal wich will be calculated
		 * \param t : the intensity wich will be calculated
		 * 
		 * \return true if the light intersect the cylinder, false either
		 */ 
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

		Vector O;/*!< the origin points of the cylinder*/ 
		double R;/*!< the ray of the cylinder*/
		double H;/*!< the height of the cylinder*/
};

/*! \class Sphere
 * \brief Allows to create Sphere shape and to calculate intersections with him
 */
class Sphere : public Shape {
public:
	/*!
	* \brief Constructor
	* Constructor of the sphere class
	* 
	* \param origin : the origin point of the sphere
	* \param ray : the ray of the sphere
	* \param color : the color of the sphere, representing by a vector, RGB
	* \param is_mirror : true if the sphere is a mirror, false either
	* \param is_transparent : true if the sphere is transparent, false either
	* \param spec : the brightness of the sphere, representing by a double
	* \param phong_expo : the phong exponent of the sphere, representing by a double
	*/
	Sphere(const Vector &origin, double ray, const Vector& color, bool is_mirror, bool is_transparent, double spec = 0, double phong_expo = 1000) : O(origin), R(ray), Shape(color, is_mirror, is_transparent, spec, phong_expo) {};
	
	/*!
	* \brief Calculate the intersection of a ray with the sphere if it exist
	* Calculate the P intersection point, N the normal intersection and t the intensity with the d ray
	* 
	* \param d : the ray within the sphere intersect
	* \param P : the intersection point wich will be calculated
	* \param N : the normal wich will be calculated
	* \param t : the intensity wich will be calculated
	* 
	* \return true if the light intersect the sphere, false either
	*/ 
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
	Vector O; /*!< the origin points of the sphere*/ 
	double R;/*!< the ray of the sphere*/ 
};

