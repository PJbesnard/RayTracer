/*!
 * \file Vector.h
 * \brief Define all objects that can be placed in a scene
 * \author Pierre-Jean Besnard & Louis Billaut
 * \version 1.0
 */
#include <random>
#include <iostream>
#include <math.h>
#include <vector>

static std::uniform_real_distribution<double> generator_rand(0, 1); /*!< a random distribution on the range [min, max) with equal probability*/
static std::default_random_engine generator;  /*!< a random number generator*/

/*! \class Vector
 * \brief Allows to create vectors and apply operations to them
 */
class Vector;

/*!
* \brief Allows to do additions on vectors
* Method that define + for vectors
* 
* \param x : the x vector
* \param y : the y vector
* \return the vector which is the result of the operation
*/
Vector operator+(const Vector& x, const Vector &y);

/*!
* \brief Allows to do substraction on vectors
* Method that define substraction for vectors
* 
* \param x : the x vector
* \param y : the y vector
* \return the vector which is the result of the operation
*/
Vector operator-(const Vector& x, const Vector &y);

/*!
* \brief Allows to get the negative vector
* Method that define -vector
* 
* \param x : the x vector
* \return the vector which is the result of the negative vector
*/
Vector operator-(const Vector& x);

/*!
* \brief Allows to do multiplications on vectors
* Method that define multiplications for vectors
* 
* \param x : the x value
* \param y : the y vector
* \return the vector which is the result of the operation
*/
Vector operator*(double x, const Vector &y);

/*!
* \brief Allows to do multiplications on vectors
* Method that define multiplications for vectors
*
* \param y : the y vector
* \param x : the x value
* \return the vector which is the result of the operation
*/
Vector operator*(const Vector &y, double x);

/*!
* \brief Allows to do multiplications on vectors
* Method that define multiplications for vectors
* 
* \param x : the x vector
* \param y : the y vector
* \return the vector which is the result of the operation
*/
Vector operator*(const Vector &x, const Vector &y);

/*!
* \brief Allows to do divisions on vectors
* Method that define divisions for vectors
* 
* \param x : the x vector
* \param y : the y value
* \return the vector which is the result of the operation
*/
Vector operator/(const Vector &x, double y);

/*!
* \brief Allows to do divisions on vectors
* Method that define divisions for vectors
* 
* \param x : the x vector
* \param y : the y vector
* \return the vector which is the result of the operation
*/
Vector operator/(const Vector& x, const Vector &y);

/*!
* \brief Allows to do scalar product on vectors
* Method that define scalar product for vectors
* 
* \param x : the x vector
* \param y : the y vector
* \return the vector which is the result of the operation
*/
double dot(const Vector& x, const Vector& y);

/*!
* \brief Allows to do vectorial product on vectors
* Method that define vectorial product for vectors
* 
* \param x : the x vector
* \param y : the y vector
* \return the vector which is the result of the operation
*/
Vector cross(const Vector& x, const Vector& y);

/*!
* \brief Allows to do apply a random cos on a vector
* 
* \param N : the N vector
* \return the vector which is the result of the operation
*/
Vector random_cos(const Vector& N);

/*!
* \brief Allows to do apply a random phong on a vector
* 
* \param R : the R vector
* \param phong_expo : the phong exponent 
* \return the vector which is the result of the operation
*/
Vector random_phong(const Vector& R, double phong_expo);

/*! \class Vector
 * \brief Allows to create vectors and apply operations to them
 */
class Vector {
public:
	/*!
	* \brief Constructor
	* Constructor of Vector class
	*
	* \param x : x position of the vector
	* \param y : y position of the vector
	* \param z : z position of the vector
	*/
	Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}

	/*!
	 * \brief Allows to get the x, y or z position of the vector
	 * Method that define [] for vector, [0] to get x, [1] to get y, [2] to get z
	 * 
	 * \param i : index of the desired coordinate
	 * \return the desired coordinate
	 */
	const double& operator[](int i) const { 
		return coord[i]; 
	}

	/*!
	 * \brief Allows to get the reflection of the vector regarding to the normal
	 * 
	 * \param N : the normal vector
	 * \return the reflected vector regarding to the normal
	 */
	Vector reflect(const Vector& N) const{
		Vector result = *this - 2 * dot(*this, N) * N;
		return result;
	}

	/*!
	 * \brief Calculate the normal of this vector
	 * 
	 * \return normal of this vector
	 */
	double getNorm() {
		return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
	}

	/*!
	 * \brief Normalize this vector
	 */
	void normalize() {
		double norm = sqrt(getNorm());
		coord[0] /= norm;
		coord[1] /= norm;
		coord[2] /= norm;
		
	}

	/*!
	 * \brief Get this normalized vector 
	 * 
	 * \return this normalized vector 
	 */
	Vector getNormalizedCopy() {
		Vector result(*this);
		result.normalize();
		return result;
	}

private:
	double coord[3]; /*!< coordinates of the vector */
};

/*! \class Ray
 * \brief defines a ray
 */
class Ray {
public:
	/*!
	 * \brief Constructor
	 * 
	 * Constructor of the Ray class
	 */
	Ray() {};
	/*!
	 * \brief Constructor
	 * 
	 * Constructor of the Ray class
	 * 
	 * \param o : origin of the ray, representing by a vector
	 * \param d : direction of the ray, representing by a vector
	 */
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
    Vector origin, direction; /*!< origin and direction of the ray */
};

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
		Shape(const Vector& color, bool is_mirror, bool is_transparent, double spec, double phong_expo): albedo(color), is_mirror(is_mirror), is_transparent(is_transparent), phong_expo(phong_expo), spec(spec) {};
		
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
		Vector albedo; /*!< color of the vector, RGB */ 
		bool is_mirror; /*!< define if the shape is a mirror or not */ 
		bool is_transparent; /*!< define if the shape is transparent or not */ 
		double phong_expo; /*!< the phong exponent of the shape */ 
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
			double beta = detb / detm; // barycentric coordinate to B

			double g12 = b11;
			double g22 = b21;
			double detg = m11 * g22 - m12 * g12;
			double gamma = detg / detm; // barycentric coordinate to C

			double alpha = 1 - beta - gamma;
			if(alpha < 0 || alpha > 1) return false;
			if(beta < 0 || beta > 1) return false;
			if(gamma < 0 || gamma > 1) return false;
			if(alpha + beta + gamma > 1) return false;

			return true;
		}
		Vector A, B, C; /*!< the A, B, and C points of the triangle */ 
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

		Vector A, B, C, D; /*!< the A, B, C and D points of the rectangle */ 
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

		Vector O; /*!< the origin points of the cylinder*/ 
		double R; /*!< the ray of the cylinder*/
		double H; /*!< the height of the cylinder*/
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

		// intersection behind camera
		if (t2 < 0) return false;
		if (t1 > 0) t = t1;
		else t = t2;
		P = d.origin + t * d.direction; // intersection point between ray and sphere
		N = (P - O).getNormalizedCopy(); 
		return true;
	}
	Vector O; /*!< the origin points of the sphere*/ 
	double R; /*!< the ray of the sphere*/ 
};

/*! \class Scene
 * \brief Allows to create a scene, add object to it and calculate intersection with them
 */
class Scene {
public:
	/*!
	 * \brief Constructor
	 * Constructor of the scene class
	 */
	Scene() {};

	/*!
	* \brief Calculate the intersection of a ray with objects of the scene
	* Calculate the P intersection point, N the normal intersection and id the id of the intersected object
	* 
	* \param d : the ray within the sphere intersect
	* \param P : the intersection point wich will be calculated
	* \param N : the normal wich will be calculated
	* \param id : the id of the object wich will be intersected
	* \param min : the minimal value of intensity
	* 
	* \return true if the light intersect the sphere, false either
	*/ 
	bool intersection(const Ray& d, Vector& P, Vector& N, int& id, double& min) const{
		bool is_inter = false;
		min = 1E99;
		for(int i = 0; i < shapes.size(); i++){
			double t;
			Vector lP, lN;
			bool got_inter = shapes[i] -> intersection(d, lP, lN, t);
			if (got_inter){
				is_inter = true;
				if (min > t){
					min = t;
					P = lP;
					N = lN;
					id = i;
				}
			}
		}
		return is_inter;
	}

	/*!
	 * \brief Add sphere object to the scene
	 * 
	 * \param s : the sphere wich will be added to the scene
	 */
	void addSphere(const Sphere* s) {
		shapes.push_back(s);
	}

	/*!
	 * \brief Add triangle object to the scene
	 * 
	 * \param s : the triangle wich will be added to the scene
	 */
	void addTriangle(const Triangle* s) {
		shapes.push_back(s);
	}

	/*!
	 * \brief Add cylinder object to the scene
	 * 
	 * \param s : the cylinder wich will be added to the scene
	 */
	void addCylinder(const Cylinder* s) {
		shapes.push_back(s);
	}

	/*!
	 * \brief Add rectangle object to the scene
	 * 
	 * \param s : the rectangle wich will be added to the scene
	 */
	void addRectangle(const Rectangle* s) {
		shapes.push_back(s);
	}

	std::vector<const Shape*> shapes; /*!< objects of the scene*/ 
	Sphere* light; /*!< light of the scene represented by a sphere*/ 
	double light_intensity; /*!< light intensity of the scene*/ 
};