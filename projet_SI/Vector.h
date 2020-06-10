/*!
 * \file Vector.h
 * \brief Define Vectors and Rays
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
* \brief Allows to do apply a random phong on a vector
* 
* \param R : the R vector
* \param phong_expo : the phong exponent 
* \return the vector which is the result of the operation
*/
Vector random_phong(const Vector& R, double phong_expo);

/*!
* \brief Allows to do apply a random cos on a vector
* 
* \param N : the N vector
* \return the vector which is the result of the operation
*/
Vector random_cos(const Vector& N);

/*!
* \brief Allows to do vectorial product on vectors
* Method that define vectorial product for vectors
* 
* \param x : the x vector
* \param y : the y vector
* \return the vector which is the result of the operation
*/
Vector cross(const Vector& x, const Vector& y);

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
	double coord[3];/*!< coordinates of the vector */
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
    Vector origin, direction;/*!< origin and direction of the ray */
};

