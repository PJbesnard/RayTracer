/*!
 * \file Scene.h
 * \brief Define the scene
 * \author Pierre-Jean Besnard & Louis Billaut
 * \version 1.0
 */
#include "Shape.h"

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
	Sphere* light;/*!< light of the scene represented by a sphere*/ 
	double light_intensity; /*!< light intensity of the scene*/ 
};