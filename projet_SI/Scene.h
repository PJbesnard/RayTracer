#include "Shape.h"

class Scene {
public:
	Scene() {};

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

	void addSphere(const Sphere* s) {
		shapes.push_back(s);
	}
	void addTriangle(const Triangle* s) {
		shapes.push_back(s);
	}

	void addCylinder(const Cylinder* s) {
		shapes.push_back(s);
	}
	void addRectangle(const Rectangle* s) {
		shapes.push_back(s);
	}



	// Tous les objets qui composent une scene
	std::vector<const Shape*> shapes;
	Sphere* light;
	double light_intensity;
};