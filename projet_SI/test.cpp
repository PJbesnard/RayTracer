#include <stdio.h>
#include <iostream>
#include <vector>
#include "Vector.h"



void save_img(const char* filename, const unsigned char* pixels, int W, int H) {

	FILE* f = fopen(filename, "wb");
	fprintf(f, "P3\n");  
	fprintf(f, "%d %d\n", W, H);
	fprintf(f, "255\n");  

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			int r = pixels[(i * W + j) * 3 + 0];
			fprintf(f, "%d ", r);
			int g = pixels[(i * W + j) * 3 + 1];
			fprintf(f, "%d ", g);
			int b = pixels[(i * W + j) * 3 + 2];
			fprintf(f, "%d ", b);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}


bool intersection(const Ray& d, const Sphere& s, Vector& P, Vector& N) {
	double a = 1;
	double b = 2 * dot(d.direction, d.origin - s.O);
	double c = (d.origin - s.O).getNorm2() - s.R * s.R;

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
	
	double t;
	if (t1 > 0) {
		t = t1;
	}
	else {
		t = t2;
	}

	P = d.origin + t * d.direction; // point d'intersection entre notre rayon lancé et la sphere

	N = (P - s.O).getNormalized(); 

	return true;
}

 double carre(double n){
 	return n * n;
 }


// const enmpeche de faire une copie locale
// On fait une fonction recursive pour calculer la reflexion des surfaces
// La fonction renvoie la couleur du pixel obtenu en envoyant un rayon R dans la scene S
Vector getColor(const Ray& r, const Scene& s, int recursion = 0){
	if (recursion >= 5){
		return Vector(0, 0, 0);
	}

	int id;
	double t;
	Vector P, N;
	bool has_inter = s.intersection(r, P, N, id, t);
	Vector intensite_pixel(0,0,0); // trois cmposantes du coup on prend un vecteur (A CHANGER)
			
	if (has_inter) {
		
		if (s.spheres[id].is_mirror){
			Vector direction_mirroir = r.direction - 2 * dot(N, r.direction) * N;
			Ray rayon_mirroir(P + 0.001*N, direction_mirroir);

			intensite_pixel = getColor(rayon_mirroir, s, recursion + 1);
		} 
		else 
			if (s.spheres[id].is_transparent){
			double n1 = 1;
			double n2 = 1.3;
			Vector normale_pour_transparence(N);
			if (dot(r.direction, N) > 0){ //On sort de la shere
				n1 = 1.3;
				n2 = 1;
				normale_pour_transparence = -N;
			}

			double radical = 1 - carre(n1 / n2) * (1 - carre(dot(normale_pour_transparence, r.direction)));
			if (radical > 0){
				Vector direction_refracte = (n1 / n2) * (r.direction - dot(r.direction, normale_pour_transparence) * normale_pour_transparence) - normale_pour_transparence * sqrt(radical);
				Ray rayon_refracte(P - 0.001 * normale_pour_transparence, direction_refracte);
				intensite_pixel = getColor(rayon_refracte, s, recursion + 1);
			}

		} else {
			Ray rayLight(P + 0.01*N , (s.position_lumiere - P).getNormalized()); // On envoi un rayon du point vers la lumiere
			// Au dessus on ajoute un truc à P pour décoller un peu l'origine du rayon et éviter le bruit sur l'image
			Vector N_light, P_light;
			int sphere_id_light;
			double t_light;
			bool intersection_light = s.intersection(rayLight, P_light, N_light, sphere_id_light, t_light); // On regarde si il y a une intersection entre le point qui va vers la lumiere et la lumiere
			// Pour au dessus, creer une structure de retour qui contient tout ces parametres, c'est beaucoup plus propre
			double d_light2 = (s.position_lumiere - P).getNorm2();

			

			if (intersection_light && ((t_light * t_light) < d_light2)){
				intensite_pixel = Vector(0, 0, 0);
			}
			else{
				intensite_pixel = s.spheres[id].albedo * s.intensite_lumiere * std::max(0., dot((s.position_lumiere - P).getNormalized(), N)) / d_light2;
			}
		}


	}

	return intensite_pixel;
}

int main(int argc, char const *argv[]){
	int H = 1024;
	int W = 1024;
	double fov = 90 * M_PI / 180; // angle de vue

	Sphere s1(Vector(30, 0, -55), 20, Vector(1, 0, 0), true, false); // (coord par rapport à la caméra) puis rayons
	Sphere s2(Vector(0, -2000-20, 0), 2000, Vector(1, 1, 1), false, false); // (coord par rapport à la caméra) puis rayons
	Sphere s3(Vector(0, 2000+100, 0), 2000, Vector(1, 1, 1), false, false); // (coord par rapport à la caméra) puis rayons
	Sphere s4(Vector(-2000-50, 0, 0), 2000, Vector(0, 1, 0), false, false); // (coord par rapport à la caméra) puis rayons
	Sphere s5(Vector(2000+50, 0, 0), 2000, Vector(0, 0, 1), false, false); // (coord par rapport à la caméra) puis rayons
	Sphere s6(Vector(0, 0, -2000-100), 2000, Vector(0, 1, 1), false, false); // (coord par rapport à la caméra) puis rayons
	Sphere s7(Vector(-30, 0, -55), 20, Vector(1, 0, 0), false, true); // (coord par rapport à la caméra) puis rayons

	Scene s;
	s.addSphere(s1);
	s.addSphere(s2);
	s.addSphere(s3);
	s.addSphere(s4);
	s.addSphere(s5);
	s.addSphere(s6);
	s.addSphere(s7);

	s.position_lumiere = Vector(15, 70, -30);//Notre lumiere

	s.intensite_lumiere = 1000000000;// plus la lumiere est intense plus ça va briller

	std::vector<unsigned char> image(W * H *3);

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector direction(j - W / 2, i - H / 2, - W / (2 * tan(fov / 2)));
			direction.normalize();

			Ray r(Vector(0, 0, 0), direction);
		
			Vector color = getColor(r, s);

			image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0], 1/2.2))); // rouge 
			image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1], 1/2.2))); // vert
			image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2],1/2.2))); // bleu
			// Le power permet de faire une correction gamma, ça permet de rendre les couleures claires plus basses car l'oeil humain y est plus sensible
		}
	}


	save_img("output.ppm", &image[0], W, H);

	std::cout << "hello" << std::endl;

	return 0;
}


