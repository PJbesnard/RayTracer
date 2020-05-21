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

int main(int argc, char const *argv[]){
	int H = 1024;
	int W = 1024;
	double fov = 90 * M_PI / 180; // angle de vue

	//Sphere s1(Vector(0, 0, -55), 20, Vector(1, 0, 0)); // (coord par rapport à la caméra) puis rayons
	Sphere s2(Vector(0, -2000-20, 0), 2000, Vector(1, 1, 1)); // (coord par rapport à la caméra) puis rayons
	Sphere s3(Vector(0, 2000+100, 0), 2000, Vector(1, 1, 1)); // (coord par rapport à la caméra) puis rayons
	Sphere s4(Vector(-2000-50, 0, 0), 2000, Vector(0, 1, 0)); // (coord par rapport à la caméra) puis rayons
	Sphere s5(Vector(2000+50, 0, 0), 2000, Vector(0, 0, 1)); // (coord par rapport à la caméra) puis rayons
	Sphere s6(Vector(0, 0, -2000-100), 2000, Vector(0, 1, 1)); // (coord par rapport à la caméra) puis rayons
	Sphere s7(Vector(0, 2000+50, 0), 2000, Vector(1, 0, 1)); // (coord par rapport à la caméra) puis rayons

	Scene s;
	//s.addSphere(s1);
	s.addSphere(s2);
	s.addSphere(s3);
	s.addSphere(s4);
	s.addSphere(s5);
	s.addSphere(s6);
	s.addSphere(s7);

	Triangle t(Vector(60, 0, -70), Vector(30, 20, -80), Vector(10, 10, -100), Vector(0, 0, 1));
	s.addTriangle(t);

	Cylindre c(Vector(20, -20, -55), 5, 10, Vector(1, 0, 0));
	s.addCylindre(c);

	Vector position_lumiere(15, 60, -40); //Notre lumiere

	double intensite_lumiere = 1000000; // plus la lumiere est intense plus ça va briller

	std::vector<unsigned char> image(W * H *3);

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector direction(j - W / 2, i - H / 2, - W / (2 * tan(fov / 2)));
			direction.normalize();

			Ray r(Vector(0, 0, 0), direction);
			int id;
			Vector P, N;
			bool has_inter = s.intersection(r, P, N, id);

			Vector intensite_pixel = 0; // trois cmposantes du coup on prend un vecteur (A CHANGER)
			if (has_inter) {
				intensite_pixel = s.objects[id]->albedo * (intensite_lumiere * std::max(0., dot((position_lumiere - P).getNormalized(), N))) / (position_lumiere - P).getNorm2();
			}

			image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel[0])); // rouge 
			image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel[1])); // vert
			image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel[2])); // bleu
		}
	}


	save_img("output.ppm", &image[0], W, H);

	std::cout << "hello" << std::endl;

	return 0;
}


