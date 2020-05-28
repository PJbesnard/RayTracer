#include <stdio.h>
#include <iostream>
#include <vector>
#include <GL/glut.h>
#include "include/CLI11.hpp"
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

int H = 1024;
int W = 1024;
//std::vector<unsigned char> image(W * H *3);
unsigned char image[1024*1024*3];
double cam[3] = {0, 0, 0};
double ang[2] = {0, 0};
Scene s;

void display() {
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(W,H,GL_RGB,GL_UNSIGNED_BYTE, image);
	glutSwapBuffers();
}

void intersect(){
	double fov = 90 * M_PI / 180; // angle de vue

	Vector position_lumiere(15, 60, -40); //Notre lumiere

	double intensite_lumiere = 1000000; // plus la lumiere est intense plus ça va briller

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector direction((j - W / 2) + ang[0], (i - H / 2) + ang[1] , - W / (2 * tan(fov / 2)));
			direction.normalize();

			//Ray r(Vector(0, 0, 0), direction);
			Ray r(Vector(cam[0], cam[1], cam[2]), direction);
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
	glutDisplayFunc(display);
}

void intersect1(){
	double fov = 90 * M_PI / 180; // angle de vue

	Vector position_lumiere(15, 60, -40); //Notre lumiere

	double intensite_lumiere = 1000000; // plus la lumiere est intense plus ça va briller

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector direction((j - W / 2), (i - H / 2), - W / (2 * tan(fov / 2)));
			direction.normalize();

			Ray r(Vector(0, 0, 0), direction);
			int id;
			Vector P, N;
			bool has_inter = s.intersection(r, P, N, id);

			Vector intensite_pixel(0,0,0); // trois cmposantes du coup on prend un vecteur (A CHANGER)
			if (has_inter) {
				intensite_pixel = s.objects[id]->albedo;
			}

			image[((H - i - 1) * W + j) * 3 + 0] = intensite_pixel[0] * 255; // rouge 
			image[((H - i - 1) * W + j) * 3 + 1] = intensite_pixel[1] * 255; // vert
			image[((H - i - 1) * W + j) * 3 + 2] = intensite_pixel[2] * 255; // bleu
		}
	}
}

void keyBoard(unsigned char key, int x, int y){
	bool isMoving = false;
	switch(key){
		case 'z':
		case 'Z':
			std::cout << "on bouge au dessus" << std::endl;
			isMoving = true;
			cam[1] += -10;
			break;

		case 'q':
		case 'Q':
			std::cout << "on bouge à gauche" << std::endl;
			isMoving = true;
			cam[0] += -10;
			break;

		case 's':
		case 'S':
			std::cout << "on bouge en bas" << std::endl;
			isMoving = true;
			cam[1] += 10;
			break;

		case 'd':
		case 'D':
			std::cout << "on bouge à droite " << cam[0] << " " << cam[1] << " " << cam[2] << std::endl;
			isMoving = true;
			cam[0] += 10;
			break;

		case 'a':
		case 'A':
			std::cout << "on recule " << cam[0] << " " << cam[1] << " " << cam[2] << std::endl;
			isMoving = true;
			cam[2] += -10;
			break;

		case 'e':
		case 'E':
			std::cout << "on avance " << cam[0] << " " << cam[1] << " " << cam[2] << std::endl;
			isMoving = true;
			cam[2] += 10;
			break;

		case 'j':
		case 'J':
			std::cout << "on tourne l'angle a gauche " << ang[0] << " " << ang[1] << " " << ang[2] << std::endl;
			isMoving = true;
			ang[0] += -100;
			break;

		case 'l':
		case 'L':
			std::cout << "on tourne l'angle a droite " << ang[0] << " " << ang[1] << " " << ang[2] << std::endl;
			isMoving = true;
			ang[0] += 100;
			break;

		case 'i':
		case 'I':
			std::cout << "on leve la tete " << ang[0] << " " << ang[1] << " " << ang[2] << std::endl;
			isMoving = true;
			ang[1] += -100;
			break;

		case 'k':
		case 'K':
			std::cout << "on baisse la tete " << ang[0] << " " << ang[1] << " " << ang[2] << std::endl;
			isMoving = true;
			ang[1] += 100;
			break;

		default:
			break;
	}
	if(isMoving){
		intersect();
		glutPostRedisplay();
	}
}

void createWindow(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE);
	glutInitWindowSize(W, H);
	glutCreateWindow("first");
	gluLookAt(-100., -200., -100., -100., -100., -100., -100., -51., -100.);
}

int main(int argc, char *argv[]){
	CLI::App app{"In computer graphics, ray tracing is a rendering technique for generating an image by tracing the path of light as pixels in an image plane and simulating the effects of its encounters with virtual objects.\n This project allows to create scenes using raytracing."};
	int number = 1;
    app.add_option("-n,--number", number, "Wich level number of the project you want to use. Must be between 1 and 3.");
	int pixelsampling;
	app.add_option("-p, --pixelsampling", pixelsampling, "Number of pixel sampling that you want.");
	std::string file = "defaut";
	app.add_option("-i,--input", file, "File name that you want to use. This file needs to be in json format.");
	std::string outputName;
	app.add_option("-o,--output", outputName, "File name of the output image")->required();
	CLI11_PARSE(app, argc, argv);
	
	Sphere s1(Vector(0, 0, -55), 20, Vector(1, 0, 0)); // (coord par rapport à la caméra) puis rayons
	Sphere s2(Vector(0, -2000-20, 0), 2000, Vector(1, 1, 1)); // (coord par rapport à la caméra) puis rayons
	Sphere s3(Vector(0, 2000+100, 0), 2000, Vector(1, 1, 1)); // (coord par rapport à la caméra) puis rayons
	Sphere s4(Vector(-2000-50, 0, 0), 2000, Vector(0, 1, 0)); // (coord par rapport à la caméra) puis rayons
	Sphere s5(Vector(2000+50, 0, 0), 2000, Vector(0, 0, 1)); // (coord par rapport à la caméra) puis rayons
	Sphere s6(Vector(0, 0, -2000-100), 2000, Vector(0, 1, 1)); // (coord par rapport à la caméra) puis rayons
	Sphere s7(Vector(0, 2000+50, 0), 2000, Vector(1, 0, 1)); // (coord par rapport à la caméra) puis rayons

	//s.addSphere(s1);
	s.addSphere(s2);
	s.addSphere(s3);
	s.addSphere(s4);
	s.addSphere(s5);
	s.addSphere(s6);
	s.addSphere(s7);

	Triangle t(Vector(40, 0, -70), Vector(40, 40, -70), Vector(10, 30, -70), Vector(0, 0, 1));
	//s.addTriangle(t);

	Cylindre c(Vector(-20, 10, -55), 20, 100, Vector(1, 0, 0));
	s.addCylindre(c);

	Rectangle r(Vector(10, 20, -70), Vector(20, 10, -70), Vector(10, 30, -70), Vector(20, 30, -70), Vector(1, 0, 0));
	s.addRectangle(r);

	switch(number) {
		case 1:
			intersect1();
			std::cout << "image " << outputName.c_str() << " created" << std::endl;
			break;
		case 2:
			//TODO
			break;
		case 3:
			//TODO
			break;
		default:
			std::cout << "Number " << number << " isn't a valid number. Must be 1, 2, or 3 \n";
			return 0;
	}
	
	save_img(outputName.c_str(), &image[0], W, H);

	/*
    //glutDisplayFunc(display);
	glutKeyboardFunc(keyBoard);
	glutMainLoop();*/

	return 0;
}