#include <stdio.h>
#include <iostream>
#include <vector>
#include <GL/glut.h>
#include "include/CLI11.hpp"
#include <fstream>
#include <string>
#include <jsoncpp/json/json.h>
#include "Color.h"

int H = 1024;
int W = 1024;
unsigned char image[1024*1024*3];
unsigned char imageGL[1024 * 1024 * 3];
double cam[3] = {0, 0, 0};
double ang[2] = {0, 0};

Scene s;
Vector light_position_value;
int fov_cam = 15;

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

void display() {
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(W,H,GL_RGB,GL_UNSIGNED_BYTE, imageGL);
    glFlush();
}

void createImage1(int i, int j, double fov, bool has_inter, int id){
	Vector pixel_intensity(0,0,0); // trois cmposantes du coup on prend un vecteur (A CHANGER)
	if (has_inter) {
		pixel_intensity = s.shapes[id]->albedo;
	}
	image[((H - i - 1) * W + j) * 3 + 0] = pixel_intensity[0] * 255; // rouge 
	image[((H - i - 1) * W + j) * 3 + 1] = pixel_intensity[1] * 255; // vert
	image[((H - i - 1) * W + j) * 3 + 2] = pixel_intensity[2] * 255; // bleu
}

void createImage2(int i, int j, double fov, int light_intensity, bool has_inter, int id, Vector& P, Vector& N) {
	Vector pixel_intensity = 0; // trois cmposantes du coup on prend un vecteur (A CHANGER)
	if (has_inter) {
		pixel_intensity = s.shapes[id]->albedo * (light_intensity * std::max(0., dot((light_position_value - P).getNormalizedCopy(), N))) / (light_position_value - P).getNorm();
	}
	image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., pixel_intensity[0])); // rouge 
	image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., pixel_intensity[1])); // vert
	image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., pixel_intensity[2])); // bleu
}

void createImageGL(int i, int j, double fov, int light_intensity, bool has_inter, int id, Vector& P, Vector& N){
	Vector pixel_intensity = 0; // trois cmposantes du coup on prend un vecteur (A CHANGER)
	if (has_inter) {
		pixel_intensity = s.shapes[id]->albedo * (light_intensity * std::max(0., dot((light_position_value - P).getNormalizedCopy(), N))) / (light_position_value - P).getNorm();
	}
	imageGL[((H - (H - i - 1) - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., pixel_intensity[0])); // rouge 
	imageGL[((H - (H - i - 1) - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., pixel_intensity[1])); // vert
	imageGL[((H - (H - i - 1) - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., pixel_intensity[2])); // bleu
}

void createImage3(int i, int j, double fov, double nb_sampling){
	// methode de Box Muller (anti-alliasing) ->PIXEL SAMPLING
	// créé un nombre aléatoire qui suit une gaussienne, ça permet d'envoyer un rayon pas forcement au centre
	double r1 = generator_rand(generator);
	double r2 = generator_rand(generator);
	double dx = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
	double dy = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2);
	Vector direction(j - W / 2 + 0.5 * dx , i - H / 2 + 0.5 * dy, - W / (2 * tan(fov / 2)));
	direction.normalize();
	Ray r(Vector(0, 0, 0), direction);
	Vector color(0., 0., 0.);
	Color colorPixel;
	for (int i = 0; i < nb_sampling; i++){
		color = color + (colorPixel.getPixelColor(r, s) / nb_sampling);
	} 
	image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0], 1 / 2.2))); // rouge 
	image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1], 1 / 2.2))); // vert
	image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2], 1 / 2.2))); // bleu
	// Le power permet de faire une correction gamma, ça permet de rendre les couleures claires plus basses car l'oeil humain y est plus sensible
}

void intersect(int option, double pixelsampling) {
	double fov = fov_cam * M_PI / 180; // angle de vue
	int light_intensity = s.light_intensity / 1000;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector direction((j - W / 2) + ang[0], (i - H / 2) + ang[1], - W / (2 * tan(fov / 2)));
			direction.normalize();
 			Ray r(Vector(cam[0], cam[1], cam[2]), direction);
			int id;
			Vector P, N;
			double min;
			bool has_inter = s.intersection(r, P, N, id, min);
			switch(option) {
				case 1:
					createImage1(i, j, fov, has_inter, id);
					break;
				case 2:
					createImage2(i, j, fov, light_intensity, has_inter, id, P, N);
					break;
				case 3:
					createImageGL(i, j, fov, light_intensity, has_inter, id, P, N);
					break;
				case 4:
					createImage3(i, j, fov, pixelsampling);
					break;
				default:
					return;
			}
		}
	}
	if(option == 3){
		glutDisplayFunc(display);
	}
}

void keyBoard(unsigned char key, int x, int y){
	bool isMoving = false;
	switch(key){
		case 'z':
		case 'Z':
			isMoving = true;
			cam[1] += 10;
			break;
		case 'q':
		case 'Q':
			isMoving = true;
			cam[0] += -10;
			break;
		case 's':
		case 'S':
			isMoving = true;
			cam[1] += -10;
			break;
		case 'd':
		case 'D':
			isMoving = true;
			cam[0] += 10;
			break;
		case 'a':
		case 'A':
			isMoving = true;
			cam[2] += -10;
			break;
		case 'e':
		case 'E':
			isMoving = true;
			cam[2] += 10;
			break;
		case 'j':
		case 'J':
			isMoving = true;
			ang[0] += -100;
			break;
		case 'l':
		case 'L':
			isMoving = true;
			ang[0] += 100;
			break;
		case 'i':
		case 'I':
			isMoving = true;
			ang[1] += 100;
			break;
		case 'k':
		case 'K':
			isMoving = true;
			ang[1] += -100;
			break;
		default:
			break;
	}
	if(isMoving){
		intersect(3, 1);
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

void keyboard_help() {
	std::cout << "you can use the keyboard keys to move around the image. Keys are as follows: " << std::endl;
	std::cout << "\t z, move up" << std::endl;
	std::cout << "\t s, move down" << std::endl;
	std::cout << "\t d, move right" << std::endl;
	std::cout << "\t q, move left" << std::endl;
	std::cout << "\t a, move forward" << std::endl;
	std::cout << "\t e, move backward" << std::endl;
	std::cout << "Following keys allow to change camera direction: " << std::endl;
	std::cout << "\t i, look up" << std::endl;
	std::cout << "\t k, look down" << std::endl;
	std::cout << "\t l, look right" << std::endl;
	std::cout << "\t j, look left" << std::endl;
}

int main(int argc, char *argv[]){
	CLI::App app{"In computer graphics, ray tracing is a rendering technique for generating an image by tracing the path of light as pixels in an image plane and simulating the effects of its encounters with virtual objects.\n This project allows to create scenes using raytracing."};
	int number = 1;
    app.add_option("-n,--number", number, "Wich level number of the project you want to use. Must be between 1 and 3.");
	double pixelsampling = 1;
	app.add_option("-p, --pixelsampling", pixelsampling, "Number of pixel sampling that you want.");
	std::string file = "defaut";
	app.add_option("-i,--input", file, "File name that you want to use. This file needs to be in json format.");
	std::string outputName;
	app.add_option("-o,--output", outputName, "File name of the output image")->required();
	CLI11_PARSE(app, argc, argv);
	
	Reader r;
	r.read_scene_file(file.c_str(), s, light_position_value, fov_cam);

	switch(number) {
		case 1:
			intersect(1, pixelsampling);
			break;
		case 2:
			createWindow(argc, argv);
			keyboard_help();
			intersect(2, pixelsampling);
			intersect(3, pixelsampling);
			break;
		case 3:
			intersect(4, pixelsampling);
			break;
		default:
			std::cout << "Number " << number << " isn't a valid number. Must be 1, 2, or 3 \n";
			return 0;
	}
	
	save_img(outputName.c_str(), &image[0], W, H);
	std::cout << "image " << outputName.c_str() << " created" << std::endl;

	if(number == 2){
		glutKeyboardFunc(keyBoard);
		glutMainLoop();
	}

	return 0;
}