#include <stdio.h>
#include <iostream>
#include <vector>
#include <GL/glut.h>
#include "include/CLI11.hpp"
#include <fstream>
#include <string>
#include <jsoncpp/json/json.h>
#include "Raytracer.h"

int H = 1024;
int W = 1024;
unsigned char image[1024 * 1024 * 3];
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

void keyBoard(unsigned char key, int x, int y){
	SimpleRayTracer rayTracer3(2);
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
		rayTracer3.createImage(s, fov_cam, imageGL, H, W, ang, cam, light_position_value);
		glutPostRedisplay();
	}
}

void createWindow(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE);
	glutInitWindowSize(W, H);
	glutCreateWindow("RayTracer");
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


	FlatPaintingRayTracer rayTracer1;
	SimpleRayTracer rayTracer2(1);
	SimpleRayTracer rayTracer3(2);
	ComplexRayTracer rayTracer4;

	switch(number) {
		case 1:
			rayTracer1.createImage(s, fov_cam, image, H, W, ang, cam);
			break;
		case 2:
			createWindow(argc, argv);
			keyboard_help();
			rayTracer2.createImage(s, fov_cam, image, H, W, ang, cam, light_position_value);
			rayTracer3.createImage(s, fov_cam, imageGL, H, W, ang, cam, light_position_value);
			glutDisplayFunc(display);
			break;
		case 3:
			rayTracer4.createImage(s, fov_cam, image, H, W, ang, cam, pixelsampling);
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