#include <stdio.h>
#include <iostream>
#include <vector>
#include <GL/glut.h>
#include "include/CLI11.hpp"
#include "Vector.h"
#include <fstream>
#include <string>
#include <jsoncpp/json/json.h>

int H = 1024;
int W = 1024;
unsigned char image[1024*1024*3];
unsigned char imageGL[1024 * 1024 * 3];
double cam[3] = {0, 0, 0};
double ang[2] = {0, 0};
Scene s;
Vector light_position_value;
int fov_cam;

/*
void read_scene_file(const char* filename){
	std::ifstream file(filename);
	Json::Reader reader;
	Json::Value obj;
	reader.parse(file, obj);

	const Json::Value& spheres = obj["spheres"];
	const Json::Value& rectangles = obj["rectangles"];
	const Json::Value& triangles = obj["triangles"];
	const Json::Value& cylindres = obj["cylindres"];

	const Json::Value& intensite_lumiere = obj["intensite_lumiere"];
	const Json::Value& position_lum = obj["position_lumiere"];
	const Json::Value& fov = obj["fov"];

	for (int i = 0; i < spheres.size(); i++){
		bool mirror = (spheres[i]["mirror"].asInt() == 1) ? true: false; 
		bool glass = (spheres[i]["glass"].asInt() == 1) ? true: false;
		Vector axe(spheres[i]["axe"][0].asInt(), spheres[i]["axe"][1].asInt(), spheres[i]["axe"][2].asInt());
		int rayon = spheres[i]["rayon"].asInt();
		Vector color(spheres[i]["couleur"][0].asInt(), spheres[i]["couleur"][1].asInt(), spheres[i]["couleur"][2].asInt());
		double spec = spheres[i]["spec"].asDouble();
		Sphere* s5 = new Sphere(axe, rayon, color, mirror, glass, spec);
		s.addSphere(s5);
		if (i == 0){
			s.addSphere(s5);
			s.light = s5;
			continue;
		} 
	}
	
	s.light_intensity = intensite_lumiere.asInt();
	light_position = position_lum[0].asInt(), position_lum[1].asInt(), position_lum[2].asInt();
	fov_cam = fov.asInt();


	for (int i = 0; i < rectangles.size(); i++){
		Vector vec1(rectangles[i]["c1"][0].asInt(), rectangles[i]["c1"][1].asInt(), rectangles[i]["c1"][2].asInt());
		Vector vec2(rectangles[i]["c2"][0].asInt(), rectangles[i]["c2"][1].asInt(), rectangles[i]["c2"][2].asInt());
		Vector vec3(rectangles[i]["c3"][0].asInt(), rectangles[i]["c3"][1].asInt(), rectangles[i]["c3"][2].asInt());
		Vector vec4(rectangles[i]["c4"][0].asInt(), rectangles[i]["c4"][1].asInt(), rectangles[i]["c4"][2].asInt());
		Vector color(rectangles[i]["couleur"][0].asInt(), rectangles[i]["couleur"][1].asInt(), rectangles[i]["couleur"][2].asInt());
		bool mirror = rectangles[i]["mirror"].asInt() == 1;
		bool glass = rectangles[i]["glass"].asInt() == 1;
		Rectangle* s6 = new Rectangle(vec1, vec2, vec3, vec4, color, mirror, glass, rectangles[i]["spec"].asDouble());
		s.addRectangle(s6);
	}

	for (int i = 0; i < triangles.size(); i++){
		Vector vec1(triangles[i]["c1"][0].asInt(), triangles[i]["c1"][1].asInt(), triangles[i]["c1"][2].asInt());
		Vector vec2(triangles[i]["c2"][0].asInt(), triangles[i]["c2"][1].asInt(), triangles[i]["c2"][2].asInt());
		Vector vec3(triangles[i]["c3"][0].asInt(), triangles[i]["c3"][1].asInt(), triangles[i]["c3"][2].asInt());
		Vector color(triangles[i]["couleur"][0].asInt(), triangles[i]["couleur"][1].asInt(), triangles[i]["couleur"][2].asInt());
		bool mirror = triangles[i]["mirror"].asInt() == 1;
		bool glass = triangles[i]["glass"].asInt() == 1; 
		Triangle* s7 = new Triangle(vec1, vec2, vec3, color, mirror, glass, triangles[i]["spec"].asDouble());
		s.addTriangle(s7);
	}

	for (int i = 0; i < cylindres.size(); i++){ 
		bool mirror = cylindres[i]["mirror"].asInt() == 1;
		bool glass = cylindres[i]["glass"].asInt() == 1;
		Vector color(cylindres[i]["couleur"][0].asInt(), cylindres[i]["couleur"][1].asInt(), cylindres[i]["couleur"][2].asInt());
		Vector origin(cylindres[i]["axe"][0].asInt(), cylindres[i]["axe"][1].asInt(), cylindres[i]["axe"][2].asInt());
		Cylinder* s8 = new Cylinder(origin, cylindres[i]["rayon"].asInt(), cylindres[i]["hauteur"].asInt(), color, mirror, glass, cylindres[i]["spec"].asDouble());
		s.addCylinder(s8);
	}
	
	return;
}*/

void get_rectangles_from_file(Json::Value obj){
	const Json::Value& rectangles = obj["rectangles"];
	for (int i = 0; i < rectangles.size(); i++){
		Vector vec1(rectangles[i]["c1"][0].asInt(), rectangles[i]["c1"][1].asInt(), rectangles[i]["c1"][2].asInt());
		Vector vec2(rectangles[i]["c2"][0].asInt(), rectangles[i]["c2"][1].asInt(), rectangles[i]["c2"][2].asInt());
		Vector vec3(rectangles[i]["c3"][0].asInt(), rectangles[i]["c3"][1].asInt(), rectangles[i]["c3"][2].asInt());
		Vector vec4(rectangles[i]["c4"][0].asInt(), rectangles[i]["c4"][1].asInt(), rectangles[i]["c4"][2].asInt());
		Vector color(rectangles[i]["color"][0].asInt(), rectangles[i]["color"][1].asInt(), rectangles[i]["color"][2].asInt());
		bool mirror = rectangles[i]["mirror"].asInt() == 1;
		bool glass = rectangles[i]["glass"].asInt() == 1;
		Rectangle* s6 = new Rectangle(vec1, vec2, vec3, vec4, color, mirror, glass, rectangles[i]["spec"].asDouble());
		s.addRectangle(s6);
	}
}

void get_triangles_from_file(Json::Value obj){
	const Json::Value& triangles = obj["triangles"];
	for (int i = 0; i < triangles.size(); i++){
		Vector vec1(triangles[i]["c1"][0].asInt(), triangles[i]["c1"][1].asInt(), triangles[i]["c1"][2].asInt());
		Vector vec2(triangles[i]["c2"][0].asInt(), triangles[i]["c2"][1].asInt(), triangles[i]["c2"][2].asInt());
		Vector vec3(triangles[i]["c3"][0].asInt(), triangles[i]["c3"][1].asInt(), triangles[i]["c3"][2].asInt());
		Vector color(triangles[i]["color"][0].asInt(), triangles[i]["color"][1].asInt(), triangles[i]["color"][2].asInt());
		bool mirror = triangles[i]["mirror"].asInt() == 1;
		bool glass = triangles[i]["glass"].asInt() == 1; 
		Triangle* s7 = new Triangle(vec1, vec2, vec3, color, mirror, glass, triangles[i]["spec"].asDouble());
		s.addTriangle(s7);
	}
}

void get_cylinders_from_file(Json::Value obj){
	const Json::Value& cylinders = obj["cylinders"];
	for (int i = 0; i < cylinders.size(); i++){ 
		bool mirror = cylinders[i]["mirror"].asInt() == 1;
		bool glass = cylinders[i]["glass"].asInt() == 1;
		Vector color(cylinders[i]["color"][0].asInt(), cylinders[i]["color"][1].asInt(), cylinders[i]["color"][2].asInt());
		Vector origin(cylinders[i]["axis"][0].asInt(), cylinders[i]["axis"][1].asInt(), cylinders[i]["axis"][2].asInt());
		Cylinder* s8 = new Cylinder(origin, cylinders[i]["rayon"].asInt(), cylinders[i]["height"].asInt(), color, mirror, glass, cylinders[i]["spec"].asDouble());
		s.addCylinder(s8);
	}
}

void get_spheres_from_file(Json::Value obj){
	const Json::Value& spheres = obj["spheres"];
	for (int i = 0; i < spheres.size(); i++){
		bool mirror = (spheres[i]["mirror"].asInt() == 1) ? true: false; 
		bool glass = (spheres[i]["glass"].asInt() == 1) ? true: false;
		Vector axis(spheres[i]["axis"][0].asInt(), spheres[i]["axis"][1].asInt(), spheres[i]["axis"][2].asInt());
		int rayon = spheres[i]["rayon"].asInt();
		Vector color(spheres[i]["color"][0].asInt(), spheres[i]["color"][1].asInt(), spheres[i]["color"][2].asInt());
		double spec = spheres[i]["spec"].asDouble();
		Sphere* s5 = new Sphere(axis, rayon, color, mirror, glass, spec);
		s.addSphere(s5);
		if (i == 0){
			s.addSphere(s5);
			s.light = s5;
			continue;
		} 
	}
}

void get_light_and_cam_from_file(Json::Value obj){
	const Json::Value& light_intensity_value = obj["light_intensity"];
	const Json::Value& light_position = obj["light_position"];
	const Json::Value& fov = obj["fov"];
	s.light_intensity = light_intensity_value.asDouble();
	light_position_value = Vector(light_position[0].asInt(), light_position[1].asInt(), light_position[2].asInt());
	fov_cam = fov.asInt();
}
void read_scene_file(const char* filename){
	std::ifstream file(filename);
	Json::Reader reader;
	Json::Value obj;
	reader.parse(file, obj);
	get_spheres_from_file(obj);
	get_cylinders_from_file(obj);
	get_triangles_from_file(obj);
	get_rectangles_from_file(obj);
	get_light_and_cam_from_file(obj);
	/*light_position_value = Vector(15, 60, -40);
	s.light_intensity = 100000000000;*/

}

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

 double carre(double n){
 	return n * n;
 }

double Phong_BRDF(const Vector& wi, const Vector& wo, const Vector& N, double phong_expo){
	Vector reflechi = wo.reflect(N);
	double lobe = std::pow(dot(reflechi, wi), phong_expo) * (phong_expo + 2) / (2. * M_PI);
	return lobe;
}

Vector getColor(const Ray& r, const Scene& s, int recursion = 0);

// const enmpeche de faire une copie locale
// On fait une fonction recursive pour calculer la reflexion des surfaces
// La fonction renvoie la couleur du pixel obtenu en envoyant un rayon R dans la scene S
Vector getColor(const Ray& r, const Scene& s, int recursion){
	if (recursion >= 5){
		return Vector(0, 0, 0);
	}

	int id;
	double t;
	Vector P, N;
	bool has_inter = s.intersection(r, P, N, id, t);
	Vector intensite_pixel(0,0,0); // trois cmposantes du coup on prend un vecteur (A CHANGER)

	if (has_inter) {
		//On check  la nature de la forme que croise le rayon
		if (s.shapes[id] -> is_mirror){
			Vector direction_mirroir = r.direction.reflect(N);
			Ray rayon_mirroir(P + 0.001 * N, direction_mirroir);
			intensite_pixel = getColor(rayon_mirroir, s, recursion + 1);
		} 
		//calcule intensité pixel si transparant
		else if (s.shapes[id] -> is_transparent){
			double n1 = 1;
			double n2 = 1.3;
			Vector normale_pour_transparence(N);
			Ray new_ray;
			bool entering = true;
			if (dot(r.direction, N) > 0){ //On sort de l'objet
				n1 = 1.3;
				n2 = 1;
				normale_pour_transparence = -N;
				entering = false;
			}

			double radical = 1 - carre(n1 / n2) * (1 - carre(dot(normale_pour_transparence, r.direction)));
			if (radical > 0){
				Vector direction_refracte = (n1 / n2) * (r.direction - dot(r.direction, normale_pour_transparence) * normale_pour_transparence) - normale_pour_transparence * sqrt(radical);
				Ray rayon_refracte(P - 0.001 * normale_pour_transparence, direction_refracte);
				
				// Utilisation du coef de shnik
				double R0 = carre((n1 - n2) / (n1 + n2));
				double R;
				if (entering){
					R = R0 + (1 - R0) * std::pow(1 + dot(r.direction, N), 5);
				}
				else {
					R = R0 + (1 - R0) * std::pow(1 - dot(direction_refracte, N), 5);
				}

				if (generator_rand(generator) < R){
					new_ray = Ray(P + 0.001 * normale_pour_transparence, r.direction.reflect(N));

				}
				else{
					new_ray = Ray(P - 0.001 * normale_pour_transparence, direction_refracte);
				}
				intensite_pixel = getColor(new_ray, s, recursion + 1);

			}
			// cas d'une reflection totale
			else {
				new_ray = Ray(P + 0.001 * normale_pour_transparence, r.direction.reflect(N));
			}

		} 
		// Ajout éclairage direct
		//ni transparent ni miroir
		else {
			// On envoi un rayon vers la sphere
			Vector axeOP = (P - s.light-> O).getNormalizedCopy();
			Vector dir_aleatoire = random_cos(axeOP);
			Vector point_aleatoire = dir_aleatoire * s.light -> R + s.light -> O;
			Vector wi = (point_aleatoire - P).getNormalizedCopy();
			double d_light2 = (point_aleatoire - P).getNorm();
			Vector Np = dir_aleatoire; // meme chose

			Ray rayLight(P + 0.01*N , wi); // On envoi un rayon du point vers la lumiere
			// Au dessus on ajoute un truc à P pour décoller un peu l'origine du rayon et éviter le bruit sur l'image
			Vector N_light, P_light;
			int sphere_id_light;
			double t_light;
			bool intersection_light = s.intersection(rayLight, P_light, N_light, sphere_id_light, t_light); // On regarde si il y a une intersection entre le point qui va vers la lumiere et la lumiere

			// On relance un rayon vers la lumierre pour savoir s'il y a intersection, on ajoute 0.99 car la lumiere est une sphere
			if (intersection_light && ((t_light * t_light) < d_light2 * 0.99)){
				intensite_pixel = Vector(0, 0, 0); 
			}
			else{
				intensite_pixel = (s.light_intensity / (4 * M_PI * d_light2) * std::max(0., dot(N, wi)) * dot(Np, -wi) / dot(axeOP, dir_aleatoire)) * (M_PI) * ((1. - s.shapes[id] -> spec) * s.shapes[id] -> albedo / M_PI + Phong_BRDF(wi, r.direction, N, s.shapes[id] -> phong_expo) * s.shapes[id] -> spec * s.shapes[id] -> albedo);
			}

		// Ajout éclairage indirect 
			Vector random_direction;
			double p = 1 - s.shapes[id] -> spec;
			bool sample_diffuse;
			Vector R = r.direction.reflect(N);
			// On genere soit un echantillon diffut
			if (generator_rand(generator) < p){
				sample_diffuse = true;
				random_direction = random_cos(N);
			}
			// Soit une echantillon speculaire
			else { 
				sample_diffuse = false;
				random_direction = random_phong(R, s.shapes[id] -> phong_expo);
				if ((dot(random_direction, N) < 0) || (dot(random_direction, R) < 0)) {
					return Vector(0., 0., 0.);
				}
			}
			Ray random_ray(P + 0.001 * N, random_direction);

			double proba_phong = (s.shapes[id] -> phong_expo + 1) / (2. * M_PI) * std::pow(dot(R, random_direction), s.shapes[id] -> phong_expo);
			double proba_globale = p * dot(N, random_direction) / (2. * M_PI) + (1. - p) * proba_phong;
			
			if (sample_diffuse){
				intensite_pixel = intensite_pixel + getColor(random_ray, s, recursion + 1) * s.shapes[id] -> albedo * dot(N, random_direction) / M_PI / proba_globale;
			}
			else{
				intensite_pixel = intensite_pixel + getColor(random_ray, s, recursion + 1) * dot(N, random_direction) * Phong_BRDF(random_direction, r.direction, N, s.shapes[id] -> phong_expo) * s.shapes[id] -> spec * s.shapes[id] -> albedo / proba_globale;
			}
		}
	}
	return intensite_pixel;
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
	for (int i = 0; i < nb_sampling; i++){
		color = color + (getColor(r, s) / nb_sampling);
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
	
	read_scene_file(file.c_str());

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