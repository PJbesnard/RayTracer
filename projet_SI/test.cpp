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
//std::vector<unsigned char> image(W * H *3);
unsigned char image[1024*1024*3];
unsigned char imageGL[1024 * 1024 * 3];
double cam[3] = {0, 0, 0};
double ang[2] = {0, 0};
Scene s;
Vector position_lumiere;
int fov_cam;


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
		double ks = spheres[i]["spec"].asDouble();
		Sphere* s5 = new Sphere(axe, rayon, color, mirror, glass, ks);
		s.addSphere(s5);
		if (i == 0){
			s.addSphere(s5);
			s.lumiere = s5;
			continue;
		} 
	}
	
	s.intensite_lumiere = intensite_lumiere.asInt();
	position_lumiere = position_lum[0].asInt(), position_lum[1].asInt(), position_lum[2].asInt();
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
		Cylindre* s8 = new Cylindre(origin, cylindres[i]["rayon"].asInt(), cylindres[i]["hauteur"].asInt(), color, mirror, glass, cylindres[i]["spec"].asDouble());
		s.addCylindre(s8);
	}
	
	return;

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
	/*
    unsigned char image2[1024 * 1024 * 3];
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            image2[((H - i - 1) * W + j) * 3 + 0] = image[((H - (H - i - 1) - 1) * W + j) * 3 + 0]; // rouge 
            image2[((H - i - 1) * W + j) * 3 + 1] = image[((H - (H - i - 1) - 1) * W + j) * 3 + 1]; // vert
            image2[((H - i - 1) * W + j) * 3 + 2] = image[((H - (H - i - 1) - 1) * W + j) * 3 + 2]; // bleu
        }
    }*/
    glDrawPixels(W,H,GL_RGB,GL_UNSIGNED_BYTE, imageGL);
    glFlush();
}

 double carre(double n){
 	return n * n;
 }


double Phong_BRDF(const Vector& wi, const Vector& wo, const Vector& N, double phong_exposant){
	Vector reflechi = wo.reflect(N);
	double lobe = std::pow(dot(reflechi, wi), phong_exposant) * (phong_exposant + 2) / (2. * M_PI);
	return lobe;
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

		//On check  la nature de la forme que croise le rayon
		if (s.objects[id] -> is_mirror){
			Vector direction_mirroir = r.direction.reflect(N);
			Ray rayon_mirroir(P + 0.001 * N, direction_mirroir);

			intensite_pixel = getColor(rayon_mirroir, s, recursion + 1);
		} 
		else if (s.objects[id] -> is_transparent){
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

				if (uniform(generator) < R){
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
		else {
			// On envoi un rayon vers la sphere
			Vector axeOP = (P - s.lumiere-> O).getNormalized();
			Vector dir_aleatoire = random_cos(axeOP);
			Vector point_aleatoire = dir_aleatoire * s.lumiere -> R + s.lumiere -> O;
			Vector wi = (point_aleatoire - P).getNormalized();
			double d_light2 = (point_aleatoire - P).getNorm2();
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
				intensite_pixel = (s.intensite_lumiere / (4 * M_PI * d_light2) * std::max(0., dot(N, wi)) * dot(Np, -wi) / dot(axeOP, dir_aleatoire)) * (M_PI) * ((1. - s.objects[id] -> ks) * s.objects[id] -> albedo / M_PI + Phong_BRDF(wi, r.direction, N, s.objects[id] -> phong_exposant) * s.objects[id] -> ks * s.objects[id] -> albedo);
			}

		// Ajout éclairage indirect 
			Vector random_direction;
			double p = 1 - s.objects[id] -> ks;
			bool sample_diffuse;
			Vector R = r.direction.reflect(N);
			// On genere soit un echantillon diffut
			if (uniform(generator) < p){
				sample_diffuse = true;
				random_direction = random_cos(N);
			}
			// Soit une echantillon speculaire
			else { 
				sample_diffuse = false;
				random_direction = random_phong(R, s.objects[id] -> phong_exposant);
				if ((dot(random_direction, N) < 0) || (dot(random_direction, R) < 0)) {
					return Vector(0., 0., 0.);
				}
			}
			Ray random_ray(P + 0.001 * N, random_direction);

			double proba_phong = (s.objects[id] -> phong_exposant + 1) / (2. * M_PI) * std::pow(dot(R, random_direction), s.objects[id] -> phong_exposant);
			double proba_globale = p * dot(N, random_direction) / (2. * M_PI) + (1. - p) * proba_phong;
			
			if (sample_diffuse){
				intensite_pixel = intensite_pixel + getColor(random_ray, s, recursion + 1) * s.objects[id] -> albedo * dot(N, random_direction) / M_PI / proba_globale;
			}
			else{
				intensite_pixel = intensite_pixel + getColor(random_ray, s, recursion + 1) * dot(N, random_direction) * Phong_BRDF(random_direction, r.direction, N, s.objects[id] -> phong_exposant) * s.objects[id] -> ks * s.objects[id] -> albedo / proba_globale;
			}
			
		}


	}

	return intensite_pixel;
}

void intersect3(int pixelsampling){
	int nb_sampling = pixelsampling;
	double fov = fov_cam * M_PI / 180; // angle de vue
#pragma cmp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			// methode de Box Muller (anti-alliasing) ->PIXEL SAMPLING
			// créé un nombre aléatoire qui suit une gaussienne, ça permet d'envoyer un rayon pas forcement au centre
			double r1 = uniform(generator);
			double r2 = uniform(generator);
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
	}

}

void intersect2(){
	double fov = fov_cam * M_PI / 180; // angle de vue
	int intensite_lum = s.intensite_lumiere / 1000;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector direction((j - W / 2) + ang[0], (i - H / 2) + ang[1] , - W / (2 * tan(fov / 2)));
			direction.normalize();

			//Ray r(Vector(0, 0, 0), direction);
			Ray r(Vector(cam[0], cam[1], cam[2]), direction);
			int id;
			Vector P, N;
			double min;
			bool has_inter = s.intersection(r, P, N, id, min);

			Vector intensite_pixel = 0; // trois cmposantes du coup on prend un vecteur (A CHANGER)
			if (has_inter) {
				intensite_pixel = s.objects[id]->albedo * (intensite_lum * std::max(0., dot((position_lumiere - P).getNormalized(), N))) / (position_lumiere - P).getNorm2();
			}

			image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel[0])); // rouge 
			image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel[1])); // vert
			image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel[2])); // bleu
		}
	}
}

void intersect2GL(){
	double fov = fov_cam * M_PI / 180; // angle de vue
	int intensite_lum = s.intensite_lumiere / 1000;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector direction((j - W / 2) + ang[0], (i - H / 2) + ang[1] , - W / (2 * tan(fov / 2)));
			direction.normalize();

			Ray r(Vector(cam[0], cam[1], cam[2]), direction);
			int id;
			Vector P, N;
			double min;
			bool has_inter = s.intersection(r, P, N, id, min);

			Vector intensite_pixel = 0; // trois cmposantes du coup on prend un vecteur (A CHANGER)
			if (has_inter) {
				intensite_pixel = s.objects[id]->albedo * (intensite_lum * std::max(0., dot((position_lumiere - P).getNormalized(), N))) / (position_lumiere - P).getNorm2();
			}

			imageGL[((H - (H - i - 1) - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., intensite_pixel[0])); // rouge 
			imageGL[((H - (H - i - 1) - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., intensite_pixel[1])); // vert
			imageGL[((H - (H - i - 1) - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., intensite_pixel[2])); // bleu
		}
	}
	glutDisplayFunc(display);
}

void intersect1(){
	double fov = fov_cam * M_PI / 180; // angle de vue

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector direction((j - W / 2), (i - H / 2), - W / (2 * tan(fov / 2)));
			direction.normalize();

			Ray r(Vector(0, 0, 0), direction);
			int id;
			Vector P, N;
			double min;
			bool has_inter = s.intersection(r, P, N, id, min);

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
			cam[1] += 10;
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
			cam[1] += -10;
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
			ang[1] += 100;
			break;

		case 'k':
		case 'K':
			std::cout << "on baisse la tete " << ang[0] << " " << ang[1] << " " << ang[2] << std::endl;
			isMoving = true;
			ang[1] += -100;
			break;

		default:
			break;
	}
	if(isMoving){
		intersect2GL();
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
	int pixelsampling = 1;
	app.add_option("-p, --pixelsampling", pixelsampling, "Number of pixel sampling that you want.");
	std::string file = "defaut";
	app.add_option("-i,--input", file, "File name that you want to use. This file needs to be in json format.");
	std::string outputName;
	app.add_option("-o,--output", outputName, "File name of the output image")->required();
	CLI11_PARSE(app, argc, argv);
	
	read_scene_file(file.c_str());

	switch(number) {
		case 1:
			intersect1();
			break;
		case 2:
			createWindow(argc, argv);
			intersect2();
			intersect2GL();
			break;
		case 3:
			intersect3(pixelsampling);
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