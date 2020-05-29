#include <stdio.h>
#include <iostream>
#include <vector>
#include "Vector.h"
#include <fstream>
#include <string>
#include <jsoncpp/json/json.h>


// Si mode 1 ou 2 lumiere statique et si 3 lumiere = premiere sphere 
void read_scene_file(const char* filename, Scene& s){
	std::ifstream file(filename);
	Json::Reader reader;
	Json::Value obj;
	reader.parse(file, obj);

	const Json::Value& spheres = obj["spheres"];
	const Json::Value& rectangles = obj["rectangles"];
	const Json::Value& triangles = obj["triangles"];
	const Json::Value& cylindres = obj["cylindres"];
	const Json::Value& lumiere = obj["intensite_lumiere"];



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
	
	s.intensite_lumiere = lumiere.asInt();

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

// const enmpeche de faire une copie locale
// On fait une fonction recursive pour calculer la reflexion des surfaces
// La fonction renvoie la couleur du pixel obtenu en envoyant un rayon R dans la scene S

int max(int a, int b){
	if (a > b){
		return a;
	}
	return b;
}

int main(int argc, char const *argv[]){
	
	Scene s = Scene();
	read_scene_file(argv[1], s);

	//Scene s;
	int H = 1024;
	int W = 1024;
	double fov = 90 * M_PI / 180; // angle de vue

	int nb_sampling = 1.0;

	std::vector<unsigned char> image(W * H *3);

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


	save_img("output.ppm", &image[0], W, H);

	std::cout << "hello" << std::endl;

	return 0;
}


