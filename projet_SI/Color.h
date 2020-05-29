#include <stdio.h>
#include <iostream>
#include <vector>
#include "reader.h"
#include <fstream>
#include <string>
class Color{
    private:
        double carre(double n){
            return n * n;
        }
        
        double PhongRef(const Vector& omi, const Vector& omo, const Vector& N, double phong_expo){
            Vector reflected_vec = omo.reflect(N);
            return (phong_expo + 2) / (2. * M_PI) * std::pow(dot(reflected_vec, omi), phong_expo);
        }

        void getColorMirror(const Ray& r, const Scene& s, int recursion, Vector& A, Vector& B, Vector* pixel_color){
            Vector mirror_dir = r.direction.reflect(A);
            Ray mirror_ray(B + 0.001 * A, mirror_dir);
            *pixel_color = getPixelColor(mirror_ray, s, recursion + 1);
        }

        void getColorGlass(const Ray& r, const Scene& s, int recursion, Vector& A, Vector& B, Vector* pixel_color){
            Vector glass_norm(A);
            double c1 = 1, c2 = 1.3;
            Ray new_ray;
            bool entering = true;
            if (0 < dot(r.direction, A)){ //On sort de l'objet
                c2 = 1;
                c1 = 1.3;
                glass_norm = -A;
                entering = false;
            }
            double rad = 1 - carre(c1 / c2) * (1 - carre(dot(glass_norm, r.direction)));
            if (0 < rad){
                Vector refracted_dir = (c1 / c2) * (r.direction - dot(r.direction, glass_norm) * glass_norm) - glass_norm * sqrt(rad);
                Ray refracted_ray(B - 0.001 * glass_norm, refracted_dir);
                // Utilisation du coef de shnik
                double R;
                double Rz = carre((c1 - c2) / (c1 + c2));
                if (entering)R = Rz + (1 - Rz) * std::pow(1 + dot(r.direction, A), 5);
                else R = Rz + (1 - Rz) * std::pow(1 - dot(refracted_dir, A), 5);
                if (generator_rand(generator) < R) new_ray = Ray(B + 0.001 * glass_norm, r.direction.reflect(A));
                else new_ray = Ray(B - 0.001 * glass_norm, refracted_dir);
                *pixel_color = getPixelColor(new_ray, s, recursion + 1);
            }
            // cas d'une reflection totale
            else new_ray = Ray(B + 0.001 * glass_norm, r.direction.reflect(A));
        }

        // Ajout éclairage direct
        void getDirectLightColor(const Ray& r, const Scene& s, int recursion, Vector& A, Vector& B, Vector* pixel_color, int id){
            // On envoi un rayon vers la sphere
            Vector axeOP = (B - s.light-> O).getNormalizedCopy();
            Vector rand_dir = random_cos(axeOP);
            Vector rand_point = rand_dir * s.light -> R + s.light -> O;
            Vector omi = (rand_point - B).getNormalizedCopy();
            Vector Np = rand_dir; // meme chos
            double d_light = (rand_point - B).getNorm();
            Ray rayLight(B + 0.01 * A , omi); // On envoi un rayon du point vers la lumiere
            // Au dessus on ajoute un truc à P pour décoller un peu l'origine du rayon et éviter le bruit sur l'image
            Vector A_light, B_light;
            int sphere_id_light;
            double t_light;
            bool intersection_light = s.intersection(rayLight, B_light, A_light, sphere_id_light, t_light); // On regarde si il y a une intersection entre le point qui va vers la lumiere et la lumiere
            // On relance un rayon vers la lumierre pour savoir s'il y a intersection, on ajoute 0.99 car la lumiere est une sphere
            if (intersection_light && ((t_light * t_light) < d_light * 0.99)) *pixel_color = Vector(0, 0, 0); 
            else *pixel_color = (s.light_intensity / (4 * M_PI * d_light) * std::max(0., dot(A, omi)) * dot(Np, -omi) / dot(axeOP, rand_dir)) * (M_PI) * ((1. - s.shapes[id] -> spec) * s.shapes[id] -> albedo / M_PI + PhongRef(omi, r.direction, A, s.shapes[id] -> phong_expo) * s.shapes[id] -> spec * s.shapes[id] -> albedo);
        }

        // Ajout éclairage indirect
        int getIndirectLightColor(const Ray& r, const Scene& s, int recursion, Vector& N, Vector& P, Vector* pixel_color, int id){
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
                if ((dot(random_direction, N) < 0) || (dot(random_direction, R) < 0)) return 1;
            }
            Ray random_ray(P + 0.001 * N, random_direction);
            double proba_phong = (s.shapes[id] -> phong_expo + 1) / (2. * M_PI) * std::pow(dot(R, random_direction), s.shapes[id] -> phong_expo);
            double proba_globale = p * dot(N, random_direction) / (2. * M_PI) + (1. - p) * proba_phong;	
            if (sample_diffuse) *pixel_color = *pixel_color + getPixelColor(random_ray, s, recursion + 1) * s.shapes[id] -> albedo * dot(N, random_direction) / M_PI / proba_globale;
            else *pixel_color = *pixel_color + getPixelColor(random_ray, s, recursion + 1) * dot(N, random_direction) * PhongRef(random_direction, r.direction, N, s.shapes[id] -> phong_expo) * s.shapes[id] -> spec * s.shapes[id] -> albedo / proba_globale;
            return 0;	
        }

    public:
        Color(){};
        // const enmpeche de faire une copie locale
        // On fait une fonction recursive pour calculer la reflexion des surfaces
        // La fonction renvoie la couleur du pixel obtenu en envoyant un rayon R dans la scene S
        Vector getPixelColor(const Ray& r, const Scene& s, int recursion = 0){
            if (recursion >= 5) return Vector(0, 0, 0);
            Vector pixel_color(0,0,0); // trois cmposantes du coup on prend un vecteur (A CHANGER)
            int id;
            double t;
            Vector P, N;
            bool has_inter = s.intersection(r, P, N, id, t);
            if (has_inter) {
                //On check  la nature de la forme que croise le rayon
                if (s.shapes[id] -> is_mirror) getColorMirror(r, s, recursion, N, P, &pixel_color);
                //calcule intensité pixel si transparant
                else if (s.shapes[id] -> is_transparent) getColorGlass(r, s, recursion, N, P, &pixel_color);
                //ni transparent ni miroir
                else {
                    getDirectLightColor(r, s, recursion, N, P, &pixel_color, id);
                    int bad_ray = getIndirectLightColor(r, s, recursion, N, P, &pixel_color, id);			
                    if (bad_ray) return Vector(0., 0., 0.);
                }
            }
            return pixel_color;
        }
};

