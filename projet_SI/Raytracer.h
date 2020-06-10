/*!
 * \file Raytracer.h
 * \brief Define all reaytracers that can be used in this project
 * \author Pierre-Jean Besnard & Louis Billaut
 * \version 1.0
 */
#include <stdio.h>
#include <iostream>
#include <vector>
#include "Reader.h"
#include <fstream>
#include <string>

/*! \class ComplexRayTracer
 * \brief Allow to create a RayTracer wich use mirror and transparents objects, direct and indirect color and Phong exponent
 */
class ComplexRayTracer{
    private:
        /*!
        * \brief Calculate square of a number
        * 
        * \param n : the number wich will be squared
        * \return the squared number
        */
        double carre(double n){
            return n * n;
        }
        
        /*!
         * \brief Calculate the Phong reflection
         * 
         * \param omi : the first vector
         * \param omo : the second vector
         * \param N : the normal vector
         * \param phong_expo : the phong exponent 
         * \return the Phong reflection
         */
        double PhongRef(const Vector& omi, const Vector& omo, const Vector& N, double phong_expo){
            Vector reflected_vec = omo.reflect(N);
            return (phong_expo + 2) / (2. * M_PI) * std::pow(dot(reflected_vec, omi), phong_expo);
        }

        /*!
         * \brief Get the color of mirror objects
         * 
         * \param r : the ray light
         * \param s : the scene
         * \param recursion : the number of recursions
         * \param A : the A vector
         * \param B : the B vector
         * \param pixel_color : the pixel color
         */
        void getColorMirror(const Ray& r, const Scene& s, int recursion, Vector& A, Vector& B, Vector* pixel_color){
            Vector mirror_dir = r.direction.reflect(A);
            Ray mirror_ray(B + 0.001 * A, mirror_dir);
            *pixel_color = getPixelColor(mirror_ray, s, recursion + 1);
        }

        /*!
         * \brief Get the color of glass objects
         * 
         * \param r : the ray light
         * \param s : the scene
         * \param recursion : the number of recursions
         * \param A : the A vector
         * \param B : the B vector
         * \param pixel_color : the pixel color
         */
        void getColorGlass(const Ray& r, const Scene& s, int recursion, Vector& A, Vector& B, Vector* pixel_color){
            Vector glass_norm(A);
            double c1 = 1, c2 = 1.3;
            Ray new_ray;
            bool entering = true;
            if (0 < dot(r.direction, A)){ //Leaving object
                c2 = 1;
                c1 = 1.3;
                glass_norm = -A;
                entering = false;
            }
            double rad = 1 - carre(c1 / c2) * (1 - carre(dot(glass_norm, r.direction)));
            if (0 < rad){
                Vector refracted_dir = (c1 / c2) * (r.direction - dot(r.direction, glass_norm) * glass_norm) - glass_norm * sqrt(rad);
                Ray refracted_ray(B - 0.001 * glass_norm, refracted_dir);
               // Shnik coef utilisation
                double R;
                double Rz = carre((c1 - c2) / (c1 + c2));
                if (entering)R = Rz + (1 - Rz) * std::pow(1 + dot(r.direction, A), 5);
                else R = Rz + (1 - Rz) * std::pow(1 - dot(refracted_dir, A), 5);
                if (generator_rand(generator) < R) new_ray = Ray(B + 0.001 * glass_norm, r.direction.reflect(A));
                else new_ray = Ray(B - 0.001 * glass_norm, refracted_dir);
                *pixel_color = getPixelColor(new_ray, s, recursion + 1);
            }
            //total reflection case
            else new_ray = Ray(B + 0.001 * glass_norm, r.direction.reflect(A));
        }

        //Add direct lighting
        /*!
         * \brief get the direct light color of objects
         * 
         * \param r : the ray light
         * \param s : the scene
         * \param recursion : the number of recursions
         * \param A : the A vector
         * \param B : the B vector
         * \param pixel_color : the pixel color
         * \param id : id of the object
         */
        void getDirectLightColor(const Ray& r, const Scene& s, int recursion, Vector& A, Vector& B, Vector* pixel_color, int id){
            // On envoi un rayon vers la sphere
            Vector axeOP = (B - s.light-> O).getNormalizedCopy();
            Vector rand_dir = random_cos(axeOP);
            Vector rand_point = rand_dir * s.light -> R + s.light -> O;
            Vector omi = (rand_point - B).getNormalizedCopy();
            Vector Np = rand_dir; //same thing
            double d_light = (rand_point - B).getNorm();
            Ray rayLight(B + 0.01 * A , omi); //We send a ray from the point to the light
            //Above, allow to loosen the origin of the ray a little and avoid noise on the image
            Vector A_light, B_light;
            int sphere_id_light;
            double t_light;
            bool intersection_light = s.intersection(rayLight, B_light, A_light, sphere_id_light, t_light); // On regarde si il y a une intersection entre le point qui va vers la lumiere et la lumiere
            //We send again a ray towards the light to know if there is an intersection, we add 0.99 because the light is a sphere
            if (intersection_light && ((t_light * t_light) < d_light * 0.99)) *pixel_color = Vector(0, 0, 0); 
            else *pixel_color = (s.light_intensity / (4 * M_PI * d_light) * std::max(0., dot(A, omi)) * dot(Np, -omi) / dot(axeOP, rand_dir)) * (M_PI) * ((1. - s.shapes[id] -> spec) * s.shapes[id] -> color / M_PI + PhongRef(omi, r.direction, A, s.shapes[id] -> phong_expo) * s.shapes[id] -> spec * s.shapes[id] -> color);
        }

        //Add indirect lighting
        /*!
         * \brief get the indirect light color of objects
         * 
         * \param r : the ray light
         * \param s : the scene
         * \param recursion : the number of recursions
         * \param A : the A vector
         * \param B : the B vector
         * \param pixel_color : the pixel color
         * \param id : id of the object
         */
        int getIndirectLightColor(const Ray& r, const Scene& s, int recursion, Vector& N, Vector& P, Vector* pixel_color, int id){
            Vector random_direction;
            double p = 1 - s.shapes[id] -> spec;
            bool sample_diffuse;
            Vector R = r.direction.reflect(N);
            //We generate either a diffuse sample
            if (generator_rand(generator) < p){
                sample_diffuse = true;
                random_direction = random_cos(N);
            }
            //Either a specular sample
            else { 
                sample_diffuse = false;
                random_direction = random_phong(R, s.shapes[id] -> phong_expo);
                if ((dot(random_direction, N) < 0) || (dot(random_direction, R) < 0)) return 1;
            }
            Ray random_ray(P + 0.001 * N, random_direction);
            double proba_phong = (s.shapes[id] -> phong_expo + 1) / (2. * M_PI) * std::pow(dot(R, random_direction), s.shapes[id] -> phong_expo);
            double proba_globale = p * dot(N, random_direction) / (2. * M_PI) + (1. - p) * proba_phong;	
            if (sample_diffuse) *pixel_color = *pixel_color + getPixelColor(random_ray, s, recursion + 1) * s.shapes[id] -> color * dot(N, random_direction) / M_PI / proba_globale;
            else *pixel_color = *pixel_color + getPixelColor(random_ray, s, recursion + 1) * dot(N, random_direction) * PhongRef(random_direction, r.direction, N, s.shapes[id] -> phong_expo) * s.shapes[id] -> spec * s.shapes[id] -> color / proba_globale;
            return 0;	
        }

        //We do a recursive function to calculate the reflection of surfaces
        //The function returns the color of the pixel obtained by sending a ray R in the scene S
        /*!
         * \brief get the color of pixels of the image
         * 
         * \param r : the ray light
         * \param s : the scene
         * \param recursion : number of recursions
         */
        Vector getPixelColor(const Ray& r, const Scene& s, int recursion = 0){
            if (recursion >= 5) return Vector(0, 0, 0);
            Vector pixel_color(0,0,0); 
            int id;
            double t;
            Vector P, N;
            bool has_inter = s.intersection(r, P, N, id, t);
            if (has_inter) {
                //We check the nature of the shape crossed by the ray
                if (s.shapes[id] -> is_mirror) getColorMirror(r, s, recursion, N, P, &pixel_color);
                //calculate pixel intensity if transparent
                else if (s.shapes[id] -> is_transparent) getColorGlass(r, s, recursion, N, P, &pixel_color);
                //neither transparent neither mirror
                else {
                    getDirectLightColor(r, s, recursion, N, P, &pixel_color, id);
                    int bad_ray = getIndirectLightColor(r, s, recursion, N, P, &pixel_color, id);           
                    if (bad_ray) return Vector(0., 0., 0.);
                }
            }
            return pixel_color;
        }

        /*!
         * \brief define the color of pixels of the image
         * 
         * \param s : the scene
         * \param i : the i index
         * \param j : the j index
         * \param fov : the fov of the camera
         * \param nb_sampling : number of pixel sampling
         * \param image : the image
         * \param H : the height
         * \param W : the width
         */
        void setPixelColor(const Scene& s, int i, int j, double fov, double nb_sampling, unsigned char image[1024*1024*3], int H, int W){
            //Box Muller's method (anti-alliasing) -> PIXEL SAMPLING
            //created a random number which follows a Gaussian, it allows to send a ray not necessarily in the center
            double r1 = generator_rand(generator), r2 = generator_rand(generator);
            double dx = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2), dy = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2);
            Vector direction(j - W / 2 + 0.5 * dx , i - H / 2 + 0.5 * dy, - W / (2 * tan(fov / 2)));
            direction.normalize();
            Ray r(Vector(0, 0, 0), direction);
            Vector color(0., 0., 0.);
            ComplexRayTracer colorPixel;
            for (int i = 0; i < nb_sampling; i++){
                color = color + (colorPixel.getPixelColor(r, s) / nb_sampling);
            } 
            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0], 1 / 2.2))); // red
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1], 1 / 2.2))); // green
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2], 1 / 2.2))); // blue
            //The power allows to make a gamma correction, it makes the light colors lower because the human eye is more sensitive
        }

    public:
        /*!
         * \brief Constructor
         * Constructor of ComplexRayTracer class
         */
        ComplexRayTracer(){};
        // const prevents local copy
        // We do a recursive function to calculate the reflection of surfaces
        // The function returns the color of the pixel obtained by sending a ray R in the scene S

        /*!
         * \brief Create the image
         * 
         * \param s : the scene 
         * \param fov_cam : the fov of the camera
         * \param image : the image
         * \param H : the height
         * \param W : the width
         * \param ang : the angle of the camera
         * \param cam : the position of the camera 
         * \param nb_sampling : the number of pixel sampling
         */
        void createImage(const Scene& s, const int fov_cam, unsigned char image[1024 * 1024 * 3], int H, int W, double ang[2], double cam[3], int nb_sampling){
            double fov = fov_cam * M_PI / 180; // angle de vue
            for (int i = 0; i < H; i++) {
                for (int j = 0; j < W; j++) {
                    setPixelColor(s, i, j, fov, nb_sampling, image, H, W);
                }
            }
        }

};

/*! \class FlatPaintingRayTracer
 * \brief Create a Raytracer in flat painting
 */
class FlatPaintingRayTracer{
    private:
        /*!
         * \brief define the color of pixels of the image
         * 
         * \param s : the scene
         * \param i : the i index
         * \param j : the j index
         * \param fov : the fov of the camera
         * \param has_inter : boolean which represent if an intersection occured or not
         * \param nb_sampling : number of pixel sampling
         * \param image : the image
         * \param H : the height
         * \param W : the width
         * \param id : the id of the object
         */
        void setPixelColor(const Scene& s, int i, int j, double fov, bool has_inter, unsigned char image[1024*1024*3], int H, int W, int id){
            Vector pixel_intensity(0,0,0); // trois cmposantes du coup on prend un vecteur (A CHANGER)
            if (has_inter) pixel_intensity = s.shapes[id] -> color;
            image[((H - i - 1) * W + j) * 3 + 0] = pixel_intensity[0] * 255; // rouge 
            image[((H - i - 1) * W + j) * 3 + 1] = pixel_intensity[1] * 255; // vert
            image[((H - i - 1) * W + j) * 3 + 2] = pixel_intensity[2] * 255; // bleu
        }


    public:
        /*!
         * \brief Constructor
         * Constructor of FlatPaintingRayTracer class
         */
        FlatPaintingRayTracer(){};

        /*!
         * \brief Create the image
         * 
         * \param s : the scene 
         * \param fov_cam : the fov of the camera
         * \param image : the image
         * \param H : the height
         * \param W : the width
         * \param ang : the angle of the camera
         * \param cam : the position of the camera 
         */
        void createImage(const Scene& s, const int fov_cam, unsigned char image[1024*1024*3], int H, int W, double ang[2], double cam[3]){
            double fov = fov_cam * M_PI / 180; // angle de vue
            for (int i = 0; i < H; i++) {
                for (int j = 0; j < W; j++) {
                    Vector direction((j - W / 2) + ang[0], (i - H / 2) + ang[1], - W / (2 * tan(fov / 2)));
                    direction.normalize();
                    Ray r(Vector(cam[0], cam[1], cam[2]), direction);
                    int id;
                    Vector P, N;
                    double min;
                    bool has_inter = s.intersection(r, P, N, id, min);
                    setPixelColor(s, i, j, fov, has_inter, image, H, W, id);
                }
            }
        }

};


/*! \class SimpleRayTracer
 * \brief Create a simple Raytracer
 */
class SimpleRayTracer{
    private:
        /*!
         * \brief define the color of pixels of the image
         * 
         * \param s : the scene
         * \param i : the i index
         * \param j : the j index
         * \param fov : the fov of the camera
         * \param has_inter : boolean which represent if an intersection occured or not
         * \param nb_sampling : number of pixel sampling
         * \param image : the image
         * \param H : the height
         * \param W : the width
         * \param id : the id of the object
         * \param P : the intersection vector
         * \param N : the normal vector
         * \param light_position_value : the light position
         */
        void setPixelColor(const Scene& s, int i, int j, double fov, bool has_inter, unsigned char image[1024 * 1024 * 3], int H, int W, int id, Vector& P, Vector& N, Vector& light_position_value){
            int light_intensity = s.light_intensity / 1000;
            Vector pixel_intensity = 0; // trois cmposantes du coup on prend un vecteur (A CHANGER)
            if (has_inter) {
                pixel_intensity = s.shapes[id]->color * (light_intensity * std::max(0., dot((light_position_value - P).getNormalizedCopy(), N))) / (light_position_value - P).getNorm();
            }
            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., pixel_intensity[0])); // rouge 
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., pixel_intensity[1])); // vert
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., pixel_intensity[2])); // bleu
        }

        /*!
         * \brief define the color of pixels of the openGL image
         * 
         * \param s : the scene
         * \param i : the i index
         * \param j : the j index
         * \param fov : the fov of the camera
         * \param nb_sampling : number of pixel sampling
         * \param image : the image
         * \param H : the height
         * \param W : the width
         * \param id : the id of the object
         * \param P : the intersection vector
         * \param N : the normal vector
         * \param light_position_value : the light position
         */
        void setPixelColorGL(const Scene& s, int i, int j, double fov, bool has_inter, unsigned char imageGL[1024*1024*3], int H, int W, int id, Vector& P, Vector& N, Vector& light_position_value){
            Vector pixel_intensity = 0; // trois cmposantes du coup on prend un vecteur (A CHANGER)
            int light_intensity = s.light_intensity / 1000;
            if (has_inter) {
                pixel_intensity = s.shapes[id]->color * (light_intensity * std::max(0., dot((light_position_value - P).getNormalizedCopy(), N))) / (light_position_value - P).getNorm();
            }
            imageGL[((H - (H - i - 1) - 1) * W + j) * 3 + 0] = std::min(255., std::max(0., pixel_intensity[0])); // rouge 
            imageGL[((H - (H - i - 1) - 1) * W + j) * 3 + 1] = std::min(255., std::max(0., pixel_intensity[1])); // vert
            imageGL[((H - (H - i - 1) - 1) * W + j) * 3 + 2] = std::min(255., std::max(0., pixel_intensity[2])); // bleu
        }



    public:
        /*!
         * \brief Constructor
         * Constructor of SimpleRayTracer class
         * 
         * \param mode : the desired mode
         */
        SimpleRayTracer(int mode = 1) : mode(mode) {};
        int mode;

        /*!
         * \brief Create the image
         * 
         * \param s : the scene 
         * \param fov_cam : the fov of the camera
         * \param image : the image
         * \param H : the height
         * \param W : the width
         * \param ang : the angle of the camera
         * \param cam : the position of the camera 
         * \param light_position_value : the position of the light
         */
        void createImage(const Scene& s, const int fov_cam, unsigned char image[1024 * 1024 * 3], int H, int W, double ang[2], double cam[3], Vector& light_position_value){
            double fov = fov_cam * M_PI / 180; // angle de vue
            for (int i = 0; i < H; i++) {
                for (int j = 0; j < W; j++) {
                    Vector direction((j - W / 2) + ang[0], (i - H / 2) + ang[1], - W / (2 * tan(fov / 2)));
                    direction.normalize();
                    Ray r(Vector(cam[0], cam[1], cam[2]), direction);
                    int id;
                    Vector P, N;
                    double min;
                    bool has_inter = s.intersection(r, P, N, id, min);
                    if (mode == 1) setPixelColor(s, i, j, fov, has_inter, image, H, W, id, P, N, light_position_value);
                    else setPixelColorGL(s, i, j, fov, has_inter, image, H, W, id, P, N, light_position_value);
                }
            }
        }

};