/*!
 * \file Reader.h
 * \brief Define JSON readers
 * \author Pierre-Jean Besnard & Louis Billaut
 * \version 1.0
 */
#include <stdio.h>
#include <iostream>
#include <vector>
#include "Vector.h"
#include <fstream>
#include <string>
#include <jsoncpp/json/json.h>

/*! \class Reader
 * \brief Allows to create a JSON file reader
 */
class Reader {
	public:
		/*!
		 * \brief Read a JSON file and create a scene with objects of this file
		 * 
		 * \param filename : the file name
		 * \param scene : the scene where objects will be placed
		 * \param l : the light position
		 * \param f : the fov of the camera
		 */
		void read_scene_file(const char* filename, Scene& scene, Vector& l, int &f){
			std::ifstream file(filename);
			Json::Reader reader;
			Json::Value obj;
			reader.parse(file, obj);
			get_spheres_from_file(obj);
			get_cylinders_from_file(obj);
			get_triangles_from_file(obj);
			get_rectangles_from_file(obj);
			get_light_and_cam_from_file(obj);
			scene = s;
			l = light_position_value;
			f = fov_cam;
		}
		Scene s; /*!< the scene */
		Vector light_position_value; /*!< the light position */
		int fov_cam; /*!< the fov of the camera */

	private: 
		/*!
		 * \brief Get rectangles from the JSON file
		 * 
		 * \param obj : a JSON object
		 */
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

		/*!
		 * \brief Get triangles from the JSON file
		 * 
		 * \param obj : a JSON object
		 */
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

		/*!
		 * \brief Get cylinders from the JSON file
		 * 
		 * \param obj : a JSON object
		 */
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

		/*!
		 * \brief Get spheres from the JSON file
		 * 
		 * \param obj : a JSON object
		 */
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

		/*!
		 * \brief Get the light and the fov of the camera from the JSON file
		 * 
		 * \param obj : a JSON object
		 */
		void get_light_and_cam_from_file(Json::Value obj){
			const Json::Value& light_intensity_value = obj["light_intensity"];
			const Json::Value& light_position = obj["light_position"];
			const Json::Value& fov = obj["fov"];
			s.light_intensity = light_intensity_value.asDouble();
			light_position_value = Vector(light_position[0].asInt(), light_position[1].asInt(), light_position[2].asInt());
			fov_cam = fov.asInt();
		}
};

