#pragma once

#include <vector>
#include <algorithm>

#include <Shapes.h>
#include <FrameBuffer.h>

using namespace std;

class Scene {
public:
	vector<Shape*>* shapes;
	vector<LightSource*>* lightSources;
	vec4 backgroundColour;

	Scene() : Scene(vec4(0,0,0,1)) { }
	Scene(vec4 bgColour) {
		this->shapes = new vector<Shape*>;
		this->lightSources = new vector<LightSource*>;
		this->backgroundColour = bgColour;
	}

	~Scene() {
		for (int i = 0; i < shapes->size(); i++) {
			delete (*shapes)[i];
		}
		for (int i = 0; i < lightSources->size(); i++) {
			delete (*lightSources)[i];
		}
		delete this->shapes;
		delete this->lightSources;
	}

};

/*
PHONG MODEL:
ks, kd, ka are the sppecular, diffuse, and ambient constants,
alpha - shininess constant (for determining the size of the specular reflection)
L = object-to-light source direction vector
N = surface normal vector 
V = object-to-view/camera direction vector
R = object-to-reflection direction vector (R = 2(L.N)N -L)

Ip = k_a*i_a + sum(over all light sources m) k_d(L_m.N)i_m,d  +  k_s((R_m.V)^alpha)*(i_m,s)
*/

vec4 PhongShader(vec3 hitPosition, vec3 normal, Camera* camera, vector<LightSource*>* lightSources) {
	// constants
	const float ambientReflectionConstant = 0.1f;
	const float diffuseReflectionConstant = 0.5f;
	const float specularReflectionConstant = 0.5f;
	const int shininessConstant = 32;
	
	// ambient lighting
	vec3 ambientPortion = ambientReflectionConstant*vec3(1, 1, 1); // doesn't belong to any particular light source
	vec3 diffusePortion = vec3(0, 0, 0);
	vec3 specularPortion = vec3(0, 0, 0);
	for (int i = 0; i<lightSources->size(); i++) {
		LightSource* light = (*lightSources)[i];

		// diffuse lighting
		vec3 hitToLightDir = normalize(light->position - hitPosition);
		float dotProduct = dot(normal, hitToLightDir);
		if (dotProduct < 0.0f) dotProduct = 0.0f;
		diffusePortion += diffuseReflectionConstant*dotProduct*light->colour;

		// specular lighting
		vec3 hitToCameraDir = normalize(camera->position - hitPosition);
		vec3 reflectDir = reflect(-hitToLightDir, normal);
		float powerTerm = pow(std::max(dot(hitToCameraDir, reflectDir), 0.0f), shininessConstant);
		specularPortion += specularReflectionConstant*powerTerm*light->colour;
	}

	vec3 totalLight = ambientPortion + diffusePortion + specularPortion;

	//if (totalLight > 1.0f) totalLight = 1.0f;
	
	if (totalLight.x > 1.0f) totalLight.x = 1.0f; // this adds weird effects, like making the light too green and stuff
	if (totalLight.y > 1.0f) totalLight.y = 1.0f;
	if (totalLight.z > 1.0f) totalLight.z = 1.0f;
	
	return vec4(totalLight,1.0f);
}


Ray RayThruPixel(Camera* camera, int i, int j) {
	Ray ray;
	ray.origin = camera->position;

	float posX = ((i + 0.5f) / camera->heightPixels)*camera->heightDistance - camera->heightDistance / 2.0f;
	float posY = (1 - ((j + 0.5f) / camera->widthPixels))*camera->widthDistance - camera->widthDistance / 2.0f;

	ray.direction = (camera->screenDistance*camera->forward) + (posX*camera->side) + (posY*camera->upward);
	return ray;
}

Shape* Intersect(Ray ray, Scene* scene, vec3& closestHit) {
	vec3 hit, normal;
	bool foundHit = false;
	Shape* returnShape = nullptr;

	vector<Shape*> shapeList = *(scene->shapes);

	for (int i = 0; i < shapeList.size(); i++) {
		if (shapeList[i]->Intersect(ray, hit, normal)) {
			// if we haven't found a collision, or if the new collision is closer
			if (!foundHit || hit.z < closestHit.z) {
				closestHit = hit;
				returnShape = shapeList[i];
			}
		}
	}
	return returnShape;
}

vec4 FindColour(Scene* scene, Camera* camera, Shape* shape, vec3 hitPosition) {
	vec3 normal = shape->GetNormal(hitPosition);
	
	vector<LightSource*>* lightList= (scene->lightSources);
	//LightSource* light = lightList[0];

	vec4 lightEffect = PhongShader(hitPosition, normal, camera, lightList);
	return lightEffect * shape->colour;
}

void RayTrace(Camera* camera, Scene* scene, FrameBuffer& buffer, int width, int height) {
	vec3 closestHit;
	vec4 colour;
	for (int i = 0; i<height; i++) {
		for (int j = 0; j<width; j++) {
			Ray ray = RayThruPixel(camera, i, j);
			Shape* hitShape = Intersect(ray, scene, closestHit);
			if (hitShape == nullptr) {
				colour = scene->backgroundColour;
			}
			else {
				colour = FindColour(scene, camera, hitShape, closestHit);
			}
			buffer.colourBuffer[(i*width + j)*4 + 0] = colour.x;
			buffer.colourBuffer[(i*width + j)*4 + 1] = colour.y;
			buffer.colourBuffer[(i*width + j)*4 + 2] = colour.z;
			buffer.colourBuffer[(i*width + j)*4 + 3] = colour.w;
		}
	}

}

/*
RayTrace(Camera cam, Scene scene, int width, int height) {
	Image image = new ImagE(width,height);
	for (int i=0; i<height; i++) {
		for(int j=0; j<width; j++) {
			Ray ray = RayThruPixel(cam, i, j);
			Intersection hit = Intersect(ray, scene);
			Image[i][j] = FindColor(hit);
		}
	}

}
*/