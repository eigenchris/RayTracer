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
	const float diffuseReflectionConstant = 0.8f;
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

	ray.direction = (camera->screenDistance*camera->backward) + (posX*camera->side) + (posY*camera->upward);
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

void RayTrace(Camera* camera, Scene* scene, FrameBuffer* buffer) {
	vec3 closestHit;
	vec4 colour;
	for (int i = 0; i<buffer->height; i++) {
		for (int j = 0; j<buffer->width; j++) {
			Ray ray = RayThruPixel(camera, i, j);
			Shape* hitShape = Intersect(ray, scene, closestHit);
			if (hitShape == nullptr) {
				colour = scene->backgroundColour;
			}
			else {
				colour = FindColour(scene, camera, hitShape, closestHit);
			}
			buffer->colourBuffer[(i*buffer->width + j)*4 + 0] = colour.x;
			buffer->colourBuffer[(i*buffer->width + j)*4 + 1] = colour.y;
			buffer->colourBuffer[(i*buffer->width + j)*4 + 2] = colour.z;
			buffer->colourBuffer[(i*buffer->width + j)*4 + 3] = colour.w;
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

// U = C-A,	 V = B-A
// This formula calculates U x V; the z-coords of U and V are zero, so the result only points in the Z-direction.
// Taking, V to point "up", the result is + if  U is right of V, and - if U is left of V
//		i.e. + is c is right of V, - is c left of V
// The result is the area of the parallelogram formed by U and V
//		if we divide by 2, we get the area of the triangle formed by (a,b,c)
float edgeFunction(const vec3& a, const vec3& b, const vec3& c)
{
	return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
}

// Perspective-correct colouring and z-values are found from the following page:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/perspective-correct-interpolation-vertex-attributes
void Rasterizer(Camera* camera, Scene* scene, FrameBuffer* frameBuffer) {
	int imageWidth = frameBuffer->width;
	vector<Shape*> shapesList = *scene->shapes; //TODO: does this move the list onto the stack?

	// Set all colour values to zero; Set depth buffer all all INFINITY
	frameBuffer->ClearBuffers();

	for (int idx = 0; idx < shapesList.size(); idx++) {
		if (shapesList[idx]->tag != "Triangle") continue;
		Triangle* triangle = (Triangle*)shapesList[idx];
		vec3 v0 = camera->worldPointToPixelCoords(triangle->pA);
		vec3 v1 = camera->worldPointToPixelCoords(triangle->pB);
		vec3 v2 = camera->worldPointToPixelCoords(triangle->pC);
		// prepare C/Z ratios for perspective colouring
		vec3 c0 = triangle->cA / v0.z;
		vec3 c1 = triangle->cB / v1.z;
		vec3 c2 = triangle->cC / v2.z;

		// get bounds of rasterization box, and clip them if needed
		vec2 minV = vec2(std::min(std::min(v0.x, v1.x), v2.x), std::min(std::min(v0.y, v1.y), v2.y));
		vec2 maxV = vec2(std::max(std::max(v0.x, v1.x), v2.x), std::max(std::max(v0.y, v1.y), v2.y));
		if (minV.x < 0) minV.x = 0;
		if (minV.y < 0) minV.y = 0;
		if (maxV.x >= camera->heightPixels)	minV.x = camera->heightPixels -1;
		if (maxV.y >= camera->widthPixels) minV.y = camera->widthPixels -1;

		float screenArea = edgeFunction(v0, v1, v2);

		// loop over all pixels in bounding box
		for (int i = floor(minV.y); i <= ceil(maxV.y); i++) {
			for (int j = floor(minV.x); j <= ceil(maxV.x); j++) {
				vec3 p = vec3(j+0.5, i+0.5, 0.0f); //3rd coord is arbitrary

				// check if pixel is inside triangle
				float w0 = edgeFunction(v1, v2, p);
				float w1 = edgeFunction(v2, v0, p);
				float w2 = edgeFunction(v0, v1, p);

				if (w0 >= 0 && w1 >= 0 && w2 >= 0) {					
					// get barycentric coordinates
					w0 /= screenArea;
					w1 /= screenArea;
					w2 /= screenArea;
					
					// get true Z of point (perspective)
					//	1/Z = w0/Z0 + w1/Z1 + w2/Z2
					float z = 1/(w0/v0.z + w1/v1.z + w2/v2.z);

					// if z is father away than the current object, then ignore this pixel
					if (z >= frameBuffer->zBuffer[i*imageWidth + j]) continue;			

					// Update depth-buffer with that depth
					frameBuffer->zBuffer[i*imageWidth + j] = z;

					// get true colour of point (perspective)
					// C/Z = w0*(C0/Z0) + w1*(C1/Z1) + w2*(C2/Z2), with C being any given colour component 
					float r = z*(w0*c0.x + w1*c1.x + w2*c2.x);
					float g = z*(w0*c0.y + w1*c1.y + w2*c2.y);
					float b = z*(w0*c0.z + w1*c1.z + w2*c2.z);

					frameBuffer->colourBuffer[(i*imageWidth + j) * 4 + 0] = r;
					frameBuffer->colourBuffer[(i*imageWidth + j) * 4 + 1] = g;
					frameBuffer->colourBuffer[(i*imageWidth + j) * 4 + 2] = b;
					frameBuffer->colourBuffer[(i*imageWidth + j) * 4 + 3] = 1.0f;
				}

			}
		}
	}
}




/*
// rasterization algorithm
for (each triangle in scene) {
// STEP 1: project vertices of the triangle using perspective projection
Vec2f v0 = perspectiveProject(triangle[i].v0);
Vec2f v1 = perspectiveProject(triangle[i].v1);
Vec2f v2 = perspectiveProject(triangle[i].v2);
for (each pixel in image) {
// STEP 2: is this pixel contained in the projected image of the triangle?
if (pixelContainedIn2DTriangle(v0, v1, v2, x, y)) {
image(x,y) = triangle[i].color;
}
}
}


*/