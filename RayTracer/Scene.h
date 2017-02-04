#pragma once

#include <vector>
#include <algorithm>

#include <Camera.h>
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
Albega = reflect light / incident light
	   = % light reflected

Law of reflection: th_incident = th_reflection
	R = I - 2(N.I)N
Snell's Law: n1*sin(th1) = n2*sin(th2)
	n = refractive index = c/v
Transmission ray:
	T = eta*I + (eta*c1 -c2)N
	c1 = cos(th1) = I.N
	c2 = cos(th2) = sqrt( 1 - (eta^2)*sin(th1)^2 )
	sin(th1)^2 = 1 - cos(th1)^2 = 1 - c1*c1
	eta = n1/n2
Fresnel equations:
	R||	= ( (n2*cos(th1) - n1*cos(th2))/(n2*cos(th1) + n1*cos(th2)) )^2
	R_|_= ( (n1*cos(th1) - n2*cos(th2))/(n1*cos(th1) + n2*cos(th2)) )^2

	R = 1/2 * (R|| + R_|_)
	T = 1 - R 
*/

vec3 Reflect(const vec3& incident, const vec3& normal) {
	return incident - 2 * (dot(normal, incident))*normal;
}

vec3 Transmit(const vec3& incident, const vec3& normal, const float& n1, const float& n2) {
	float c1 = dot(normal, incident);

	// the normal vector should point *away* from the boundary we are approaching
	// there is a chance that we might approaching the boundary from the "inside"
	vec3 normalRef = normal;
	if (c1 < 0) { c1 = -c1; }
	else { normalRef = -normal; }

	float eta = n1 / n2;
	float k = 1 - eta*eta*(1 - c1*c1);
	if (k < 0) return vec3(0);	// total internal reflection (only if n2>n1)
	vec3 transmitted = eta*incident + (eta*c1 - sqrt(k))*normalRef;
	
	transmitted = normalize(transmitted);
	return transmitted;
}

float Fresnel(const vec3& incident, const vec3& transmission, const vec3& normal, const float& n1, const float& n2) {
	float cos1 = fabs(dot(incident, normal));
	float cos2 = fabs(dot(transmission, normal));

	float R_parallel = (n2*cos1 - n1*cos2) / (n2*cos1 + n1*cos2);
	float R_perpendidular = (n1*cos1 - n2*cos2) / (n1*cos1 + n2*cos2);

	float reflectionFraction = 0.5*(R_parallel*R_parallel + R_perpendidular*R_perpendidular);
	return reflectionFraction;
}

/*
PHONG MODEL:
ks, kd, ka are the sppecular, diffuse, and ambient constants,
alpha - shininess constant (for determining the size of the specular reflection)
L = object-to-light source direction vector
N = surface normal vector 
V = object-to-view/camera direction vector
R = object-to-reflection direction vector (R = 2(L.N)N -L)

Ip = k_a*i_a + sum(over all light sources m) k_d(L_m.N)i_m,d  +  k_s((R_m.V)^alpha)*(i_m,s)


Should I Tint Specular Reflections with the Objects' Color?
		NO		
if the material is delectric (i.e. non-conductor/non-metal)
A specular reflection is only a reflection of a light source off of the surface of an object.
*/



vec4 PhongShader(vec3 hitPosition, vec3 normal, Camera* camera, vector<LightSource*>* lightSources) {
	// constants
	const float ambientReflectionConstant = 0.1f;
	const float diffuseReflectionConstant = 0.8f;
	const float specularReflectionConstant = 0.5f;
	const int shininessConstant = 32;
	
	// ambient lighting
	vec3 ambientPortion = vec3(1, 1, 1); // doesn't belong to any particular light source
	vec3 diffusePortion = vec3(0, 0, 0);
	vec3 specularPortion = vec3(0, 0, 0);
	for (int i = 0; i<lightSources->size(); i++) {
		LightSource* light = (*lightSources)[i];

		// diffuse lighting
		vec3 hitToLightDir = normalize(light->position - hitPosition);
		float dotProduct = dot(normal, hitToLightDir);
		if (dotProduct < 0.0f) dotProduct = 0.0f;
		diffusePortion += dotProduct*light->colour;

		// specular lighting
		vec3 hitToCameraDir = normalize(camera->positionInWorldCoords - hitPosition);
		vec3 reflectDir = reflect(-hitToLightDir, normal);
		float powerTerm = pow(std::max(dot(hitToCameraDir, reflectDir), 0.0f), shininessConstant);
		specularPortion += powerTerm*light->colour; // case for tinting specular reflection with colour
	}

	vec3 totalLight = ambientReflectionConstant*ambientPortion + diffuseReflectionConstant*diffusePortion + specularReflectionConstant*specularPortion;

	//if (totalLight > 1.0f) totalLight = 1.0f;
	
	if (totalLight.x > 1.0f) totalLight.x = 1.0f; // this adds weird effects, like making the light too green and stuff
	if (totalLight.y > 1.0f) totalLight.y = 1.0f;
	if (totalLight.z > 1.0f) totalLight.z = 1.0f;
	
	return vec4(totalLight,1.0f);
}


Ray RayThruPixel(Camera* camera, int i, int j) {
	Ray ray;
	ray.origin = camera->positionInWorldCoords;

	float posX, posY;
	camera->pixelCoordsToCanvasPoint(i, j, posX, posY);
	ray.direction = (camera->screenDistance*-camera->backwardInWorldCoords) + (posX*camera->sideInWorldCoords) + (posY*camera->upwardInWorldCoords);

	ray.origin = camera->cameraCoordsToWorldCoordsMatrix * vec4(ray.origin,1);
	ray.direction = camera->cameraCoordsToWorldCoordsMatrix * vec4(ray.direction,1);
	// M*(P_cam - O_cam) = M*P_cam - M*O_cam = P_world - O_world

	return ray;
}

// Solve for t in the equation
//	ray.origin + ray.direction*t = targetPoint
float GetNumberOfUnitsAway(const Ray& ray, const vec3& targetPoint) {
	float answer;
	if (ray.direction.x != 0) {
		answer = (targetPoint.x - ray.origin.x) / ray.direction.x;
	}
	else if (ray.direction.y != 0) {
		answer = (targetPoint.y - ray.origin.y) / ray.direction.y;
	}
	else {
		answer = (targetPoint.z - ray.origin.z) / ray.direction.z;
	}
	if (answer < 0) answer = INFINITY;
	return answer;
}

Shape* GetClosestIntersection(Ray ray, Scene* scene, vec3& closestHit) {
	vec3 hit;
	float howManyUnitsAway, bestUnitsAway = INFINITY;
	bool foundHit = false;
	Shape* returnShape = nullptr;

	vector<Shape*> shapeList = *(scene->shapes);

	for (int i = 0; i < shapeList.size(); i++) {
		if (shapeList[i]->Intersect(ray, hit)) {
			howManyUnitsAway = GetNumberOfUnitsAway(ray, hit);
			if (howManyUnitsAway < 0 || howManyUnitsAway==INFINITY) continue; // ignore objects which are "behind" the ray, or which hit nothing
			// only save data if we've never hit before, or if the new collision is the closest yet in the forward direction
			if (!foundHit || howManyUnitsAway < bestUnitsAway) {
				foundHit = true;
				closestHit = hit;
				bestUnitsAway = howManyUnitsAway;
				returnShape = shapeList[i];				
			}
		}
	}
	return returnShape;
}

vec4 FindColour(Scene* scene, Camera* camera, Shape* shape, vec3 hitPosition) {
	vec3 normal = shape->GetNormal(hitPosition);
	
	vector<LightSource*>* lightList = (scene->lightSources);
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
			Shape* hitShape = GetClosestIntersection(ray, scene, closestHit);
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


vec4 TraceRayToColour(Ray incidentRay, Scene* scene, Camera* camera, int recursesLeft, float ambientRefractionIndex = 1.0f) {
	const float bias = 10e-3;
	vec4 colour;
	vec3 closestHit;
	vec3 reflectedRayDir, transmittedRayDir, reflectedRayPos, transmittedRayPos, biasVector;
	bool fromOutside;

	Shape* hitShape = GetClosestIntersection(incidentRay, scene, closestHit);
	if (hitShape == nullptr) return vec4(0,0,0,1); // we just return the background colour if we don't hit any objects
	vec3 normal = hitShape->GetNormal(closestHit);

	switch (hitShape->materialType) {
	case(REFLECTION_AND_REFRACTION): {
		if (recursesLeft == 0) return vec4(0,0,0,1);

		fromOutside = dot(normal, incidentRay.direction) < 0; //inefficient... I'm also doing this in the Transmit function
		float incomingRefractionIndex = fromOutside ? ambientRefractionIndex : hitShape->refractiveIndex;
		float outgoingRefractionIndex = fromOutside ? hitShape->refractiveIndex : ambientRefractionIndex;

		reflectedRayDir = Reflect(incidentRay.direction, normal);
		transmittedRayDir = Transmit(incidentRay.direction, normal, incomingRefractionIndex, outgoingRefractionIndex);
		float reflectedPercentage = Fresnel(incidentRay.direction, transmittedRayDir, normal, incomingRefractionIndex, outgoingRefractionIndex);

		biasVector = fromOutside ? bias*normal : -bias*normal; // points "up" back toward the incident ray direction
		reflectedRayPos = closestHit + biasVector;
		transmittedRayPos = closestHit - biasVector;

		vec4 reflectedColour = TraceRayToColour(Ray{ reflectedRayPos,reflectedRayDir }, scene, camera, recursesLeft - 1);
		vec4 transmittedColour = TraceRayToColour(Ray{ transmittedRayPos,transmittedRayDir }, scene, camera, recursesLeft - 1);
		return reflectedPercentage*reflectedColour + (1 - reflectedPercentage)*transmittedColour;
		}
		break;
	case(REFLECTION): {
		if (recursesLeft == 0) return vec4(0,0,0,1);
		fromOutside = dot(normal, incidentRay.direction) < 0; //inefficient... I'm also doing this in the Transmit function

		reflectedRayDir = Reflect(incidentRay.direction, normal);
		biasVector = fromOutside ? bias*normal : -bias*normal; // points "up" back toward the incident ray direction
		reflectedRayPos = closestHit + biasVector;

		return TraceRayToColour(Ray{ reflectedRayPos,reflectedRayDir }, scene, camera, recursesLeft - 1);
		}
		break;
	case(DIFFUSE):
		return FindColour(scene, camera, hitShape, closestHit);
		break;
	}
}

void RecursiveRayTrace(Camera* camera, Scene* scene, FrameBuffer* buffer, int recurseNumber) {
	for (int i = 0; i < buffer->height; i++) {
		for (int j = 0; j < buffer->width; j++) {

			//if (i == 400 && j == 500) {
			if (i == 500 && j == 250) {
				cout << "at the stop point" << endl;
			}

			Ray ray = RayThruPixel(camera, i, j);
			ray.direction = normalize(ray.direction);
			vec4 colour = TraceRayToColour(ray, scene, camera, recurseNumber);

			buffer->colourBuffer[(i*buffer->width + j) * 4 + 0] = colour.x;
			buffer->colourBuffer[(i*buffer->width + j) * 4 + 1] = colour.y;
			buffer->colourBuffer[(i*buffer->width + j) * 4 + 2] = colour.z;
			buffer->colourBuffer[(i*buffer->width + j) * 4 + 3] = colour.w;
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

					// if z doesn't fall between the near and far planes, ignore this pixel
					if (z < camera->screenDistance || camera->farplane < z) continue;

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