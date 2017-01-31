// standard libs
#include <iostream>
#include <limits>

// packages
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// my own includes
#include <FrameBuffer.h>
#include <Shapes.h>
#include <Scene.h>
#include <MyPngWriter.h>

using namespace std;
using namespace glm;

int main() {
	cout << "Hello World!" << endl;

	const int WIDTH = 512;
	const int HEIGHT = 512;

	FrameBuffer* frameBuffer = new FrameBuffer(WIDTH, HEIGHT);
	Camera* camera = new Camera(WIDTH, HEIGHT);
	Scene* scene = new Scene();

//	Sphere* sphere = new Sphere(vec3(-3, 0, -6), 1, vec4(0,0.5,1,1));
//	Sphere* sphere2 = new Sphere(vec3(1, 2, -7), 1, vec4(0.6, 0.8, 0.4, 1));
//	Sphere* sphere3 = new Sphere(vec3(2, -1, -5), 1, vec4(1.0, 0.35, 0.1, 1));
//	Triangle* triangle = new Triangle(vec3(-5, -2, -1), vec3(5, -2, -1), vec3(-5, 0, -10), vec4(0.65, 0, 0.9, 1));
//	Triangle* triangle = new Triangle(vec3(-5, -4, -3), vec3(5, -4, -3), vec3(-5, -3, -20), vec4(0.65,0,0.9,1));
//	Triangle* triangle2 =new Triangle(vec3(5, -3, -20), vec3(-5, -3, -20), vec3(5, -4, -3), vec4(0.65,0,0.9,1));
//	scene->shapes->push_back((Shape*)sphere);
//	scene->shapes->push_back((Shape*)sphere2);
//	scene->shapes->push_back((Shape*)sphere3);
//	scene->shapes->push_back((Shape*)triangle);
//	scene->shapes->push_back((Shape*)triangle2);
//	LightSource* light = new LightSource{ vec3(0, 2, -5), vec3(1, 1, 1) };
//	LightSource* light2 = new LightSource{ vec3(0, -2, -5), vec3(1, 1, 1) };
//	scene->lightSources->push_back(light);
//	scene->lightSources->push_back(light2);

/*
	Triangle* tri1 = new Triangle(vec3(0, 0, -5), vec3(-2, 0, -3), vec3(-2, 2, -3), vec4(1.0, 0.3, 0.0, 1));
	Triangle* tri2 = new Triangle(vec3(0, 0, -5), vec3(2, 2, -3), vec3(0, 2, -3), vec4(0.0, 1.0, 0.3, 1));
	Triangle* tri3 = new Triangle(vec3(0, 0, -5), vec3(2, -2, -3), vec3(-2, -2, -3), vec4(0.3, 0.0, 1.0, 1));
	scene->shapes->push_back((Shape*)tri1);
	scene->shapes->push_back((Shape*)tri2);
	scene->shapes->push_back((Shape*)tri3);
	LightSource* light = new LightSource{ vec3(0, 0, -4), vec3(1, 1, 1) };
	scene->lightSources->push_back(light);
*/

	//vec3 r = camera->worldPointToPixelCoords(vec3(0,0,-1));
	//vec3 r2 = camera->worldPointToPixelCoords(vec3(-1, 1, -1));
	//vec3 r3 = camera->worldPointToPixelCoords(vec3(1, -1, -1));

	Triangle* tri1 = new Triangle(vec3(0, 1, -3), vec3(-1, -1, -2), vec3(1, -1, -5), vec4(1.0));
	tri1->cA = vec3(1, 0, 0);
	tri1->cB = vec3(0, 1, 0);
	tri1->cC = vec3(0, 0, 1);
	Triangle* tri2 = new Triangle(vec3(0.8, 1, -3.5), vec3(-0.2, -1, -3.5), vec3(1.8, -1, -3.5), vec4(1.0));
	tri2->cA = vec3(1, 0, 0);
	tri2->cB = vec3(0, 1, 0);
	tri2->cC = vec3(0, 0, 1);

	Triangle* tri3 = new Triangle(vec3(0.8, 1, 3.5), vec3(-0.2, -1, 3.5), vec3(1.8, -1, 3.5), vec4(1.0));
	tri3->cA = vec3(0, 0.5, 1);
	tri3->cB = vec3(0, 0.5, 1);
	tri3->cC = vec3(0, 0.5, 1);

	scene->shapes->push_back((Shape*)tri1);
	scene->shapes->push_back((Shape*)tri2);
	scene->shapes->push_back((Shape*)tri3);

	//RayTrace(camera, scene, frameBuffer);
	Rasterizer(camera, scene, frameBuffer);

	MyPngWriter(frameBuffer, WIDTH, HEIGHT, "output.png");

	return 0;
}

/*
FrameBuffer frameBuffer(3, 2);

Ray ray{ vec3(0,0,0), vec3(0,1,2) };
Sphere sph(vec3(0, 0, 2), 1);
Triangle tri(vec3(-1, 0, 1), vec3(1, 0, 1), vec3(0, 1, 1));

vec3 p;
vec3 n;
sph.Intersect(ray, p, n);
tri.Intersect(ray, p, n);
*/