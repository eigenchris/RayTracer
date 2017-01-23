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

	FrameBuffer frameBuffer = FrameBuffer(WIDTH, HEIGHT);
	Camera* camera = new Camera(WIDTH, HEIGHT);
	Scene* scene = new Scene();

	Sphere* sphere = new Sphere(vec3(-3, 0, -6), 1, vec4(0,0.5,1,1));
	Sphere* sphere2 = new Sphere(vec3(1, 2, -7), 1, vec4(0.6, 0.8, 0.4, 1));
	Sphere* sphere3 = new Sphere(vec3(2, -1, -5), 1, vec4(1.0, 0.35, 0.1, 1));
	//Triangle* triangle = new Triangle(vec3(2, 1, -5), vec3(1, -1, -6), vec3(3, -1, -3), vec4(1,0,0,1));
	scene->shapes->push_back((Shape*)sphere);
	scene->shapes->push_back((Shape*)sphere2);
	scene->shapes->push_back((Shape*)sphere3);
	//scene->shapes->push_back((Shape*)triangle);
	LightSource* light = new LightSource{ vec3(0, 2, -5), vec3(1, 1, 1) };
	LightSource* light2 = new LightSource{ vec3(0, -2, -5), vec3(1, 1, 1) };
	scene->lightSources->push_back(light);
	scene->lightSources->push_back(light2);


	RayTrace(camera, scene, frameBuffer, WIDTH, HEIGHT);

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