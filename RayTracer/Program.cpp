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
#include <Camera.h>
#include <MyPngWriter.h>
#include <Bezier.h>

#include <UnitTests.h>

using namespace std;
using namespace glm;

void LoadScene1(Scene*);
void LoadScene2(Scene*);
void LoadScene3(Scene*);
void LoadScene4(Scene*);
void LoadScene5(Scene*);
void LoadScene6(Scene*);
void LoadScene7(Scene*);

/*
clock_t start = clock();
double duration = (clock() - start)/(double) CLOCKS_PER_SEC;
cout << "duration: " << duration << endl;
*/

int main() {
	cout << "Hello World!" << endl;

	const int WIDTH = 512;
	const int HEIGHT = 512;

	FrameBuffer* frameBuffer = new FrameBuffer(WIDTH, HEIGHT);
	Camera* camera = new Camera(WIDTH, HEIGHT);
	Scene* scene = new Scene();

	LoadScene7(scene);


	int recurseNumber = 5;
	//RecursiveRayTrace(camera, scene, frameBuffer, recurseNumber);

	Rasterizer(camera, scene, frameBuffer);

	MyPngWriter(frameBuffer, WIDTH, HEIGHT, "output.png");

	return 0;
}

// 3 spheres for ray caster
void LoadScene1(Scene* scene) {
	Sphere* sphere = new Sphere(vec3(-3, 0, -6), 1, vec4(0,0.5,1,1));
	Sphere* sphere2 = new Sphere(vec3(1, 2, -7), 1, vec4(0.6, 0.8, 0.4, 1));
	Sphere* sphere3 = new Sphere(vec3(2, -1, -5), 1, vec4(1.0, 0.35, 0.1, 1));
	scene->shapes->push_back((Shape*)sphere);
	scene->shapes->push_back((Shape*)sphere2);
	scene->shapes->push_back((Shape*)sphere3);
	LightSource* light = new LightSource{ vec3(0, 2, -5), vec3(1, 1, 1) };
	scene->lightSources->push_back(light);
}

// 3 triangles for ray caster
void LoadScene2(Scene* scene) {
	Triangle* triangle = new Triangle(vec3(-5, -2, -1), vec3(5, -2, -1), vec3(-5, 0, -10), vec4(0.65, 0, 0.9, 1));
	//Triangle* triangle = new Triangle(vec3(-5, -4, -3), vec3(5, -4, -3), vec3(-5, -3, -20), vec4(0.65,0,0.9,1));
	Triangle* triangle2 =new Triangle(vec3(5, -3, -20), vec3(-5, -3, -20), vec3(5, -4, -3), vec4(0.65,0,0.9,1));
	scene->shapes->push_back((Shape*)triangle);
	scene->shapes->push_back((Shape*)triangle2);
	LightSource* light = new LightSource{ vec3(0, 2, -5), vec3(1, 1, 1) };
	scene->lightSources->push_back(light);
}

// overlapping triangles for the rasterizer
void LoadScene3(Scene* scene) {
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
}

// glass ball, triangle, and floor
void LoadScene4(Scene* scene) {
	Triangle* tri1 = new Triangle(vec3(-4, -0.5, 0), vec3(4, -0.5, 0), vec3(-4, -0.5, -10), vec4(1, 0, 0, 1));
	Triangle* tri2 = new Triangle(vec3(-4, -0.5, -10), vec3(4, -0.5, 0), vec3(4, -0.5, -10), vec4(1, 0, 0, 1));
	
	Triangle* tri3 = new Triangle(vec3(1, 2, -4), vec3(0, 0.5, -4), vec3(2, 0.5, -4), vec4(0, 0, 1, 1.0));
	Triangle* tri4 = new Triangle(vec3(-2, -0.5, -10), vec3(2, -0.5, -10), vec3(0, 4, -10), vec4(0, 1, 0, 1.0));

	Sphere* sphere = new Sphere(vec3(0, 0, -4), 0.4, vec4(0, 0.5, 1, 1));
	sphere->refractiveIndex = 0.8;
	//sphere->materialType = REFLECTION_AND_REFRACTION;
	sphere->materialType = REFLECTION;

	//LightSource* light = new LightSource{ vec3(0.5, 2, -1), vec3(1, 1, 1) };
	LightSource* light = new LightSource{ vec3(0.5, 2, -3.8), vec3(1, 1, 1) };

	scene->shapes->push_back((Shape*)tri1);
	scene->shapes->push_back((Shape*)tri2);
	scene->shapes->push_back((Shape*)tri3);
	scene->shapes->push_back((Shape*)tri4);
	scene->shapes->push_back((Shape*)sphere);
	scene->lightSources->push_back(light);
}

// A sphere and lightsource at a variable distance from the screen
void LoadScene5(Scene* scene) {
	float distanceFromCamera = 2;
	
	vec3 d = vec3(0, 0, -distanceFromCamera);

	Sphere* sph = new Sphere(d, 1, vec4(0, 1, 0, 1));
	//scene->shapes->push_back((Shape*)sph);

	scene->shapes = TriangulateSphere(sph,4,4);

	LightSource* light = new LightSource{ vec3(1,1,1)+d, vec3(1, 1, 1) };
	//LightSource* light2 = new LightSource{ vec3(-2,1,1)+d, vec3(1, 1, 1) };
	scene->lightSources->push_back(light);
	//scene->lightSources->push_back(light2);
}


// Spheres to check coord system
void LoadScene6(Scene* scene) {
	Sphere* sph1 = new Sphere(vec3(2,0,-5), 1, vec4(1, 0, 0, 1));
	Sphere* sph2 = new Sphere(vec3(0,2,-5), 1, vec4(0, 1, 0, 1));
	Sphere* sph3 = new Sphere(vec3(0,0,-3), 1, vec4(0, 0, 1, 1));
	scene->shapes->push_back(sph1);
	scene->shapes->push_back(sph2);
	scene->shapes->push_back(sph3);
	
	LightSource* light = new LightSource{ vec3(1,1,1), vec3(1, 1, 1) };
	scene->lightSources->push_back(light);
}


// Simple Bezier surface
void LoadScene7(Scene* scene) {
	vec3 P2d[16];
	P2d[0] = vec3(0, 0, -5);
	P2d[1] = vec3(0, 1, -5);
	P2d[2] = vec3(0, 2, -5);
	P2d[3] = vec3(0, 3, -5);
	P2d[4] = vec3(1, 0, -5);
	P2d[5] = vec3(1, 1, -3);
	P2d[6] = vec3(1, 2, -3);
	P2d[7] = vec3(1, 3, -5);
	P2d[8] = vec3(2, 0, -5);
	P2d[9] = vec3(2, 1, -3);
	P2d[10] = vec3(2, 2, -3);
	P2d[11] = vec3(2, 3, -5);
	P2d[12] = vec3(3, 0, -5);
	P2d[13] = vec3(3, 1, -5);
	P2d[14] = vec3(3, 2, -5);
	P2d[15] = vec3(3, 3, -5);
	
	mat4 trans1 = mat4();
	trans1[3][2] = 5;
	mat4 rot = mat4();
	rot[0][0] = cos(45*M_PI/180); rot[0][2] = sin(45*M_PI/180);
	rot[2][0] = -sin(45*M_PI/180); rot[2][2] = cos(45*M_PI/180);
	mat4 trans2 = mat4();
	trans2[3][2] = -5;

	vec4 temp;
	for (int i = 0; i < 16; i++) {
		temp = vec4(P2d[i], 1.0);
		temp = trans1 * temp;
		temp = rot* temp;
		temp = trans2 * temp;
		P2d[i] = temp;
	}
		
	/*
	P2d[0] = vec3(0, 0, -3);
	P2d[1] = vec3(0, 0, -3);
	P2d[2] = vec3(0, 0, -3);
	P2d[3] = vec3(0, 0, -3);
	P2d[4] = vec3(1, 0, -4);
	P2d[5] = vec3(1, 1, -4);
	P2d[6] = vec3(1, 1, -4);
	P2d[7] = vec3(1, 0, -4);
	P2d[8] = vec3(2, 0, -5);
	P2d[9] = vec3(2, 1, -5);
	P2d[10] = vec3(2, 1, -5);
	P2d[11] = vec3(2, 0, -5);
	P2d[12] = vec3(3, 0, -6);
	P2d[13] = vec3(3, 0, -6);
	P2d[14] = vec3(3, 0, -6);
	P2d[15] = vec3(3, 0, -6);
	*/

	TriangulateBezierSurface(scene->shapes, P2d, 100, 100, vec4(0, 1, 0, 1));
	
	LightSource* light = new LightSource{ vec3(1,1,-1), vec3(1, 1, 1) };
	scene->lightSources->push_back(light);
}

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