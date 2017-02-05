#pragma once

#include <glm\glm.hpp>
#include <Shapes.h>

float Random01() {
	return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

float Random11() {
	return 2*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) -1;
}

void RandomTripleSumTo1(float& alpha, float& beta, float& gamma) {
	float r1 = Random01();
	float r2 = (1 - r1)*Random01() + r1;
	alpha = r1;
	beta = r2 - r1;
	gamma = 1 - r2;
}

vec3 RandomV3() {
	return vec3(Random11(), Random11(), Random11());
}

vec3 RandomV2() {
	return vec3(Random11(), Random11(), 0);
}

Triangle RandomTriangle() {
	return Triangle(RandomV3(), RandomV3(), RandomV3());
}

Triangle RandomTriangle2D() {
	return Triangle(RandomV2(), RandomV2(), RandomV2());
}

bool TestTriangleInterestions() {
	int numTests = 50;

	for (int i = 0; i < numTests; i++) {
		cout << "inside triangle test: " << i << endl;
		Triangle tri = RandomTriangle();
		//tri = Triangle(vec3(1, 0, 0), vec3(0, 1, 0), vec3(1, 1, 0));
		float alpha, beta, gamma;
		RandomTripleSumTo1(alpha, beta, gamma);
		vec3 point = alpha*tri.pA + beta*tri.pB + gamma*tri.pC;

		float w0 = edgeFunctionForTriangle(tri.pB, tri.pC, point, tri.normalVector);
		float w1 = edgeFunctionForTriangle(tri.pC, tri.pA, point, tri.normalVector);
		float w2 = edgeFunctionForTriangle(tri.pA, tri.pB, point, tri.normalVector);

		bool allPlus = w0 >= 0 && w1 >= 0 && w2 >= 0;
		bool allMinus = w0 <= 0 && w1 <= 0 && w2 <= 0;

		if (!(allPlus || allMinus)) {
			throw "TestTriangleInterestions failed for points inside the triangle";
		}
	}

	for (int i = 0; i < numTests; i++) {
		cout << "outside triangle test: " << i << endl;
		Triangle tri = RandomTriangle();
		float alpha, beta, gamma;
		RandomTripleSumTo1(alpha, beta, gamma);
		alpha += 1;
		//vec3 point = alpha*tri.pA + beta*tri.pB + gamma*tri.pC; //bad... pA is a vector from the origin... not in the plane of the triangle... what if it's huge?

		vec3 u = tri.pB - tri.pA;
		vec3 v = cross(normalize(u), tri.normalVector);

		vec3 point = 1.1f*u + Random01()*v;

		if (i == 47) {
			cout << "pause" << endl;
		}

		float w0 = edgeFunctionForTriangle(tri.pB, tri.pC, point, tri.normalVector);
		float w1 = edgeFunctionForTriangle(tri.pC, tri.pA, point, tri.normalVector);
		float w2 = edgeFunctionForTriangle(tri.pA, tri.pB, point, tri.normalVector);

		bool allPlus = w0 >= 0 && w1 >= 0 && w2 >= 0;
		bool allMinus = w0 <= 0 && w1 <= 0 && w2 <= 0;



		if (allPlus || allMinus) {
			throw "TestTriangleInterestions failed for points outside the triangle";
		}
	}

	return true;
}

