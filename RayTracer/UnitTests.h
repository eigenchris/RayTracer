#pragma once

#include <glm\glm.hpp>
#include <Shapes.h>
#include <Bezier.h>

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

bool BezierTest() {
	vec3 P[4];
	P[0] = vec3(0,0,0);
	P[1] = vec3(0,1,0);
	P[2] = vec3(1,1,0);
	P[3] = vec3(1,0,0);

	vec3 out1[4], out2[4];
	SplitBezierHalf(P, out1, out2);
	if( !(out1[1].x == 0 && out1[1].y == 0.5 && out1[2].x == 0.25 && out1[2].y == 0.75 && out1[3].x == 0.5 && out1[3].y == 0.75) ) {
		throw "Bezier Split doesn't work";
	}

	vec3 P2d[16];
	P2d[0] = vec3(0, 0, 0);
	P2d[1] = vec3(0, 1, 0);
	P2d[2] = vec3(0, 2, 0);
	P2d[3] = vec3(0, 3, 0);
	P2d[4] = vec3(1, 0, 0);
	P2d[5] = vec3(1, 1, 1);
	P2d[6] = vec3(1, 2, 1);
	P2d[7] = vec3(1, 3, 0);
	P2d[8] = vec3(2, 0, 0);
	P2d[9] = vec3(2, 1, 1);
	P2d[10] = vec3(2, 2, 1);
	P2d[11] = vec3(2, 3, 0);
	P2d[12] = vec3(3, 0, 0);
	P2d[13] = vec3(3, 1, 0);
	P2d[14] = vec3(3, 2, 0);
	P2d[15] = vec3(3, 3, 0);

	float t = 0.0;
	while (t <= 1.0) {
		vec3 v = Bezier3_2D(P2d, t, 0.5);
		cout << "t=" << t << "\t" << "B(" << t << ",0.5)= (" << v.x << "," << v.y << "," << v.z << ")" << endl;
		t += 0.01;
	}

	vec3 dU = dUBezier3_2D(P2d, 0.5, 0.5);
	vec3 dV = dUBezier3_2D(P2d, 0.5, 0.5);
	vec3 normal = Bezier3_2D_GetNormal(P2d, 0.5, 0.5);

	return true;
}