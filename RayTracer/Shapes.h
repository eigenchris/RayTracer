#pragma once

#include <glm\glm.hpp>


using namespace glm;
using namespace std;



struct LightSource {
	vec3 position;
	vec3 colour;
};

struct Ray {
	vec3 origin;
	vec3 direction;
};

enum MaterialType {
	REFLECTION_AND_REFRACTION,
	REFLECTION,
	DIFFUSE
};

class Shape {
public:
	string tag = "Shape";
	vec4 colour = vec4(1, 1, 1, 1);
	MaterialType materialType;
	float refractiveIndex;
	virtual bool Intersect(const Ray& ray, vec3& intersectionPoint) = 0;
	virtual vec3 GetNormal(const vec3& point) = 0;

	Shape() {
		this->refractiveIndex = 1.0;
		this->materialType = DIFFUSE;
	}
};

class Sphere : public Shape {
public:
	vec3 center;
	float radius;
	

	Sphere(vec3 center, float radius) : Sphere(center, radius, vec4(1,1,1,1)) { }
	Sphere(vec3 center, float radius, vec4 colour) : center(center), radius(radius) {
		this->colour = colour;
		this->tag = "Sphere";
	}

	// based on (P-C).(P-C) = r^2 -> (dir.dir)t^2 + 2*dir.(org-C)t + |org-C|^2 - r^2 = 0 
	bool Intersect(const Ray& ray, vec3& intersectionPoint) {
		float a, b, c, discriminant;
		a = dot(ray.direction, ray.direction);
		b = 2*dot(ray.direction, ray.origin - this->center);
		float l = length(ray.origin - this->center);
		c = l*l - this->radius*this->radius;
		discriminant = b*b - 4 * a*c;

		if (discriminant >= 0.0) {
			float sq = sqrt(discriminant);
			float twiceA = 2 * a;
			float n1 = (-b + sq)/twiceA, n2 = (-b - sq)/twiceA;
			float t = (n1 < n2) ? n1 : n2; // set num to the point with the smaller (closer) z coordinate
			if (n1 < 0) t = n2;	// if the small z coordinate is negative, it is behind us
			if (n2 < 0) t = n1;	// if the small z coordinate is negative, it is behind us
			intersectionPoint = ray.origin + ray.direction*t;
			return true;
		}
		else {	
			return false;
		}
	}

	vec3 GetNormal(const vec3& point) {
		return normalize(point - center);
	}

 };

// Assume A, B, C are listed in counter-clockwise order
// U = B-A, and V = C-A,  mean that U x V should face in the same direction as the normal
//		(U x V).N > 0
// If A, B, C are in clockwise order, then (U x V).N < 0
float edgeFunctionForTriangle(const vec3& a, const vec3& b, const vec3& c, const vec3& normal)
{
	vec3 u = b - a;
	vec3 v = c - a;
	vec3 cr = cross(u, v);
	float dd = dot(cr, normal);

	return dot(cross(b - a, c - a), normal);
}

class Triangle : public Shape {
public:
	vec3 pA, pB, pC, normalVector;
	vec3 cA, cB, cC; // colour

	Triangle(vec3 a, vec3 b, vec3 c) : Triangle(a, b, c, vec4(1, 1, 1, 1)) { }
	Triangle(vec3 a, vec3 b, vec3 c, vec4 colour) : pA(a), pB(b), pC(c) {
		vec3 v1(b - a);
		vec3 v2(c - a);
		normalVector = normalize(cross(v1, v2));
		this->colour = colour;
		this->cA = colour;
		this->cB = colour;
		this->cC = colour;

		this->tag = "Triangle";
	}	

	bool Intersect(const Ray& ray, vec3& intersectionPoint) {
		// based on (p - a).normal = 0
		float denom = dot(ray.direction, this->normalVector);
		if (denom == 0.0) return false;
		float t = dot(this->pA - ray.origin, this->normalVector) / denom;
		intersectionPoint = ray.origin + t*ray.direction;

		// if the point is listed third, an interior point will always return a positive number
		float w0 = edgeFunctionForTriangle(pB, pC, intersectionPoint, this->normalVector);
		float w1 = edgeFunctionForTriangle(pC, pA, intersectionPoint, this->normalVector);
		float w2 = edgeFunctionForTriangle(pA, pB, intersectionPoint, this->normalVector);

		bool allPlus = w0 >= 0 && w1 >= 0 && w2 >= 0;
		bool allMinus = w0 <= 0 && w1 <= 0 && w2 <= 0;

		return allPlus || allMinus;


		/*
		// ip = alpha*pA + beta*pB + gamma*pC; alpha + beta + gamma = 1;
		// ip-A = beta(B-A) + gamma(C-A)  <- three equations since it's a vector equation
		//   R  = beta*S    + gamma*T
		vec3 R = intersectionPoint - pA;
		vec3 S = pB - pA;
		vec3 T = pC - pA;

		float gamma, denom2;
		bool sIsAlongZDirection = (S.x == 0.0 && S.y == 0.0);
		bool tIsAlongZDirection = (T.x == 0.0 && T.y == 0.0);
		if (sIsAlongZDirection || tIsAlongZDirection) { // if either S or T face completely in the Z direction...
			if (sIsAlongZDirection && !tIsAlongZDirection) { // assume T isn't along Z direction
				if (T.x != 0.0) { 
					gamma = R.x / T.x; 
				} else if (T.y != 0.0) { 
					gamma = R.y / T.y; 
				} else { 
					return false; 
				}
			}
			else if (!sIsAlongZDirection && tIsAlongZDirection) { // assume S isn't along Z direction
				if (S.x != 0.0) { 
					gamma = (R.z - (R.x*S.z / S.x))/T.z;
				} else if (S.y != 0.0) { 
					gamma = (R.z - (R.y*S.z / S.y)) / T.z;
				} else { 
					return false; 
				}
			}
			else { // this isn't a triangle... 
				return false;
			}
		} else { // otherwise, proceed as normal
			float denom2 = T.y*S.x - T.x*S.y;
			if (denom2 == 0.0) return false;
			gamma = (R.y*S.x - R.x*S.y) / denom2;
		}
		if (!(0 <= gamma && gamma <= 1)) return false;
	
		float beta;
		if (S.x != 0.0) { beta = (R.x - gamma*T.x) / S.x; }
		else if (S.y != 0.0) { beta = (R.y - gamma*T.y) / S.y; }
		else if (S.z != 0.0) { beta = (R.z - gamma*T.z) / S.z; }
		else { return false; }

		if (!(0 <= beta && beta <= 1)) return false;
		if (beta + gamma > 1) return false;

		normal = this->normalVector;
		return true;
		*/
	}

	vec3 GetNormal(const vec3& point) {
		return this->normalVector;
	}

};
