#ifndef SHAPES_H
#define SHAPES_H

#include<glm\glm.hpp>

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
	virtual bool Intersect(const Ray& ray, vec3& intersectionPoint, vec3& normal) = 0;
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
	bool Intersect(const Ray& ray, vec3& intersectionPoint, vec3& normal) {
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
			normal = this->GetNormal(intersectionPoint);
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

	bool Intersect(const Ray& ray, vec3& intersectionPoint, vec3& normal) {
		// based on (p - a).normal = 0
		float denom = dot(ray.direction, this->normalVector);
		if (denom == 0.0) return false;
		float t = dot(this->pA - ray.origin, this->normalVector) / denom;
		intersectionPoint = ray.origin + t*ray.direction;

		// ip = alpha*pA + beta*pB + gamma*pC; alpha + beta + gamma = 1;
		// ip-A = beta(B-A) + gamma(C-A)  <- three equations since it's a vector equation
		//   R  = beta*S    + gamma*T
		float rx = intersectionPoint.x - pA.x;
		float ry = intersectionPoint.y - pA.y;
		float sx = pB.x - pA.x;
		float sy = pB.y - pA.y;
		float tx = pC.x - pA.x;
		float ty = pC.y - pA.y;

		float denom2 = ty*sx - tx*sy;
		if (denom2 == 0.0 || sx == 0) return false;

		float gamma = (ry*sx - rx*sy) / denom2;
		if (!(0 <= gamma && gamma <= 1)) return false;
		float beta = (rx - gamma*tx) / sx;
		if (!(0 <= beta && beta <= 1)) return false;
		if (beta + gamma > 1) return false;

		normal = this->normalVector;
		return true;
	}

	vec3 GetNormal(const vec3& point) {
		return this->normalVector;
	}

};


#endif // !SHAPES_H

