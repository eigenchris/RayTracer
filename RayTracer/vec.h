#ifndef VEC_H
#define VEC_H

#include <math.h>

class vec {

	float x,y,z,w;

	vec(float xx, float yy, float zz, float ww) : x(xx), y(yy), z(zz), w(ww) {}
	vec(float xx, float yy, float zz) : x(xx), y(yy), z(zz), w(0) {}
	vec() : x(0), y(0), z(0), w(0) {}

	vec operator +(vec u) {
		return vec(this->x + u.x, this->y + u.y, this->z + u.z, this->w + u.w);
	}
	vec operator -(vec u) {
		return vec(this->x - u.x, this->y - u.y, this->z - u.z, this->w - u.w);
	}
	vec operator -() {
		return vec(-this->x, -this->y, -this->z, -this->w);
	}
	vec operator *(vec u) {
		return vec(this->x * u.x, this->y * u.y, this->z * u.z, this->w * u.w);
	}
	vec operator *(float s) {
		return vec(this->x * s, this->y * s, this->z * s, this->w * s);
	}
	vec operator /(float s) {
		return vec(this->x / s, this->y / s, this->z / s, this->w / s);
	}

	float dot(vec u) {
		return (this->x)*(u.x) + (this->y)*(u.y) + (this->z)*(u.z) + (this->w)*(u.w);
	}
	float norm2() {
		return this->dot(*this);
	}
	float norm() {
		return sqrt(this->norm2());
	}
	void normalize() {
		*this = (*this) / this->norm();
	}

	vec cross(vec u) {
		return vec(
			this->y*u.z - this->z*u.y,
			this->z*u.x - this->x*u.z,
			this->x*u.y - this->y*u.x,
			0.0f
		);
	}

};


#endif // !VEC_H

