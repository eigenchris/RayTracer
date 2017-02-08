#pragma once

#include <glm/glm.hpp>

using namespace glm;

// B(t) = P0*(1-t)^3 + 3*P1*t*(1-t)^2 + 3*P2*(t^2)*(1-t) + P3*t^3
// B(t) = P0*(1 - 3t + 3t^2 -t^3) P1*(0 + 3t -6t^2 + 3t^3) + P2*(0 + 0t + 3t^2 - 3t^3) + P3*(0 + 0t + 0t^2 + t^3)
// B(t) =   [P0 P1 P2 P3] * M * [1 t t^2 t^3]^T
// P0 is top row, 1 is left column

// B'(t) = P0(3 + 6t -3t^2) + P1(3 - 12t + 9t^2) + P2(0 + 6t - 9t^2) + P3(0 + 0t + 3t^2)

vec3 Bezier3_1D(vec3 P[4], float t) {
	const mat4 BM = transpose(mat4(vec4(1, -3, 3, -1), vec4(0, 3, -6, 3), vec4(0, 0, 3, -3), vec4(0, 0, 0, 1)));

	vec4 T = vec4(1, t, t*t, t*t*t);
	T = BM*T;
	vec3 answer = P[0]*T.x + P[1]*T.y + P[2]*T.z + P[3]*T.w;
	return answer;
}

vec3 dBezier3_1D(vec3 P[4], float t) {
	const mat4 BM = transpose(mat4(vec4(3, 6, -3, 0), vec4(3, -12, 9, 0), vec4(0, 6, -9, 0), vec4(0, 0, 3, 0)));

	vec4 T = vec4(1, t, t*t, 0);
	T = BM*T;
	vec3 answer = P[0]*T.x + P[1]*T.y + P[2]*T.z + P[3]*T.w;
	return answer;
}

//P(u,v) = sum_i sum_j Bi(u) * Bj(v) * P_ij
vec3 Bezier3_2D(vec3* P, float u, float v)
{
	vec3 Pu[4];
	// compute 4 control points along u direction
	for (int i = 0; i < 4; ++i) {
		vec3 curveP[4];
		curveP[0] = P[i * 4];
		curveP[1] = P[i * 4 + 1];
		curveP[2] = P[i * 4 + 2];
		curveP[3] = P[i * 4 + 3];
		Pu[i] = Bezier3_1D(curveP, u);
	}
	// compute final position on the surface using v
	return Bezier3_1D(Pu, v);
}

vec3 dUBezier3_2D(vec3 P[16], float u, float v)
{
	vec3 dPu[4];
	// compute 4 control points along u direction
	for (int i = 0; i < 4; ++i) {
		vec3 curveP[4];
		curveP[0] = P[i * 4];
		curveP[1] = P[i * 4 + 1];
		curveP[2] = P[i * 4 + 2];
		curveP[3] = P[i * 4 + 3];
		dPu[i] = dBezier3_1D(curveP, u);
	}
	// compute final position on the surface using v
	return Bezier3_1D(dPu, v);
}

vec3 dVBezier3_2D(vec3 P[16], float u, float v)
{
	vec3 Pu[4];
	// compute 4 control points along u direction
	for (int i = 0; i < 4; ++i) {
		vec3 curveP[4];
		curveP[0] = P[i * 4];
		curveP[1] = P[i * 4 + 1];
		curveP[2] = P[i * 4 + 2];
		curveP[3] = P[i * 4 + 3];
		Pu[i] = Bezier3_1D(curveP, u);
	}
	// compute final position on the surface using v
	return dBezier3_1D(Pu, v);
}

vec3 Bezier3_2D_GetNormal(vec3 P[16], float u, float v) {
	vec3 dU = dUBezier3_2D(P, u, v);
	vec3 dV = dVBezier3_2D(P, u, v);
	vec3 normalVector = normalize(cross(dU, dV));
	return normalVector;
}

void SplitBezier_1D(vec3 in[4], float t, vec3 out1[4], vec3 out2[4]) {
	float oneMinusT = 1 - t;

	vec3 p01 = oneMinusT*in[0] + t*in[1];
	vec3 p12 = oneMinusT*in[1] + t*in[2];
	vec3 p23 = oneMinusT*in[2] + t*in[3];

	vec3 p01_12 = oneMinusT*p01 + t*p12;
	vec3 p12_23 = oneMinusT*p12 + t*p23;
	vec3 p0123 = oneMinusT*p01_12 + t*p12_23;

	out1[0] = in[0];	out1[1] = p01;		out1[2] = p01_12;	out1[3] = p0123;
	out2[0] = p0123;	out2[1] = p12_23;	out2[2] = p23;		out2[3] = in[3];		
}


void SplitBezierHalf(vec3 P[4], vec3 out1[4], vec3 out2[4]) {
	out1[0] = P[0];
	out1[1] = (P[0] + P[1]) * 0.5f;
	out1[2] = (P[0] + 2.0f*P[1] + P[2]) * 0.25f;
	out1[3] = out2[0] = (P[0] + 3.0f*(P[1] + P[2]) + P[3]) * 0.125f;
	out2[1] = (P[3] + 2.0f*P[2] + P[1]) * 0.25f;
	out2[2] = (P[2] + P[3]) * 0.5f;
	out2[3] = P[3];
}
