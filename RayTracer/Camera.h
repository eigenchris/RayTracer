#pragma once

#include <glm\glm.hpp>

class Camera {
private:
	Camera(vec3 pos, vec3 fwd, vec3 up, int width, int height) : widthPixels(width), heightPixels(height) {
		SetVectors(pos, fwd, up);

		widthDistance = 2;
		heightDistance = 2;
		screenDistance = 1;
		farplane = 100;
	}

public:
	vec3 positionInWorldCoords, backwardInWorldCoords, upwardInWorldCoords, sideInWorldCoords;
	const int widthPixels, heightPixels;
	float widthDistance, heightDistance, screenDistance, farplane;
	glm::mat4 cameraCoordsToWorldCoordsMatrix;
	glm::mat4 worldCoordsToCameraCoordsMatrix;

	Camera(int width, int height) : Camera(vec3(0, 0, 0), vec3(0, 0, 1), vec3(0, 1, 0), width, height) {}

	void SetVectors(vec3 pos, vec3 bkwd, vec3 up) {
		positionInWorldCoords = pos;
		backwardInWorldCoords = normalize(bkwd);
		upwardInWorldCoords = normalize(up);
		sideInWorldCoords = cross(up, backwardInWorldCoords);

		this->cameraCoordsToWorldCoordsMatrix = glm::mat4(vec4(sideInWorldCoords, 0), vec4(upwardInWorldCoords, 0), vec4(backwardInWorldCoords, 0), vec4(positionInWorldCoords, 1));
		this->worldCoordsToCameraCoordsMatrix = glm::inverse(cameraCoordsToWorldCoordsMatrix);
	}

	// third component returns z-depth value
	vec3 worldPointToProjectedCanvasPoint(vec3 pointInWorldCoords) {
		vec4 pointInCameraCoords = this->worldCoordsToCameraCoordsMatrix * vec4(pointInWorldCoords, 1.0f);
		vec3 projection;
		projection.x = this->screenDistance * pointInCameraCoords.x / -pointInCameraCoords.z;
		projection.y = this->screenDistance * pointInCameraCoords.y / -pointInCameraCoords.z;
		projection.z = -pointInCameraCoords.z; //+z poitns outward toward us... it's better to count z distance into the screen as positive.
		return projection;
	}

	// third component returns z-depth value
	vec3 worldPointToPixelCoords(vec3 pointInWorldCoords) {
		vec3 canvasCoords = this->worldPointToProjectedCanvasPoint(pointInWorldCoords);
		vec2 pixelCoords = this->canvasPointToPixelCoordsV(vec2(canvasCoords.x, canvasCoords.y));
		return vec3(pixelCoords, canvasCoords.z);
	}

	vec2 canvasPointToPixelCoordsV(const vec2& canvasPoint) {
		//float i = (1 - ((canvasPoint.y + this->heightDistance / 2.0f) / (this->heightDistance))) * this->heightPixels;
		//float j = ((canvasPoint.x + this->widthDistance / 2.0f) / (this->widthDistance)) * this->widthPixels;
		float i = (1.0f / 2.0f - (canvasPoint.y / this->heightDistance)) * this->heightPixels;
		float j = ((canvasPoint.x / this->widthDistance) + 1.0f / 2.0f) * this->widthPixels;
		return vec2(j, i);
	}

	void canvasPointToPixelCoords(float canvasX, float canvasY, int& i, int& j) {
		i = (1 - ((canvasY + this->heightDistance / 2.0f) / (this->heightDistance))) * this->heightPixels;
		j = ((canvasX + this->widthDistance / 2.0f) / (this->widthDistance)) * this->widthPixels;
	}

	void pixelCoordsToCanvasPoint(int i, int j, float& canvasX, float& canvasY) {
		canvasX = ((j + 0.5f) / this->heightPixels)*this->heightDistance - this->heightDistance / 2.0f;
		canvasY = (1 - ((i + 0.5f) / this->widthPixels))*this->widthDistance - this->widthDistance / 2.0f;
	}


};