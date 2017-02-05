#pragma once

#include <stdint.h>
#include <string.h>
#include <limits>

class FrameBuffer {
public:

	uint32_t width, height;
	float* colourBuffer;
	float* zBuffer;

	FrameBuffer(int w, int h) {
		this->width = w;
		this->height = h;

		this->colourBuffer = new float[w*h*4];
		this->zBuffer = new float[w*h];
		
		this->ClearBuffers();
	}

	~FrameBuffer() {
		delete[] this->colourBuffer;
		delete[] this->zBuffer;
	}

	void ClearColourBuffer() {
		// beware, memset works with char, not int
		memset(this->colourBuffer, 0, sizeof(float)*width*height * 4);
	}
	void ClearZBuffer() {
		std::fill(zBuffer, zBuffer + (width*height), std::numeric_limits<float>::infinity());
	}
	void ClearBuffers() {
		this->ClearColourBuffer();
		this->ClearZBuffer();
	}

};
