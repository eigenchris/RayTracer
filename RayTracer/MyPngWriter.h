#pragma once

#include <boost/gil/extension/io/png_io.hpp>

#include <FrameBuffer.h>

using namespace std;

void MyPngWriter(FrameBuffer* buffer, int width, int height, std::string filename) {
	int numPixels = height*width;
	float* colourBuffer = buffer->colourBuffer;

	unsigned char* r = new unsigned char[numPixels];  // red
	unsigned char* g = new unsigned char[numPixels];  // green
	unsigned char* b = new unsigned char[numPixels];  // blue
	unsigned char* a = new unsigned char[numPixels];  // alpha

	for (int i = 0; i < numPixels; i++) {
		r[i] = 255 * colourBuffer[4*i + 0];
		g[i] = 255 * colourBuffer[4*i + 1];
		b[i] = 255 * colourBuffer[4*i + 2];
		a[i] = 255 * colourBuffer[4*i + 3];
	}

	boost::gil::rgba8c_planar_view_t view = boost::gil::planar_rgba_view(width, height, r, g, b, a, width);
	boost::gil::png_write_view(filename, view);

	delete[] r;
	delete[] g;
	delete[] b;
	delete[] a;
}
