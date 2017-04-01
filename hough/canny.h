#ifndef _CANNY_
#define _CANNY_

#include <string>
#include <cmath>
#include "CImg.h"

using namespace cimg_library;
using namespace std;

class Canny {
	
public:

    CImg<unsigned char> img;
    CImg<unsigned char> grayscaled;
    CImg<unsigned char> edge;

	unsigned char *data; /* input image */
    int width;           
    int height;
    int *idata;          /* output for edges */
    int *magnitude;      /* edge magnitude as detected by Gaussians */
    float *xConv;        /* temporary for convolution in x direction */
    float *yConv;        /* temporary for convolution in y direction */
    float *xGradient;    /* gradients in x direction, as detected by Gaussians */
    float *yGradient;    /* gradients in y direction,a s detected by Gaussians */

	Canny(const char *filename);
	~Canny();
	CImg<unsigned char> toGrayScale(int width, int height);
	void performcanny();
    void performcannyparam(float lowThreshold, float highthreshold, 
						  float gaussiankernelradius, int gaussiankernelwidth,  
						  int contrastnormalised);
	void normalizeContrast();
	int computeGradients(float kernelRadius, int kernelWidth);
	float hypotenuse(float x, float y);
	float gaussian(float x, float sigma);
	void performHysteresis(int low, int high);
	void follow(int x1, int y1, int i1, int threshold);
	
};

#endif
