#include "canny.h"

#define ffabs(x) ( (x) >= 0 ? (x) : -(x) ) 
#define GAUSSIAN_CUT_OFF 0.005f
#define MAGNITUDE_SCALE 100.0f
#define MAGNITUDE_LIMIT 1000.0f
#define MAGNITUDE_MAX ((int) (MAGNITUDE_SCALE * MAGNITUDE_LIMIT))

Canny::Canny(const char *filename) {

    img.load_bmp(filename);

    width = img.width();
    height = img.height();

    grayscaled = toGrayScale(width, height);
    edge = CImg<unsigned char>(width, height, 1, 1);

    int size = width * height;

    data = new unsigned char[size];
    idata = new int[size];
    magnitude = new int[size];
    xConv = new float[size];
    yConv = new float[size];
    xGradient = new float[size];
    yGradient = new float[size];

    cimg_forXY(grayscaled, x, y) {
    	data[x + y * width] = grayscaled(x, y);
    }

}

Canny::~Canny() {

	delete [] data;
	delete [] idata;
	delete [] magnitude;
	delete [] xConv;
	delete [] yConv;
	delete [] xGradient;
	delete [] yGradient;
}

CImg<unsigned char> Canny::toGrayScale(int width, int height) {

	CImg<unsigned char> gray(width, height, 1, 1);

    cimg_forXY(gray, x, y) {
    	double newValue = (img(x,y,0) * 0.2126 + img(x,y,1) * 0.7152 + img(x,y,2) * 0.0722);
    	gray(x,y) = int(newValue);
    }

    return gray;
}

void Canny::performcanny() {

	performcannyparam(2.5f, 7.5f, 5.0f, 16, 0);

}

void Canny::performcannyparam(float lowThreshold, float highthreshold, 
						float gaussiankernelradius, int gaussiankernelwidth,  
						int contrastnormalised) {


	if (contrastnormalised)
	    normalizeContrast();
    computeGradients(gaussiankernelradius, gaussiankernelwidth);

    int low = (int) (lowThreshold * MAGNITUDE_SCALE + 0.5f);
	int high = (int) (highthreshold * MAGNITUDE_SCALE + 0.5f);
	performHysteresis(low, high);

    cimg_forXY(edge, x, y) {
    	edge(x,y) = int(idata[x + y * width]);
    	edge(x, y) = edge(x, y) > 0 ? 255 : 0;
    }

}

int Canny::computeGradients(float kernelRadius, int kernelWidth) {

	float *kernel;
	float *diffKernel;
	int kwidth;

	int initX;
	int maxX;
	int initY;
	int maxY;

	int x, y;
	int i;
	int flag;

	kernel = new float[kernelWidth];
	diffKernel = new float[kernelWidth];

	/* initialise the Gaussian kernel */
	for (kwidth = 0; kwidth < kernelWidth; kwidth++) {
		float g1, g2, g3;
		g1 = gaussian((float) kwidth, kernelRadius);
		if (g1 <= GAUSSIAN_CUT_OFF && kwidth >= 2) 
			break;
		g2 = gaussian(kwidth - 0.5f, kernelRadius);
		g3 = gaussian(kwidth + 0.5f, kernelRadius);
		kernel[kwidth] = (g1 + g2 + g3) / 3.0f / (2.0f * (float) 3.14 * kernelRadius * kernelRadius);
		diffKernel[kwidth] = g3 - g2;
	}

	initX = kwidth - 1;
	maxX = width - (kwidth - 1);
	initY = width * (kwidth - 1);
	maxY = width * (height - (kwidth - 1));
	
	/* perform convolution in x and y directions */
	for(x = initX; x < maxX; x++) {
		for(y = initY; y < maxY; y += width) {
			int index = x + y;
			float sumX = data[index] * kernel[0];
			float sumY = sumX;
			int xOffset = 1;
			int yOffset = width;
			while(xOffset < kwidth) {
				sumY += kernel[xOffset] * (data[index - yOffset] + data[index + yOffset]);
				sumX += kernel[xOffset] * (data[index - xOffset] + data[index + xOffset]);
				yOffset += width;
				xOffset++;
			}
			
			yConv[index] = sumY;
			xConv[index] = sumX;
		}
 	}
 
	for (x = initX; x < maxX; x++) {
		for (y = initY; y < maxY; y += width) {
			float sum = 0.0f;
			int index = x + y;
			for (i = 1; i < kwidth; i++)
				sum += diffKernel[i] * (yConv[index - i] - yConv[index + i]);
 
			xGradient[index] = sum;
		}
 
	}

	for(x = kwidth; x < width - kwidth; x++) {
		for (y = initY; y < maxY; y += width) {
			float sum = 0.0f;
			int index = x + y;
			int yOffset = width;
			for (i = 1; i < kwidth; i++) {
				sum += diffKernel[i] * (xConv[index - yOffset] - xConv[index + yOffset]);
				yOffset += width;
			}
 
			yGradient[index] = sum;
		}
 
	}
 
	initX = kwidth;
	maxX = width - kwidth;
	initY = width * kwidth;
	maxY = width * (height - kwidth);
	for(x = initX; x < maxX; x++) {
		for(y = initY; y < maxY; y += width) {
			int index = x + y;
			int indexN = index - width;
			int indexS = index + width;
			int indexW = index - 1;
			int indexE = index + 1;
			int indexNW = indexN - 1;
			int indexNE = indexN + 1;
			int indexSW = indexS - 1;
			int indexSE = indexS + 1;
			
			float xGrad = xGradient[index];
			float yGrad = yGradient[index];
			float gradMag = hypotenuse(xGrad, yGrad);

			/* perform non-maximal supression */
			float nMag = hypotenuse(xGradient[indexN], yGradient[indexN]);
			float sMag = hypotenuse(xGradient[indexS], yGradient[indexS]);
			float wMag = hypotenuse(xGradient[indexW], yGradient[indexW]);
			float eMag = hypotenuse(xGradient[indexE], yGradient[indexE]);
			float neMag = hypotenuse(xGradient[indexNE], yGradient[indexNE]);
			float seMag = hypotenuse(xGradient[indexSE], yGradient[indexSE]);
			float swMag = hypotenuse(xGradient[indexSW], yGradient[indexSW]);
			float nwMag = hypotenuse(xGradient[indexNW], yGradient[indexNW]);
			float tmp;

			flag = ( (xGrad * yGrad <= 0.0f) /*(1)*/
				? ffabs(xGrad) >= ffabs(yGrad) /*(2)*/
					? (tmp = ffabs(xGrad * gradMag)) >= ffabs(yGrad * neMag - (xGrad + yGrad) * eMag) /*(3)*/
						&& tmp > fabs(yGrad * swMag - (xGrad + yGrad) * wMag) /*(4)*/
					: (tmp = ffabs(yGrad * gradMag)) >= ffabs(xGrad * neMag - (yGrad + xGrad) * nMag) /*(3)*/
						&& tmp > ffabs(xGrad * swMag - (yGrad + xGrad) * sMag) /*(4)*/
				: ffabs(xGrad) >= ffabs(yGrad) /*(2)*/
					? (tmp = ffabs(xGrad * gradMag)) >= ffabs(yGrad * seMag + (xGrad - yGrad) * eMag) /*(3)*/
						&& tmp > ffabs(yGrad * nwMag + (xGrad - yGrad) * wMag) /*(4)*/
					: (tmp = ffabs(yGrad * gradMag)) >= ffabs(xGrad * seMag + (yGrad - xGrad) * sMag) /*(3)*/
						&& tmp > ffabs(xGrad * nwMag + (yGrad - xGrad) * nMag) /*(4)*/
				);
            if(flag){
				magnitude[index] = (gradMag >= MAGNITUDE_LIMIT) ? MAGNITUDE_MAX : (int) (MAGNITUDE_SCALE * gradMag);
			} else {
				magnitude[index] = 0;
			}
		}
	}

	delete [] kernel;
	delete [] diffKernel;

	return 0;

}

void Canny::performHysteresis(int low, int high) {
	
    for (int i = 0; i < width * height; i++)
  	    idata[i] = 0;

  	int offset = 0;
    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {
        if(idata[offset] == 0 && magnitude[offset] >= high) 
	        follow(x, y, offset, low);
	    offset++;
        }
    }
}

void Canny::follow(int x1, int y1, int i1, int threshold) {

    int x, y;
    int x0 = x1 == 0 ? x1 : x1 - 1;
    int x2 = x1 == width - 1 ? x1 : x1 + 1;
    int y0 = y1 == 0 ? y1 : y1 - 1;
    int y2 = y1 == height -1 ? y1 : y1 + 1;
		
    idata[i1] = magnitude[i1];
    for (x = x0; x <= x2; x++) {
        for (y = y0; y <= y2; y++) {
            int i2 = x + y * width;
	        if ((y != y1 || x != x1) && idata[i2] == 0 && magnitude[i2] >= threshold) 
	            follow(x, y, i2, threshold);
        }
    }

}

void Canny::normalizeContrast() {

	int histogram[256] = {0};
    int remap[256];

	for (int i = 0; i < width * height; i++) 
		histogram[data[i]]++;
		
	int sum = 0;
	int j = 0;
	int target;
    for (int i = 0; i < 256; i++) {
		sum += histogram[i];
		target = (sum*255)/(width * height);
		for (int k = j+1; k <= target; k++) 
			remap[k] = i;
		j = target;
	 }
		
    for (int i = 0; i < width * height; i++) 
		data[i] = remap[data[i]];

}

float Canny::hypotenuse(float x, float y) {

	return (float) sqrt(x*x +y*y);

}

float Canny::gaussian(float x, float sigma) {

	return (float) exp(-(x * x) / (2.0f * sigma * sigma));
	
}
