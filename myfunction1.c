#include <stdlib.h>
#include <stdbool.h>
#include "myfunction1.h"
#include "showBMP.h"
#include <string.h>
#include <stdio.h>

float THIRD = 1.0/3;

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */

void initialize_pixel_sum(pixel_sum *sum) {
	sum->red = sum->green = sum->blue = 0;
	// sum->num = 0;
	return;
}

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, float kernelScale) {

    // divide by kernel's weight
    sum.red = (int) (sum.red * kernelScale);
    sum.green = (int) (sum.green * kernelScale);
    sum.blue = (int) (sum.blue * kernelScale);

    // truncate each pixel's color values to match the range [0,255]
    current_pixel->red = (unsigned char) (sum.red < 0 ? 0 : (sum.red > 255 ? 255 : sum.red));
    current_pixel->green = (unsigned char) (sum.green < 0 ? 0 : (sum.green > 255 ? 255 : sum.green));
    current_pixel->blue = (unsigned char) (sum.blue < 0 ? 0 : (sum.blue > 255 ? 255 : sum.blue));
    return;

}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, float weight) {
	sum->red += (int) (p.red * weight);
	sum->green += (int) (p.green * weight);
	sum->blue += (int) (p.blue * weight);
	// sum->num++;
	return;
}

static pixel applyKernel1(int dim, int i, int j, pixel *src, float kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    sum.red = 0;
    sum.blue = 0;
    sum.green = 0;

    int startIndexI = (i-1)*dim;
    int startIndexJ = j-1;

    int ii;

    pixel tmp;

    for (ii = 0; ii < 3; ii++) {

        tmp = src[startIndexI + startIndexJ + ii];

        sum.red += tmp.red;
        sum.blue += tmp.blue;
        sum.green += tmp.green;
    }

    startIndexI += dim;

    for (ii = 0; ii < 3; ii++) {

        tmp = src[startIndexI + startIndexJ + ii];

        sum.red += tmp.red;
        sum.blue += tmp.blue;
        sum.green += tmp.green;
    }

    startIndexI += dim;

    for (ii = 0; ii < 3; ii++) {

        tmp = src[startIndexI + startIndexJ + ii];

        sum.red += tmp.red;
        sum.blue += tmp.blue;
        sum.green += tmp.green;
    }


    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}
static pixel applyKernel2(int dim, int i, int j, pixel *src, float kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    sum.red = 0;
    sum.blue = 0;
    sum.green = 0;

    int startIndexI = (i-1)*dim;
    int startIndexJ = j-1;

    int ii;

    pixel tmp;

    for (ii = 0; ii < 3; ii++) {

        tmp = src[startIndexI + startIndexJ + ii];

        sum.red -= tmp.red;
        sum.blue -= tmp.blue;
        sum.green -= tmp.green;
    }

    startIndexI += dim;


    tmp = src[startIndexI + startIndexJ];

    sum.red -= tmp.red;
    sum.blue -= tmp.blue;
    sum.green -= tmp.green;

    tmp = src[startIndexI + startIndexJ + 1];

    sum.red += tmp.red * 9;
    sum.blue += tmp.blue * 9;
    sum.green += tmp.green * 9;

    tmp = src[startIndexI + startIndexJ + 2];

    sum.red -= tmp.red;
    sum.blue -= tmp.blue;
    sum.green -= tmp.green;


    startIndexI += dim;

    for (ii = 0; ii < 3; ii++) {

        tmp = src[startIndexI + startIndexJ + ii];

        sum.red -= tmp.red;
        sum.blue -= tmp.blue;
        sum.green -= tmp.green;
    }

    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}
static pixel applyKernel3(int dim, int i, int j, pixel *src, float kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    sum.red = 0;
    sum.blue = 0;
    sum.green = 0;

    int startIndex = i*dim + j-1;

    pixel tmp;

    tmp = src[startIndex];

    sum.red += tmp.red;
    sum.blue += tmp.blue;
    sum.green += tmp.green;

    tmp = src[startIndex + 1];

    sum.red += tmp.red * 2;
    sum.blue += tmp.blue * 2;
    sum.green += tmp.green * 2;

    tmp = src[startIndex + 2];

    sum.red += tmp.red;
    sum.blue += tmp.blue;
    sum.green += tmp.green;

    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}
static pixel applyKernel4(int dim, int i, int j, pixel *src, float kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    sum.red = 0;
    sum.blue = 0;
    sum.green = 0;

    int startIndex = i*dim + j-1;

    pixel tmp;

    tmp = src[startIndex];

    sum.red += tmp.red * (-2);
    sum.blue += tmp.blue * (-2);
    sum.green += tmp.green * (-2);

    tmp = src[startIndex + 1];

    sum.red += tmp.red * 6;
    sum.blue += tmp.blue * 6;
    sum.green += tmp.green * 6;

    tmp = src[startIndex + 2];

    sum.red += tmp.red * (-2);
    sum.blue += tmp.blue * (-2);
    sum.green += tmp.green * (-2);

    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}
static pixel applyKernelFilter(int dim, int i, int j, pixel *src, float kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
    int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
    int min_row, min_col, max_row, max_col;

    int startIndexI = (i-1)*dim;
    int startIndexJ = j-1;

    for (int ii = 0; ii < 3; ii++) {

        pixel loop_pixel = src[startIndexI + startIndexJ + ii];
        int intensity = loop_pixel.red + loop_pixel.green + loop_pixel.blue;

        sum.red += loop_pixel.red;
        sum.blue += loop_pixel.blue;
        sum.green += loop_pixel.green;

        min_row = (intensity <= min_intensity) ? startIndexI : min_row;
        min_col = (intensity <= min_intensity) ? startIndexJ + ii : min_col;
        min_intensity = (intensity <= min_intensity) ? intensity : min_intensity;

        max_row = (intensity > max_intensity) ? startIndexI : max_row;
        max_col = (intensity > max_intensity) ? startIndexJ + ii : max_col;
        max_intensity = (intensity > max_intensity) ? intensity : max_intensity;

    }

    startIndexI += dim;

    for (int ii = 0; ii < 3; ii++) {

        pixel loop_pixel = src[startIndexI + startIndexJ + ii];
        int intensity = loop_pixel.red + loop_pixel.green + loop_pixel.blue;

        sum.red += loop_pixel.red;
        sum.blue += loop_pixel.blue;
        sum.green += loop_pixel.green;

        if (intensity <= min_intensity) {
            min_intensity = intensity;
            min_row = startIndexI;
            min_col = startIndexJ + ii;
        }
        if (intensity > max_intensity) {
            max_intensity = intensity;
            max_row = startIndexI;
            max_col = startIndexJ + ii;
        }
    }

    startIndexI += dim;

    for (int ii = 0; ii < 3; ii++) {

        pixel loop_pixel = src[startIndexI + startIndexJ + ii];
        int intensity = loop_pixel.red + loop_pixel.green + loop_pixel.blue;

        sum.red += loop_pixel.red;
        sum.blue += loop_pixel.blue;
        sum.green += loop_pixel.green;

        if (intensity <= min_intensity) {
            min_intensity = intensity;
            min_row = startIndexI;
            min_col = startIndexJ + ii;
        }
        if (intensity > max_intensity) {
            max_intensity = intensity;
            max_row = startIndexI;
            max_col = startIndexJ + ii;
        }
    }

    pixel min = src[min_row + min_col];
    pixel max = src[max_row + max_col];

    sum.red += -min.red - max.red;
    sum.blue += -min.blue - max.blue;
    sum.green += -min.green - max.green;


    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth1(int dim, pixel *src, pixel *dst, float kernelScale) {

	int i, j;
	for (i = 1 ; i < dim - 1; i++) {
        int rowStart = i * dim;
		for (j =  1 ; j < dim - 1 ; j++) {
			dst[rowStart + j] = applyKernel1(dim, i, j, src, kernelScale);
		}
	}
}
void smooth2(int dim, pixel *src, pixel *dst, float kernelScale) {

    int i, j;
    for (i = 1 ; i < dim - 1; i++) {
        int rowStart = i * dim;
        for (j =  1 ; j < dim - 1 ; j++) {
            dst[rowStart + j] = applyKernel2(dim, i, j, src, kernelScale);
        }
    }
}
void smooth3(int dim, pixel *src, pixel *dst, float kernelScale) {

    int i, j;
    for (i = 1 ; i < dim - 1; i++) {
        int rowStart = i * dim;
        for (j =  1 ; j < dim - 1 ; j++) {
            dst[rowStart + j] = applyKernel3(dim, i, j, src, kernelScale);
        }
    }
}
void smooth4(int dim, pixel *src, pixel *dst, float kernelScale) {

    int i, j;
    for (i = 1 ; i < dim - 1; i++) {
        int rowStart = i * dim;
        for (j =  1 ; j < dim - 1 ; j++) {
            dst[rowStart + j] = applyKernel4(dim, i, j, src, kernelScale);
        }
    }
}
void smoothfilter(int dim, pixel *src, pixel *dst, float kernelScale) {

    int i, j;
    for (i = 1 ; i < dim - 1; i++) {
        int rowStart = i * dim;
        for (j =  1 ; j < dim - 1 ; j++) {
            dst[rowStart + j] = applyKernelFilter(dim, i, j, src, kernelScale);
        }
    }
}

void charsToPixels(Image *charsImg, pixel* pixels) {

    size_t numPixels = m * n;
    size_t pixelSize = sizeof(pixel);

    memcpy(pixels, charsImg->data, numPixels * pixelSize);

}

void pixelsToChars(pixel* pixels, Image *charsImg) {

    size_t numPixels = m * n;
    size_t pixelSize = sizeof(pixel);

    memcpy(charsImg->data, pixels, numPixels * pixelSize);

}

void copyPixels(pixel* src, pixel* dst) {

    size_t numPixels = m * n;
    size_t pixelSize = sizeof(pixel);

    memcpy(dst, src, numPixels * pixelSize);

}

void doConvolution(Image *image, float kernelScale, bool filter, int kernelNum) {
    pixel* pixelsImg = (pixel*) image->data;
    pixel* backupOrg = malloc(m*n*sizeof(pixel));

    copyPixels(pixelsImg, backupOrg);

    if (!filter) {
        if (kernelNum == 1) {
            smooth1(m, backupOrg, pixelsImg, kernelScale);
        } else if (kernelNum == 2) {
            smooth2(m, backupOrg, pixelsImg, kernelScale);
        } else if (kernelNum == 3) {
            smooth3(m, backupOrg, pixelsImg, kernelScale);
        } else if (kernelNum == 4) {
            smooth4(m, backupOrg, pixelsImg, kernelScale);
        }
    } else {
        smoothfilter(m, backupOrg, pixelsImg, kernelScale);
    }

    free(backupOrg);
}

