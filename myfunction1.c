#include <stdlib.h>
#include <stdbool.h>
#include "myfunction1.h"
#include "showBMP.h"
#include <string.h>
#include <stdio.h>
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
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

    // divide by kernel's weight
    sum.red = sum.red / kernelScale;
    sum.green = sum.green / kernelScale;
    sum.blue = sum.blue / kernelScale;

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

/*
 *  Applies kernel for pixel at (i,j)
 */
//static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter, int kernelNum) {
//
//	pixel_sum sum;
//	pixel current_pixel;
//	int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
//	int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
//	int min_row, min_col, max_row, max_col;
//	pixel loop_pixel;
//
//	initialize_pixel_sum(&sum);
//
//    int startIndexI = i-1;
//    int startIndexJ = j-1;
//
//    int globali, globalj, locali, localj, ii;
//
//    pixel tmp;
//
//    if (kernelNum == 1 && !filter) {
//
//        int sumRed = 0;
//        int sumBlue = 0;
//        int sumGreen = 0;
//
//        for (ii = 0; ii < 9; ii++) {
//            locali = ii/3;
//            localj = ii%3;
//
//            globali = startIndexI+locali;
//            globalj = startIndexJ+localj;
//
//            tmp = src[globali*dim + globalj];
//
//            sumRed += tmp.red;
//            sumBlue += tmp.blue;
//            sumGreen += tmp.green;
//        }
//
//        sum.red = sumRed;
//        sum.blue = sumBlue;
//        sum.green = sumGreen;
//
//    } else if (kernelNum == 2 && !filter) {
//
//        int sumRed = 0;
//        int sumBlue = 0;
//        int sumGreen = 0;
//
//        for (ii = 0; ii < 9; ii++) {
//
//            locali = ii/3;
//            localj = ii%3;
//
//            globali = startIndexI+locali;
//            globalj = startIndexJ+localj;
//
//            tmp = src[globali*dim + globalj];
//
//            if (ii == 4) {
//                sumRed += tmp.red * 9;
//                sumBlue += tmp.blue * 9;
//                sumGreen += tmp.green * 9;
//            } else {
//                sumRed -= tmp.red;
//                sumBlue -= tmp.blue;
//                sumGreen -= tmp.green;
//            }
//        }
//
//        sum.red = sumRed;
//        sum.blue = sumBlue;
//        sum.green = sumGreen;
//
//    } else if (kernelNum == 3 && !filter) {
//
//        int sumRed = 0;
//        int sumBlue = 0;
//        int sumGreen = 0;
//
//        tmp = src[(startIndexI+1)*dim + startIndexJ];
//
//        sumRed += tmp.red;
//        sumBlue += tmp.blue;
//        sumGreen += tmp.green;
//
//        tmp = src[(startIndexI+1)*dim + startIndexJ + 1];
//
//        sumRed += tmp.red * 2;
//        sumBlue += tmp.blue * 2;
//        sumGreen += tmp.green * 2;
//
//        tmp = src[(startIndexI+1)*dim + startIndexJ + 2];
//
//        sumRed += tmp.red;
//        sumBlue += tmp.blue;
//        sumGreen += tmp.green;
//
//        sum.red = sumRed;
//        sum.blue = sumBlue;
//        sum.green = sumGreen;
//
//    } else if (kernelNum == 4 && !filter) {
//
//        int sumRed = 0;
//        int sumBlue = 0;
//        int sumGreen = 0;
//
//        tmp = src[(startIndexI+1)*dim + startIndexJ];
//
//        sumRed += tmp.red * (-2);
//        sumBlue += tmp.blue * (-2);
//        sumGreen += tmp.green * (-2);
//
//        tmp = src[(startIndexI+1)*dim + startIndexJ + 1];
//
//        sumRed += tmp.red * 6;
//        sumBlue += tmp.blue * 6;
//        sumGreen += tmp.green * 6;
//
//        tmp = src[(startIndexI+1)*dim + startIndexJ + 2];
//
//        sumRed += tmp.red * (-2);
//        sumBlue += tmp.blue * (-2);
//        sumGreen += tmp.green * (-2);
//
//        sum.red = sumRed;
//        sum.blue = sumBlue;
//        sum.green = sumGreen;
//
//    }
//
//    if (filter) {
//        for (int ii = 0; ii < 9; ii++) {
//            int globali = startIndexI + ii/3;
//            int globalj = startIndexJ + ii%3;
//
//            pixel loop_pixel = src[globali * dim + globalj];
//            int intensity = loop_pixel.red + loop_pixel.green + loop_pixel.blue;
//
//            sum.red += loop_pixel.red;
//            sum.blue += loop_pixel.blue;
//            sum.green += loop_pixel.green;
//
//
//            if (intensity <= min_intensity) {
//                min_intensity = intensity;
//                min_row = globali;
//                min_col = globalj;
//            }
//            if (intensity > max_intensity) {
//                max_intensity = intensity;
//                max_row = globali;
//                max_col = globalj;
//            }
//        }
//
//        pixel min = src[min_row*dim + min_col];
//        pixel max = src[max_row*dim + max_col];
//
//        sum.red += -min.red - max.red;
//        sum.blue += -min.blue - max.blue;
//        sum.green += -min.green - max.green;
//
//    }
//
//	// assign kernel's result to pixel at [i,j]
//	assign_sum_to_pixel(&current_pixel, sum, kernelScale);
//	return current_pixel;
//}

static pixel applyKernel1(int dim, int i, int j, pixel *src, int kernelSize, int kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    sum.red = 0;
    sum.blue = 0;
    sum.green = 0;

    int startIndexI = i-1;
    int startIndexJ = j-1;

    int globali, globalj, locali, localj, ii;

    pixel tmp;

//    int sumRed = 0;
//    int sumBlue = 0;
//    int sumGreen = 0;

    for (ii = 0; ii < 9; ii++) {
        locali = ii/3;
        localj = ii%3;

        globali = startIndexI+locali;
        globalj = startIndexJ+localj;

        tmp = src[globali*dim + globalj];
//
//        sumRed += tmp.red;
//        sumBlue += tmp.blue;
//        sumGreen += tmp.green;

        sum.red += tmp.red;
        sum.blue += tmp.blue;
        sum.green += tmp.green;
    }

//    sum.red = sumRed;
//    sum.blue = sumBlue;
//    sum.green = sumGreen;

    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}
static pixel applyKernel2(int dim, int i, int j, pixel *src, int kernelSize,  int kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    sum.red = 0;
    sum.blue = 0;
    sum.green = 0;

    int startIndexI = i-1;
    int startIndexJ = j-1;

    int globali, globalj, locali, localj, ii;

    pixel tmp;

//    int sumRed = 0;
//    int sumBlue = 0;
//    int sumGreen = 0;

    for (ii = 0; ii < 9; ii++) {

        locali = ii/3;
        localj = ii%3;

        globali = startIndexI+locali;
        globalj = startIndexJ+localj;

        tmp = src[globali*dim + globalj];

        if (ii == 4) {
            sum.red += tmp.red * 9;
            sum.blue += tmp.blue * 9;
            sum.green += tmp.green * 9;
        } else {
            sum.red -= tmp.red;
            sum.blue -= tmp.blue;
            sum.green -= tmp.green;
        }
    }

//    sum.red = sumRed;
//    sum.blue = sumBlue;
//    sum.green = sumGreen;


    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}
static pixel applyKernel3(int dim, int i, int j, pixel *src, int kernelSize, int kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    sum.red = 0;
    sum.blue = 0;
    sum.green = 0;

    int startIndexI = i-1;
    int startIndexJ = j-1;

    pixel tmp;

//    int sumRed = 0;
//    int sumBlue = 0;
//    int sumGreen = 0;

    tmp = src[(startIndexI+1)*dim + startIndexJ];

    sum.red += tmp.red;
    sum.blue += tmp.blue;
    sum.green += tmp.green;

    tmp = src[(startIndexI+1)*dim + startIndexJ + 1];

    sum.red += tmp.red * 2;
    sum.blue += tmp.blue * 2;
    sum.green += tmp.green * 2;

    tmp = src[(startIndexI+1)*dim + startIndexJ + 2];

    sum.red += tmp.red;
    sum.blue += tmp.blue;
    sum.green += tmp.green;

//    sum.red = sumRed;
//    sum.blue = sumBlue;
//    sum.green = sumGreen;
//

    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}
static pixel applyKernel4(int dim, int i, int j, pixel *src, int kernelSize, int kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    sum.red = 0;
    sum.blue = 0;
    sum.green = 0;

    int startIndexI = i-1;
    int startIndexJ = j-1;

    pixel tmp;

//    int sumRed = 0;
//    int sumBlue = 0;
//    int sumGreen = 0;

    tmp = src[(startIndexI+1)*dim + startIndexJ];

    sum.red += tmp.red * (-2);
    sum.blue += tmp.blue * (-2);
    sum.green += tmp.green * (-2);

    tmp = src[(startIndexI+1)*dim + startIndexJ + 1];

    sum.red += tmp.red * 6;
    sum.blue += tmp.blue * 6;
    sum.green += tmp.green * 6;

    tmp = src[(startIndexI+1)*dim + startIndexJ + 2];

    sum.red += tmp.red * (-2);
    sum.blue += tmp.blue * (-2);
    sum.green += tmp.green * (-2);

//    sum.red = sumRed;
//    sum.blue = sumBlue;
//    sum.green = sumGreen;

    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}
static pixel applyKernelFilter(int dim, int i, int j, pixel *src, int kernelSize, int kernelScale) {

    pixel_sum sum;
    pixel current_pixel;

    int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
    int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
    int min_row, min_col, max_row, max_col;

    int startIndexI = i-1;
    int startIndexJ = j-1;

    for (int ii = 0; ii < 9; ii++) {
        int globali = startIndexI + ii/3;
        int globalj = startIndexJ + ii%3;

        pixel loop_pixel = src[globali * dim + globalj];
        int intensity = loop_pixel.red + loop_pixel.green + loop_pixel.blue;

        sum.red += loop_pixel.red;
        sum.blue += loop_pixel.blue;
        sum.green += loop_pixel.green;

        if (intensity <= min_intensity) {
            min_intensity = intensity;
            min_row = globali;
            min_col = globalj;
        }
        if (intensity > max_intensity) {
            max_intensity = intensity;
            max_row = globali;
            max_col = globalj;
        }
    }

    pixel min = src[min_row*dim + min_col];
    pixel max = src[max_row*dim + max_col];

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
void smooth1(int dim, pixel *src, pixel *dst, int kernelSize, int kernelScale) {

	int i, j;
	for (i = kernelSize / 2 ; i < dim - kernelSize / 2; i++) {
        int rowStart = i * dim;
		for (j =  kernelSize / 2 ; j < dim - kernelSize / 2 ; j++) {
			dst[rowStart + j] = applyKernel1(dim, i, j, src, kernelSize, kernelScale);
		}
	}
}
void smooth2(int dim, pixel *src, pixel *dst, int kernelSize, int kernelScale) {

    int i, j;
    for (i = kernelSize / 2 ; i < dim - kernelSize / 2; i++) {
        int rowStart = i * dim;
        for (j =  kernelSize / 2 ; j < dim - kernelSize / 2 ; j++) {
            dst[rowStart + j] = applyKernel2(dim, i, j, src, kernelSize, kernelScale);
        }
    }
}
void smooth3(int dim, pixel *src, pixel *dst, int kernelSize, int kernelScale) {

    int i, j;
    for (i = kernelSize / 2 ; i < dim - kernelSize / 2; i++) {
        int rowStart = i * dim;
        for (j =  kernelSize / 2 ; j < dim - kernelSize / 2 ; j++) {
            dst[rowStart + j] = applyKernel3(dim, i, j, src, kernelSize, kernelScale);
        }
    }
}
void smooth4(int dim, pixel *src, pixel *dst, int kernelSize, int kernelScale) {

    int i, j;
    for (i = kernelSize / 2 ; i < dim - kernelSize / 2; i++) {
        int rowStart = i * dim;
        for (j =  kernelSize / 2 ; j < dim - kernelSize / 2 ; j++) {
            dst[rowStart + j] = applyKernel4(dim, i, j, src, kernelSize, kernelScale);
        }
    }
}
void smoothfilter(int dim, pixel *src, pixel *dst, int kernelSize, int kernelScale) {

    int i, j;
    for (i = kernelSize / 2 ; i < dim - kernelSize / 2; i++) {
        int rowStart = i * dim;
        for (j =  kernelSize / 2 ; j < dim - kernelSize / 2 ; j++) {
            dst[rowStart + j] = applyKernelFilter(dim, i, j, src, kernelSize, kernelScale);
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

void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter, int kernelNum) {

	pixel* pixelsImg = malloc(m*n*sizeof(pixel));
	pixel* backupOrg = malloc(m*n*sizeof(pixel));

	charsToPixels(image, pixelsImg);
	copyPixels(pixelsImg, backupOrg);

    if (!filter) {
        if (kernelNum == 1) {
            smooth1(m, backupOrg, pixelsImg, kernelSize, kernelScale);
        } else if (kernelNum == 2) {
            smooth2(m, backupOrg, pixelsImg, kernelSize, kernelScale);
        } else if (kernelNum == 3) {
            smooth3(m, backupOrg, pixelsImg, kernelSize, kernelScale);
        } else if (kernelNum == 4) {
            smooth4(m, backupOrg, pixelsImg, kernelSize, kernelScale);
        }
    } else {
        smoothfilter(m, backupOrg, pixelsImg, kernelSize, kernelScale);
    }
//	smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter, kernelNum);

	pixelsToChars(pixelsImg, image);

	free(pixelsImg);
	free(backupOrg);
}

