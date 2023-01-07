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
static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	pixel_sum sum;
	pixel current_pixel;
	int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
	int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
	int min_row, min_col, max_row, max_col;
	pixel loop_pixel;

	initialize_pixel_sum(&sum);

    int startIndexI = i-1;
    int startIndexJ = j-1;

    int globali, globalj, locali, localj, ii;

//    for (ii = 0; ii < 9; ii++) {
//        locali = ii/3;
//        localj = ii%3;
//
//        globali = startIndexI+locali;
//        globalj = startIndexJ+localj;
//
//        sum_pixels_by_weight(&sum, src[globali*dim + globalj], kernel[locali][localj]);
//    }

    for (ii = 0; ii < 9; ii++) {
        int locali = ii/3;
        int localj = ii%3;

        int globali = startIndexI+locali;
        int globalj = startIndexJ+localj;

        int kernel_val;
        pixel src_pixel;
        memcpy(&kernel_val, &kernel[locali][localj], sizeof(int));
        memcpy(&src_pixel, &src[globali*dim + globalj], sizeof(pixel));

        sum_pixels_by_weight(&sum, src_pixel, kernel_val);
    }

    if (filter) {
        for (ii = 0; ii < 9; ii++) {
            int locali = ii/3;
            int localj = ii%3;

            int globali = startIndexI+locali;
            int globalj = startIndexJ+localj;

            pixel loop_pixel;
            memcpy(&loop_pixel, &src[globali * dim + globalj], sizeof(pixel));

            int intensity;
            intensity = ((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue);
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

        // filter out min and max
        sum_pixels_by_weight(&sum, src[min_row * dim + min_col], -1);
        sum_pixels_by_weight(&sum, src[max_row * dim + max_col], -1);
    }

//	if (filter) {
//        for (ii = 0; ii < 9; ii++) {
//            locali = ii/3;
//            localj = ii%3;
//
//            globali = startIndexI+locali;
//            globalj = startIndexJ+localj;
//
//            loop_pixel = src[globali * dim + globalj];
//
//            int intensity = ((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue);
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
//		// filter out min and max
//        sum_pixels_by_weight(&sum, src[min_row * dim + min_col], -1);
//        sum_pixels_by_weight(&sum, src[max_row * dim + max_col], -1);
//	}

	// assign kernel's result to pixel at [i,j]
	assign_sum_to_pixel(&current_pixel, sum, kernelScale);
	return current_pixel;
}

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	int i, j;
	for (i = kernelSize / 2 ; i < dim - kernelSize / 2; i++) {
        int rowStart = i * dim;
		for (j =  kernelSize / 2 ; j < dim - kernelSize / 2 ; j++) {
			dst[rowStart + j] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
		}
	}
}

void charsToPixels(Image *charsImg, pixel* pixels) {

    size_t numPixels = m * n;
    size_t pixelSize = sizeof(pixel);

    memcpy(pixels, charsImg->data, numPixels * pixelSize);

//    int numChars = 3*m*n;
//    char pix[numChars];
//    memcpy(pix, charsImg->data, numChars);
//
//	int row, col;
//	for (row = 0 ; row < m ; row++) {
//        int rowStart = row * n;
//		for (col = 0 ; col < n ; col++) {
//            int index = rowStart + col;
//            int imIndex = 3 * index;
//            pixels[index].red = pix[imIndex];
//            pixels[index].green = pix[imIndex + 1];
//            pixels[index].blue = pix[imIndex + 2];
//		}
//	}

}

void pixelsToChars(pixel* pixels, Image *charsImg) {

    size_t numPixels = m * n;
    size_t pixelSize = sizeof(pixel);

    memcpy(charsImg->data, pixels, numPixels * pixelSize);

//    int numChars = sizeof(pixel)*m*n;
//    pixel pix[numChars];
//    memcpy(pix, pixels, numChars);
//
//    int row, col;
//	for (row = 0 ; row < m ; row++) {
//        int rowStart = row * n;
//		for (col = 0 ; col < n ; col++) {
//            int index = rowStart + col;
//            int imIndex = 3 * index;
//			image->data[imIndex] = pix[index].red;
//			image->data[imIndex + 1] = pix[index].green;
//			image->data[imIndex + 2] = pix[index].blue;
//		}
//	}
}

void copyPixels(pixel* src, pixel* dst) {

    size_t numPixels = m * n;
    size_t pixelSize = sizeof(pixel);

    memcpy(dst, src, numPixels * pixelSize);

//	int row, col;
//	for (row = 0 ; row < m ; row++) {
//        int rowStart = row * n;
//		for (col = 0 ; col < n ; col++) {
//            int index = rowStart + col;
//			dst[index].red = src[index].red;
//			dst[index].green = src[index].green;
//			dst[index].blue = src[index].blue;
//		}
//	}
}

void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	pixel* pixelsImg = malloc(m*n*sizeof(pixel));
	pixel* backupOrg = malloc(m*n*sizeof(pixel));

	charsToPixels(image, pixelsImg);
	copyPixels(pixelsImg, backupOrg);

	smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter);

	pixelsToChars(pixelsImg, image);

	free(pixelsImg);
	free(backupOrg);
}

