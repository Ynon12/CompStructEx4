
#ifndef __myfunc1_h
#define __myfunc1_h

#include "myutil.h"
#include "readBMP.h"

void initialize_pixel_sum(pixel_sum *sum);
void doConvolution(Image *image, float kernelScale, bool filter, int kernelNum);

#endif

