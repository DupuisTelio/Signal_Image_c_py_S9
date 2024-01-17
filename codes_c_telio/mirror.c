
#include <mirror.h>


void img_mirror_horizontal(double * img_in, long width, long height,
													 double * img_out)
{
	int i, j;
	// première version triviale
	/*
	for (i = 0; i < width; i++) {
		for (j = 0; j < height; j++) {
			img_out[width * j + width - 1 - i] = img_in[width * j + i];
		}
	}
	*/


	// deuxième version arithmétique pointeurs
	img_out += width - 1;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			*(img_out--) = *img_in++;
		}
		img_out += 2 * width;
	}

}