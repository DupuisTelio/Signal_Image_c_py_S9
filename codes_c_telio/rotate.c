

#include <rotate.h>

/*
void* myLibInit();
void* myLibAddImg(void*, double* img, int width, int height);
void* myLibRun(void*);
double* myLibGetOut(void*);
void* myLibRelease(void*);
*/

typedef struct _RGBA_
{
	unsigned char r;
	unsigned char g;
	unsigned char b;
	unsigned char a;
} RGBA;


void img_rotate(void * img_in, long width, long height,
								long type_data,
								void * img_out)
{
	if (type_data == 0) {
		int i, j;
		double* img_data_in = (double*)img_in, *img_data_out = (double*)img_out;

		for (j = 0; j < height; j++) {
			for (i = 0; i < width; i++) 
				img_data_out[i*height + height-1-j] = img_data_in[j * width + i];
				//img_data_out[(width - 1 - i) * height + j] = img_data_in[j * width + i];
				//img_data_out[i * height + j] = img_data_in[j * width + i];
		}
		/*for (j = 0; j < height; j++) {
			for (i = 0; i < width; i++)
				img_data_out[i* height + j] = img_data_in[(height-1-i)*width+j];
		}*/

	}

	else if (type_data == 1) {
		int i, j;
		RGBA* img_data_in = (RGBA*)img_in, * img_data_out = (RGBA*)img_out;

		for (j = 0; j < height; j++) {
			for (i = 0; i < width; i++)
				img_data_out[i * height + height - 1 - j] = img_data_in[j * width + i];
		}
	}

	if (type_data == 0 || type_data == 1) {
		int i, j,k,t;
		unsigned char * img_data_in = (unsigned char*)img_in, * img_data_out = (unsigned char*)img_out;

		if (type_data == 0) t = 8; else t = 4;

		for (j = 0; j < height; j++) 
			for (i = 0; i < width; i++)
				for (k=0;k<t;k++)
					img_data_out[(i * height + height - 1 - j)*t+k] = img_data_in[(j * width + i)*t+k];
	}
}                