

#include <cut.h>

typedef struct _RGBA_
{
	unsigned char r;
	unsigned char g;
	unsigned char b;
	unsigned char a; 
} RGBA;

void img_extract_area(unsigned char * img_color_in, long width, long height,
											long x, long y, long dx, long dy,
											unsigned char * img_color_out)
{
	int i, j;
	//for (j=0;j<dy;j++){
	//	for (i = 0; i < dx; i++) {
	//		img_color_out[(j * dx + i) * 4] = img_color_in[((y+j)*width+x+i)*4];
	//		img_color_out[(j * dx + i) * 4+1] = img_color_in[((y + j) * width + x + i) * 4+1]; // à compléter prochaine séance
	//		img_color_out[(j * dx + i) * 4+2] = img_color_in[((y + j) * width + x + i) * 4+2];
	//	}
	//}


	// version optimisée
	RGBA* img_rgba_in= (RGBA *)img_color_in, * img_rgba_out = (RGBA*)img_color_out;

	for (j = 0; j < dy; j++) {
		for (i = 0; i < dx; i++) {
			img_rgba_out[j * dx + i] = img_rgba_in[(y + j) * width + x + i];
		}
	}

	// version optimisée ++
	img_rgba_in += y*width+x;
	for (j = 0; j < dy; j++) {
		for (i = 0; i < dx; i++) *(img_rgba_out++) = *(img_rgba_in++);
		img_rgba_in += width-dx;
	}
}