
#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#include <filter2d.h>


void filter2d_mean3(double* img_in, long width, long height, double* img_out)
{

	double H[3][3] = { {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0},
					   {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0},
					   {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0} };

	for (int i = 0; i < width * height; i++) img_out[i] = img_in[i];

	for (long h = 1; h < height - 1; h++) {
		for (long w = 1; w < width - 1; w++) {
			double pix = 0.0;
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					pix += H[l][k] * img_in[(w - k + 1) + width * (h - l + 1)];
				}
			}
			img_out[w + width * h] = pix;
		}
	}

}


void sobel_norm(double* img_in, long width, long height, double* img_out)
{

	double Sx[3][3] = { {1.0 / 8.0, 0.0 / 8.0, -1.0 / 8.0},
					   {2.0 / 8.0, 0.0 / 8.0, -2.0 / 8.0},
					   {1.0 / 8.0, 0.0 / 8.0, -1.0 / 8.0} };



	int k, l, i;
	double sx, sy;

	for (i = 0; i < width * height; i++) img_out[i] = img_in[i];

	for (long h = 1; h < height - 1; h++) {
		for (long w = 1; w < width - 1; w++) {
			sx = sy = 0.0;
			for (k = -1; k < 1; k++) {
				for (l = -1; l < 1; l++) {
					sx += img_in[(w - k) + width * (h - l)] * Sx[l + 1][k + 1];
					sy += img_in[(w - k) + width * (h - l)] * Sx[k + 1][l + 1];
				}
			}
			img_out[w + width * h] = sqrt(sx * sx + sy * sy);
		}
	}

}


static double* gaussian2d_create(double sigma, long* size)
{
	double* filter = NULL, s = 0.0, v;

	int t, i, j;

	t = 1 + 2 * (int)ceil(3.0 * sigma);
	*size = t;
	filter = (double*)calloc(t * t, sizeof(double));

	/*
	for (i = 0; i < t ; i++) {
		x = (double)(i - t / 2) / sigma;
		for (j = 0; j < t ; j++) {
			y = (double)(j - t / 2) / sigma;
			v = exp(  -(x*x + y*y)/2.0 );
			s + = v;
			filter[j * t + i] = v;
		}
	} */


	for (i = -t / 2; i < t / 2; i++) {
		for (j = -t / 2; j < t / 2; j++) {
			v = exp(-(double)(i * i + j * j) / (2.0 * sigma * sigma));
			s += v;
			filter[(j + t / 2) * t + i + t / 2] = v;
		}
	}

	for (i = 0; i < t; i++)
		for (j = 0; j < t; j++)
			filter[j * t + i] /= s;

	return filter;
}


static void convolution2d(double* img_in, long width, long height, double* mask2d, long tx, long ty, double* img_out)
{
	int x, y, k, l, i;
	double s;

	for (i = 0; i < width * height; i++) img_out[i] = 0.0;

	for (x = tx / 2; x < width - tx / 2; x++) {
		for (y = ty / 2; y < height - ty / 2; y++) {
			s = 0.0;

			for (k = -tx / 2; k <= tx / 2; k++) {
				for (l = -ty / 2; l <= ty / 2; l++) {
					s += img_in[(y - l) * width + x - k] * mask2d[(l + ty / 2) * tx + k + tx / 2];
				}
			}

			/*
			for (k = 0 ; k < tx ; k++) {
				for (l = 0; l < ty; l++) {
					s += img_in[(y - (l-ty/2) * width + x - (k-tx/2)] * mask2d[l * tx + k];
				}
			}*/

			img_out[y * width * x] = s;
		}
	}
}

static double** gaussian2d_create_matrix(double sigma, long* size)
{
	return NULL;
}


static void convolution2d_by_matrix(double* img_in, long width, long height,
	double** mask2d, long tx, long ty,
	double* img_out)
{
}


void filter2d_gaussian(double* img_in, long width, long height, double sigma, double* img_out)
{
	double* filter = NULL;
	long size;

	filter = gaussian2d_create(sigma, &size);
	convolution2d(img_in, width, height, filter, size, size, img_out);
	free(filter);
}


void img_get_raw(double* img, long width, long height,
	long no,
	double* v)
{
}


void img_set_raw(double* img, long width, long height,
	long no,
	double* v)
{
}


void img_get_column(double* img, long width, long height,
	long no,
	double* v)
{
}


void img_set_column(double* img, long width, long height,
	long no,
	double* v)
{
}


static double* gaussian1d_create(double sigma, long* size)
{
	return NULL;
}


static void convolution1d(double* v_in, long size,
	double* mask1d, long t,
	double* v_out)
{
}


void filter2d_gaussian_fast(double* img_in, long width, long height,
	double sigma,
	double* img_out)
{
}


typedef double (*PROC)(double*, long);


static void filter2d_generic(double* img_in, long width, long height,
	long tx, long ty, PROC proc,
	double* img_out)
{
}


void filter2d_by_method(double* img_in, long width, long height,
	long tx, long ty, long method,
	double* img_out)
{
}


