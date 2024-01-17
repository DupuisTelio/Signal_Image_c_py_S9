
#include <stdio.h>

#include <filter2d.h>

#include <math.h>

#include <stdlib.h>


void filter2d_mean3(double * img_in, long width, long height, double * img_out)
{
	int x, y, k, l,i;
	double val = 1.0 / 9.0;
	double s;
	double H[3][3] = { {val, val, val},
						{val,val,val},
						{val,val,val} };

	//Recopie des bords
	for (i = 0; i < width * height; i++) img_out[i] = img_in[i];	//Pas de soucis à l'écrire comme ça, on perd pas de temps d'execution

	//Convolution sur le "coeur"
	for (x = 1; x < width - 1; x++) 
		for (y = 1; y < height - 1; y++)
		{
			s = 0.0;
			for (k=-1;k<=1;k++)		//Fait apparaitre la symétrie du filtre, homogene avec ce qu'on voit d'habitude
				for (l=-1;l<=1;l++)	//Fait apparaitre la dimension du filtre
					s += img_in[(y - l) * width + (x - k)] * H[l+1][k+1];	
			// Alternative
			/*
			for (k = 0; k < 3; k++)		//Fait apparaitre la symétrie du filtre, homogene avec ce qu'on voit d'habitude
				for (l = 0; l < 3; l++)	//Fait apparaitre la dimension du filtre
					s += img_in[(y - (l-1)) * width + (x - (k-1)] * H[l + 1][k + 1];
			*/
			
			img_out[y * width + x] = s;

		}

}


void sobel_norm(double * img_in, long width, long height, double * img_out)
{
	int x, y, k, l, i;
	double fact = 1.0 / 8.0;
	double sx,sy;
	double Sx[3][3] = { {fact * 1.0, 0.0, fact * -1.0},
						{fact * 2.0,0.0, fact * -2.0 },
						{fact * 1.0,0.0, fact * -1.0} };
	/*double Sy[3][3] = {{fact * 1.0, fact * 2.0, fact * 1.0},
						{0.0, 0.0 ,0.0},
						{fact * -1.0,fact * -2.0,fact * -1.0} };
	*/
	//On peut utilser directemznt transpose
	//Recopie des bords
	for (i = 0; i < width * height; i++) img_out[i] = img_in[i];	//Pas de soucis à l'écrire comme ça, on perd pas de temps d'execution

	//Convolution sur le "coeur"
	for (x = 1; x < width - 1; x++)
		for (y = 1; y < height - 1; y++)
		{
			sx = sy = 0.0;
			for (k = -1; k <= 1; k++)		//Fait apparaitre la symétrie du filtre, homogene avec ce qu'on voit d'habitude
				for (l = -1; l <= 1; l++)	//Fait apparaitre la dimension du filtre
				{
					sx += img_in[(y - l) * width + (x - k)] * Sx[l + 1][k + 1];
					sy += img_in[(y - l) * width + (x - k)] * Sx[k + 1][l + 1];
				}
			// Alternative
			/*
			for (k = 0; k < 3; k++)		//Fait apparaitre la symétrie du filtre, homogene avec ce qu'on voit d'habitude
				for (l = 0; l < 3; l++)	//Fait apparaitre la dimension du filtre
					s += img_in[(y - (l-1)) * width + (x - (k-1)] * H[l + 1][k + 1];
			*/

			img_out[y * width + x] = sqrt(sx * sx + sy * sy);
		}
}


static double * gaussian2d_create(double sigma, long * size)
{

	// Calcul de la taille du filtre
	/*int M = 1 + 2 * (3 * sigma);
	int x, y;
	int val_fin = (floor(3 * sigma))/2 + 1;
	double H[M][M];
	for (x = - val_fin ; x<val_fin; x++)
		for (y = - val_fin; y<val_fin; y++)
	*/

	double * filter = NULL, s = 0.0,v;
	int t,i,j,c;
	t = 1 + 2 * (int) ceil(3.0 * sigma);
	c = t / 2;	// Division entière car c est un int -> c est inutile au final on pourrait mettre directement t/2 dans les calculs
	*size = t;
	filter = (double*)calloc(t * t, sizeof(double));
	// Première méthode
	/*
	double x,y;
	for (i = 0; i < t; i++)
		x = double (i-t/2)/sigma;
		for (j = 0; j < t; j++)
			{
				y = double (j-t/2)/sigma;
				v = exp(-(double)(x*x + y*y) / 2.0 * sigma * sigma);
				s += v;
				filter[j * t + i] = v;
			}
	*/
	// Deuxième version 
	for (i=-t/2;i<=t/2;i++)
		for (j = -t/2; j <= t/2; j++) 
			{
				v = exp(-(double)(i*i + j*j) / (2.0 * sigma * sigma));
				s += v;
				filter[(j + t / 2) * t + i + t / 2] = v;
			}

	for (i = 0; i < t; i++)
		for (j = 0; j < t; j++)
			filter[j * t + i] /= s;
	return filter;
}


static void convolution2d(double * img_in, long width, long height,
									        double * mask2d, long tx, long ty,
									        double * img_out)
{
	int x, y, k, l, i;
	double s;

	//Problème des bords
	for (i = 0; i < width * height; i++) img_out[i] = 0.0;	//Pas de soucis à l'écrire comme ça, on perd pas de temps d'execution

	//Convolution sur le "coeur"
	for (x = tx/2; x < width - tx / 2; x++)
		for (y = ty / 2; y < height - ty / 2; y++)
		{
			s = 0.0;
			for (k = -tx / 2; k <= tx / 2; k++)
				for (l = -ty / 2; l <= ty / 2; l++)
					s += img_in[(y - l) * width + (x - k)] * mask2d[(l + ty / 2) * tx + (k + tx / 2)];
			// Alternative
			/*
			for (k = 0; k < tx; k++)
				for (l = 0; l < ty; l++)
					s += img_in[(y - (l-(ty/2))) * width + (x - (k-tx/2)] * mask2d[l * tx + k];
			*/

			img_out[y * width + x] = s;
		}
}

static double ** gaussian2d_create_matrix(double sigma, long * size)
{

	double** filter = NULL, s = 0.0, v;
	int t, i, j, c;
	t = 1 + 2 * (int)ceil(3.0 * sigma);
	c = t / 2;	// Division entière car c est un int -> c est inutile au final on pourrait mettre directement t/2 dans les calculs
	*size = t;

	// Première allocation
	/*filter = (double **)calloc(t, sizeof(double *));
	for (j = 0; j < t ; j++) filter[j]= (double *)calloc(t, sizeof(double));*/

	// Deuxième allocation en mieux 
	filter = (double**)calloc(t, sizeof(double*));
	filter[0]=(double*)calloc(t*t, sizeof(double));
	for (j = 1; j < t; j++) filter[j] = filter[0]+j*t; // filter[j] = filter[j-1]+t

	for (i = -t / 2; i <= t / 2; i++)
		for (j = -t / 2; j <= t / 2; j++)
		{
			v = exp(-(double)(i * i + j * j) / (2.0 * sigma * sigma));
			s += v;
			filter[j + t / 2] [i + t / 2] = v;
		}

	for (i = 0; i < t; i++)
		for (j = 0; j < t; j++)
			filter[j][i] /= s;
	return filter;
}


static void convolution2d_by_matrix(double * img_in, long width, long height,
																		double ** mask2d, long tx, long ty,
									                  double * img_out)
{
	int x, y, k, l, i;
	double s;

	//Problème des bords
	for (i = 0; i < width * height; i++) img_out[i] = 0.0;	//Pas de soucis à l'écrire comme ça, on perd pas de temps d'execution

	//Convolution sur le "coeur"
	for (x = tx / 2; x < width - tx / 2; x++)
		for (y = ty / 2; y < height - ty / 2; y++)
		{
			s = 0.0;
			for (k = -tx / 2; k <= tx / 2; k++)
				for (l = -ty / 2; l <= ty / 2; l++)
					s += img_in[(y - l) * width + (x - k)] * mask2d[l + ty / 2] [k + tx / 2];
			img_out[y * width + x] = s;
		}
}


void filter2d_gaussian(double * img_in, long width, long height,
											 double sigma,
											 double * img_out)
{	
	//// cas vectoriel
	//double* filter2d = NULL;
	//long size;
	//filter2d = gaussian2d_create(sigma, &size);
	//convolution2d(img_in, width, height, filter2d, size, size, img_out);
	//free(filter2d);
	//// après un free on peut faire ça
	////filter2d = NULL;


	// cas matriciel
	double ** filter2d = NULL;
	long size, j;
	filter2d = gaussian2d_create_matrix(sigma, &size);

	// sans la meilleure allocation
	/*convolution2d_by_matrix(img_in, width, height, filter2d, size, size, img_out);
	for (j = 0; j < size; j++) free(filter2d[j]) ;*/

	// avec la meilleure allocation
	convolution2d_by_matrix(img_in, width, height, filter2d[0], size, size, img_out);
	free(filter2d[0]);


	free(filter2d);

}


void img_get_raw(double * img, long width, long height,  // lecture
		             long no,
						     double * v)
{
	int i;

	for (i = 0; i < width; i++) v[i] = img[no * width + i];
}


void img_set_raw(double * img, long width, long height,  //écriture
								 long no,
								 double * v)
{
	int i;

	for (i = 0; i < width; i++) img[no * width + i] = v[i];
}


void img_get_column(double * img, long width, long height, 
									  long no,
										double * v)
{
	int j;

	for (j = 0; j < height; j++) v[j] = img[no + j*width];
}


void img_set_column(double * img, long width, long height,
										long no,
										double * v)
{
	int j;

	for (j = 0; j < height; j++) img[no  + j*width] = v[j] ;
}




static double * gaussian1d_create(double sigma, long * size)
{
	double* filter = NULL, s = 0.0, v;
	int t, i;

	t = 1 + 2 * (int)ceil(3.0 * sigma);
	*size = t;
	filter = (double*)calloc(t, sizeof(double));
	
	for (i = -t / 2; i <= t / 2; i++)
	{
		v = exp(-(double)(i*i) / (2.0 * sigma * sigma));
		s += v;
		filter[ i + t / 2] = v;
	}

	for (i = 0; i < t; i++) filter[i] /= s;
	return filter;
}


static void convolution1d(double * v_in, long size,
													double * mask1d, long t,
													double * v_out)
{
	int x, k, i;
	double s;

	//Problème des bords
	for (i = 0; i < size; i++) v_out[i] = 0.0;

	//Convolution sur le "coeur"
	for (x = t/ 2; x < size - t/ 2; x++)
	{
		s = 0.0;
		for (k = -t/ 2; k <= t/ 2; k++)
			s += v_in[x - k] * mask1d[k + t/2];
		v_out[x] = s;
	}
}

// écriture d'une MACRO en c 
#define MAX(a,b) (a > b ? a : b) // si a > b alors a, sinon b

// #define MAX(a,b) ((a)+(b)+fabs((a)-(b)))/2.0) //-> plus lent car calculs et pas juste une comparaison 


void filter2d_gaussian_fast(double * img_in, long width, long height,
														double sigma,
														double * img_out)
{	
	int t,i,j;
	double* filter = NULL, *v_in = NULL, * v_out = NULL;

	// Création de la gaussienne 1d
	filter = gaussian1d_create(sigma, &t);

	// allocation v_in et v_out
	v_in = (double*)calloc(MAX(width,height), sizeof(double)); // MAX(width,height) pour l'utiliser pour lignes et colonne
	v_out = (double*)calloc(MAX(width, height), sizeof(double));

	// parcours en lignes
	for (j = 0; j < height; j++) {
		// - lecture lignes : img_in ->v_in
		img_get_raw(img_in, width, height, j, v_in);
		// - conv 1d
		convolution1d(v_in, width, filter, t, v_out);
		// - écriture : v_out-> img_out
		img_set_raw(img_out, width, height, j, v_out);
	}

	// parcours en colonne
	for (i = 0; i < width; i++) {
		// - lecture colonne : img_in ->v_in
		img_get_column(img_out, width, height, i, v_in);
		// - conv 1d
		convolution1d(v_in, height, filter, t, v_out);
		// - écriture : v_out-> img_out
		img_set_column(img_out, width, height, i, v_out);
	}

	// libération v_in et v_out
	free(v_in);
	free(v_out);

	// libération de la gaussienne 1d
	free(filter);
}




typedef double (* PROC)(double *, long);


static int compare(void const* a, void const* b)
{
	double* pa = (double*)a, * pb = (double*)b;
	if (*pa > *pb) return 1;
	if (*pa < *pb) return -1;
	return 0;
}

static double v_min(double* v, long size) 
{	
	double min = v[0];
	int i;
	for (i = 1; i < size; i++) {
		if (min > v[i]) min = v[i];
	}
	return min;
}

static double v_max(double* v, long size)
{
	double max = v[0];
	int i;
	for (i = 1; i < size; i++) {
		if (max < v[i]) max = v[i];
	}
	return max;
}


static double v_mean(double* v, long size)
{
	double somme = 0.0;
	int i;
	for (i = 0; i < size; i++) somme += v[i];
	return somme/(double )size;
}

static double v_median(double* v, long size)
{
	qsort((void*)v, (size_t)(size), sizeof(double), compare);
	if (size % 2 ==0) return (v[size/2]+ v[size/2+1])/2;
	else return v[(size+1)/2];
}



static void filter2d_generic(double * img_in, long width, long height,
							  						 long tx, long ty, PROC proc,
														 double * img_out)
{
	int x, y, k, l, i, idx;
	double* v = NULL;

	v = (double*)calloc(tx * ty, sizeof(double));

	//Problème des bords
	for (i = 0; i < width * height; i++) img_out[i] = img_in[i];	//Pas de soucis à l'écrire comme ça, on perd pas de temps d'execution

	//Convolution sur le "coeur"
	for (x = tx / 2; x < width - tx / 2; x++)
		for (y = ty / 2; y < height - ty / 2; y++)
		{
			idx = 0;
			for (k = -tx / 2; k <= tx / 2; k++)
				for (l = -ty / 2; l <= ty / 2; l++)
					v[idx++] = img_in[(y + l) * width + x - k];

			img_out[y * width + x] = proc(v,tx*ty );
		}
	free(v);
}


void filter2d_by_method(double* img_in, long width, long height,
	long tx, long ty, long method,
	double* img_out)
{	
	/*if (method == 0)
		filter2d_generic(img_in, width, height, tx, ty, v_min, img_out);
	if (method == 1)
		filter2d_generic(img_in, width, height, tx, ty, v_max, img_out);
	if (method == 2)
		filter2d_generic(img_in, width, height, tx, ty, v_mean, img_out);
	if (method == 3)
		filter2d_generic(img_in, width, height, tx, ty, v_median, img_out);*/

	PROC v_proc[4] = { v_min, v_max, v_mean, v_median };

	filter2d_generic(img_in, width, height, tx, ty, v_proc[method], img_out);
}




