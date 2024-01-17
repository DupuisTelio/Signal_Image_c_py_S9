
#include <stdlib.h>
#include <stdio.h>



#include <median2d.h>



static int compare(void const * a, void const * b)
{
	double* pa = (double*)a, * pb = (double*)b;
	if (*pa > *pb) return 1; 
	if (*pa < *pb) return -1;
	return 0;
}


void filter2d_median(double * img_in, long width, long height,
										 long tx, long ty,
										 double * img_out)
{

		int x, y, k, l, i,idx;
		double * v =NULL;

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
						v[idx++]= img_in[(y +l)*width + x - k];

				qsort((void*)v, (size_t)(tx * ty), sizeof(double), compare);

				img_out[y * width + x] = v[tx*ty/2];
			}
	
}