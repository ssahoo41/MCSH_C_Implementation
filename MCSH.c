# include <stdio.h>
# include <stdlib.h>
# include <string.h>
#include <assert.h>
#include <time.h> 
#include <mpi.h>
# include <math.h>
/* BLAS, LAPACK, LAPACKE routines */
//#ifdef USE_MKL
    // #define MKL_Complex16 double complex
    //#include <mkl.h>
//#else
//    #include <cblas.h>
//#endif

#include "MCSHHelper.h"
#include "MCSH.h"


void MaxwellCartesianSphericalHarmonics(const double *x, const double *y, const double *z, const int l, const char *n, const double rCutoff, double *result, const int size)
{
	double *r = calloc( size, sizeof(double));
	double *x_hat = calloc( size, sizeof(double));
	double *y_hat = calloc( size, sizeof(double));
	double *z_hat = calloc( size, sizeof(double));


	getRArray(x, y, z, r, size);
	divideVector(x, r, x_hat, size);
	divideVector(y, r, y_hat, size);
	divideVector(z, r, z_hat, size);

	// double *uncutResult = malloc( size * sizeof(double));
	// double *result = malloc( size * sizeof(double));

	int i;
	// printf("\n============\n");
	// for (i = 0; i < size; i++)
	// {
	// 	printf("x: %10f \t y: %10f \t z: %10f \t x_hat: %10f \t y_hat: %10f \t z_hat: %10f \t r: %10f\n", x[i],y[i],z[i],x_hat[i],y_hat[i],z_hat[i], r[i]);
	// }

	// for (i = 0; i < size; i++)
	// {
	// 	double r_calc = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
	// 	printf("x: %10f \t y: %10f \t z: %10f \t r: %10f \t r_calc: %././././././10f \t rCutoff: %10f \t r>rCutoff%d\n", x[i],y[i],z[i], r[i],  r_calc, rCutoff, r[i]>rCutoff);
	// }
	switch (l) 
	{
		case 0:
			//uncutResult = 1;
			// printf(" l = 0\n");

			if (strcmp(n, "000") == 0) 
			{
				for ( i = 0; i < size; i++)
				{
					result[i] = 1.0;

				}
			} 
			else
			{
				printf("\nWARNING: n is not valid %s \n", n);
			}
			break;

		case 1:
			// printf(" l = 1 \t");
			if (strcmp(n, "100") == 0) 
			{
				// printf(" n = 100 \n");
				for ( i = 0; i < size; i++)
				{
					result[i] = x_hat[i];

				}
			} 
			else if (strcmp(n, "010") == 0)
			{
				for ( i = 0; i < size; i++)
				{
					result[i] = y_hat[i];
				}
			}
			else if (strcmp(n, "001") == 0)
			{
				for ( i = 0; i < size; i++)
				{
					result[i] = z_hat[i];
				}
			}
			else
			{
				printf("\nWARNING: n is not valid %s \n", n);
			}
			break;

		case 2:
			if (strcmp(n, "200") == 0) 
			{
				// result = 3.0 * x_hat * x_hat - 1.0;
				double *temp = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 0, 3.0, temp, size);
				addScalarVector(temp, -1.0, result, size);
				free(temp);
			} 
			else if (strcmp(n, "020") == 0)
			{
				// result = 3.0 * y_hat * y_hat - 1.0;
				double *temp = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 0, 3.0, temp, size);
				addScalarVector(temp, -1.0, result, size);
				free(temp);
			}
			else if (strcmp(n, "002") == 0)
			{
				// result = 3.0 * z_hat * z_hat - 1.0;
				double *temp = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 2, 3.0, temp, size);
				addScalarVector(temp, -1.0, result, size);
				free(temp);
			}
			else if (strcmp(n, "110") == 0)
			{
				// result = 3.0 * x_hat * y_hat;
				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 0, 3.0, result, size);
			}
			else if (strcmp(n, "101") == 0)
			{
				// result = 3.0 * x_hat * z_hat;
				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 1, 3.0, result, size);
			}
			else if (strcmp(n, "011") == 0)
			{
				// result = 3.0 * y_hat * z_hat;
				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 1, 3.0, result, size);
			}
			else
			{
				printf("\nWARNING: n is not valid %s \n", n);
			}
			break;

		case 3:
			if (strcmp(n, "300") == 0) 
			{
				// result = 15.0 * x_hat * x_hat * x_hat - 9.0 * x_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 3, 0, 0, 15.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 0, 9.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			} 
			else if (strcmp(n, "030") == 0)
			{
				// result = 15.0 * y_hat * y_hat * y_hat - 9.0 * y_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 3, 0, 15.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 0, 9.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "003") == 0)
			{
				// result = 15.0 * z_hat * z_hat * z_hat - 9.0 * z_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 3, 15.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 1, 9.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "210") == 0)
			{
				// result = 15.0 * x_hat * x_hat * y_hat - 3.0 * y_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 2, 1, 0, 15.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 0, 3.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "120") == 0)
			{
				// result = 15.0 * x_hat * y_hat * y_hat - 3.0 * x_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 1, 2, 0, 15.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 0, 3.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "201") == 0)
			{
				// result = 15.0 * x_hat * x_hat * z_hat - 3.0 * z_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 1, 15.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 1, 3.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "102") == 0)
			{
				// result = 15.0 * x_hat * z_hat * z_hat - 3.0 * x_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 2, 15.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 0, 3.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "021") == 0)
			{
				// result = 15.0 * y_hat * y_hat * z_hat - 3.0 * z_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 1, 15.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 1, 3.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "012") == 0)
			{
				// result = 15.0 * y_hat * z_hat * z_hat - 3.0 * y_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 2, 15.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 0, 3.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "111") == 0)
			{
				// result = 15.0 * x_hat * y_hat * z_hat;
				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 1, 15.0, result, size);
			}
			else
			{
				printf("\nWARNING: n is not valid %s \n", n);
			}
			break;

		case 4:
			if (strcmp(n, "400") == 0) 
			{
				// result = 105.0 * x_hat * x_hat * x_hat * x_hat - 90.0 * x_hat * x_hat + 9.0;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				double *temp3 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 4, 0, 0, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 0, 90.0, temp2, size);
				subtractVector(temp1, temp2, temp3, size);
				addScalarVector(temp3, 9.0, result, size);
				free(temp1);
				free(temp2);
				free(temp3);
			} 
			else if (strcmp(n, "040") == 0)
			{
				// result = 105.0 * y_hat * y_hat * y_hat * y_hat - 90.0 * y_hat * y_hat + 9.0;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				double *temp3 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 4, 0, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 0, 90.0, temp2, size);
				subtractVector(temp1, temp2, temp3, size);
				addScalarVector(temp3, 9.0, result, size);
				free(temp1);
				free(temp2);
				free(temp3);
			}
			else if (strcmp(n, "004") == 0)
			{
				// result = 105.0 * z_hat * z_hat * z_hat * z_hat - 90.0 * z_hat * z_hat + 9.0;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				double *temp3 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 4, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 2, 90.0, temp2, size);
				subtractVector(temp1, temp2, temp3, size);
				addScalarVector(temp3, 9.0, result, size);
				free(temp1);
				free(temp2);
				free(temp3);
			}
			else if (strcmp(n, "310") == 0)
			{
				// result = 105.0 * x_hat * x_hat * x_hat * y_hat - 45.0 * x_hat * y_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 3, 1, 0, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 0, 45.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "130") == 0)
			{
				// result = 105.0 * x_hat * y_hat * y_hat * y_hat - 45.0 * x_hat * y_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 1, 3, 0, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 0, 45.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "301") == 0)
			{
				// result = 105.0 * x_hat * x_hat * x_hat * z_hat - 45.0 * x_hat * z_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 3, 0, 1, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 1, 45.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "103") == 0)
			{
				// result = 105.0 * x_hat * z_hat * z_hat * z_hat - 45.0 * x_hat * z_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 3, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 1, 45.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "031") == 0)
			{
				// result = 105.0 * y_hat * y_hat * y_hat * z_hat - 45.0 * y_hat * z_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 3, 1, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 1, 45.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "013") == 0)
			{
				// result = 105.0 * y_hat * z_hat * z_hat * z_hat - 45.0 * y_hat * z_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 3, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 1, 45.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "220") == 0)
			{
				// result = 105.0 * x_hat * x_hat * y_hat * y_hat - 15.0 * x_hat * x_hat - 15.0 * y_hat * y_hat + 3.0;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				double *temp3 = malloc( size * sizeof(double));
				double *temp4 = malloc( size * sizeof(double));
				double *temp5 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 2, 2, 0, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 0, 15.0, temp2, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 0, 15.0, temp3, size);
				subtractVector(temp1, temp2, temp4, size);
				subtractVector(temp4, temp3, temp5, size);
				addScalarVector(temp5, 3.0, result, size);
				free(temp1);
				free(temp2);
				free(temp3);
				free(temp4);
				free(temp5);
			}
			else if (strcmp(n, "202") == 0)
			{
				// result = 105.0 * x_hat * x_hat * z_hat * z_hat - 15.0 * x_hat * x_hat - 15.0 * z_hat * z_hat + 3.0;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				double *temp3 = malloc( size * sizeof(double));
				double *temp4 = malloc( size * sizeof(double));
				double *temp5 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 2, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 0, 15.0, temp2, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 2, 15.0, temp3, size);
				subtractVector(temp1, temp2, temp4, size);
				subtractVector(temp4, temp3, temp5, size);
				addScalarVector(temp5, 3.0, result, size);
				free(temp1);
				free(temp2);
				free(temp3);
				free(temp4);
				free(temp5);
			}
			else if (strcmp(n, "022") == 0)
			{
				// result = 105.0 * y_hat * y_hat * z_hat * z_hat - 15.0 * y_hat * y_hat - 15.0 * z_hat * z_hat + 3.0;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				double *temp3 = malloc( size * sizeof(double));
				double *temp4 = malloc( size * sizeof(double));
				double *temp5 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 2, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 0, 15.0, temp2, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 2, 15.0, temp3, size);
				subtractVector(temp1, temp2, temp4, size);
				subtractVector(temp4, temp3, temp5, size);
				addScalarVector(temp5, 3.0, result, size);
				free(temp1);
				free(temp2);
				free(temp3);
				free(temp4);
				free(temp5);
			}
			else if (strcmp(n, "211") == 0)
			{
				// result = 105.0 * x_hat * x_hat * y_hat * z_hat - 15.0 * y_hat * z_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 2, 1, 1, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 1, 15.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "121") == 0)
			{
				// result = 105.0 * x_hat * y_hat * y_hat * z_hat - 15.0 * x_hat * z_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 1, 2, 1, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 1, 15.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else if (strcmp(n, "112") == 0)
			{
				// result = 105.0 * x_hat * y_hat * z_hat * z_hat - 15.0 * x_hat * y_hat;
				double *temp1 = malloc( size * sizeof(double));
				double *temp2 = malloc( size * sizeof(double));
				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 2, 105.0, temp1, size);
				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 0, 15.0, temp2, size);
				subtractVector(temp1, temp2, result, size);
				free(temp1);
				free(temp2);
			}
			else
			{
				printf("\nWARNING: n is not valid %s \n", n);
			}
			break;

		default:
			printf("\nWARNING: l is not valid %d \n", l);
			break;
	}

	//int i;
	for (i = 0; i < size; i++)
	{	
		//printf("before: %10f \t",result[i]);
		if (r[i] > rCutoff)
		{
			result[i] = 0.0;
		}
		//printf("%10f \t %10f \n", uncutResult[i], result[i]);
		//printf("after: %10f \n",result[i]);
	}

	// for (i = 0; i < size; i++)
	// {
	// 	printf("x: %10f \t y: %10f \t z: %10f \t r: %10f \t rCutoff: %f \t l: %d \t n: %s \t result: %10f\n", x[i],y[i],z[i], r[i], rCutoff, l, n, result[i]);
	// }


	//free(uncutResult);
	free(r);
	free(x_hat);
	free(y_hat);
	free(z_hat);

	// return result;
}


void calculateStencil(const int stencilDimX, const int stencilDimY, const int stencilDimZ, const double hx, const double hy, const double hz, 
					  const double rCutoff, const int l, const char *n, const int radialFunctionType, const int radialFunctionOrder, 
					  const double *U, const int accuracy, double *stencil)
{
	int pixelEvalArrSize = accuracy * accuracy * accuracy;
	printf("%d", pixelEvalArrSize);
	double dv = calcDv(hx, hy, hz, accuracy,U);
	printf("%f", dv);
	double *refX = calloc( pixelEvalArrSize, sizeof(double));
	double *refY = calloc( pixelEvalArrSize, sizeof(double));
	double *refZ = calloc( pixelEvalArrSize, sizeof(double));

	getCentralCoords(hx, hy, hz, accuracy, refX, refY, refZ);

	//char stencilRefXFilename[128];
	//snprintf(stencilRefXFilename, 128, "stencil_RefX_%f_%d_%s_Type%d_%d.csv", rCutoff, l, n, radialFunctionType, radialFunctionOrder);
	//writeMatToFile(stencilRefXFilename, refX, accuracy, accuracy, accuracy);

	//char stencilRefYFilename[128];
	//snprintf(stencilRefYFilename, 128, "stencil_RefY_%f_%d_%s_Type%d_%d.csv", rCutoff, l, n, radialFunctionType, radialFunctionOrder);
	//writeMatToFile(stencilRefYFilename, refY, accuracy, accuracy, accuracy);

	//char stencilRefZFilename[128];
	//snprintf(stencilRefZFilename, 128, "stencil_RefZ_%f_%d_%s_Type%d_%d.csv", rCutoff, l, n, radialFunctionType, radialFunctionOrder);
	//writeMatToFile(stencilRefZFilename, refZ, accuracy, accuracy, accuracy);

	int centerX = (stencilDimX - 1)/2;
    	int centerY = (stencilDimY - 1)/2;
    	int centerZ = (stencilDimZ - 1)/2;

	double *tempXArr = calloc( pixelEvalArrSize, sizeof(double));
	double *tempYArr = calloc( pixelEvalArrSize, sizeof(double));
	double *tempZArr = calloc( pixelEvalArrSize, sizeof(double));
	double *tempMCSHResult = calloc( pixelEvalArrSize, sizeof(double));
	double *tempRadialResult = calloc( pixelEvalArrSize, sizeof(double));
	double xOffset, yOffset, zOffset;
	int i, j, k, index = 0, m;
	for (k = 0; k < stencilDimZ; k++){
		for ( j = 0; j < stencilDimY; j++) {
			for ( i = 0; i < stencilDimX; i++) {
				//printf("start %d %d %d\n", i,j,k);
				xOffset = (i-centerX) * hx;
				yOffset = (j-centerY) * hy;
				zOffset = (k-centerZ) * hz;
				//index = k * stencilDimX * stencilDimY + j * stencilDimX + i;

				addScalarVector(refX, xOffset, tempXArr, pixelEvalArrSize);
				addScalarVector(refY, yOffset, tempYArr, pixelEvalArrSize);
				addScalarVector(refZ, zOffset, tempZArr, pixelEvalArrSize);

				applyU2(tempXArr, tempYArr, tempZArr, U, pixelEvalArrSize);
				MaxwellCartesianSphericalHarmonics(tempXArr, tempYArr, tempZArr, l, n, rCutoff, tempMCSHResult, pixelEvalArrSize);

				if (radialFunctionType == 2)
				{
					LegendrePolynomial(tempXArr, tempYArr, tempZArr, radialFunctionOrder, rCutoff, tempRadialResult, pixelEvalArrSize);
					multiplyVector(tempMCSHResult, tempRadialResult, tempMCSHResult, pixelEvalArrSize);
				}
				//stencil[index] = cblas_dasum(pixelEvalArrSize, tempMCSHResult, 1) * dv;
				//printf("start %d %d %d \t sum: %f \t sumAbs: %f \t dv: %f\n", i,j,k, sumArr(tempMCSHResult, pixelEvalArrSize), sumAbsArr(tempMCSHResult, pixelEvalArrSize), dv);
				stencil[index] = sumArr(tempMCSHResult, pixelEvalArrSize) * dv;
				index++;
			}
		}
	}
	free(refX);
	free(refY);
	free(refZ);
	free(tempXArr);
	free(tempYArr);
	free(tempZArr);
	free(tempMCSHResult);
}


void calcStencilAndConvolveAndAddResult(const double *image, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
										const double rCutoff, const int l, const char *n, const int radialFunctionType, const int radialFunctionOrder, 
										const double *U, const int accuracy, double *convolveResult, const int stencilIndex)
{
	double start_t, end_stencil_t, end_convolve_t, end_convolve2_t; 
	//time(&start_t); 
	start_t = MPI_Wtime();

	int pixelEvalArrSize = accuracy * accuracy * accuracy;

	int stencilDimX, stencilDimY, stencilDimZ;
	GetDimensionsPlane(hx, hy, hz, rCutoff, U, &stencilDimX, &stencilDimY, &stencilDimZ);

	double *stencil = calloc( stencilDimX * stencilDimY * stencilDimZ, sizeof(double));
	calculateStencil(stencilDimX, stencilDimY, stencilDimZ, hx, hy, hz, 
					 rCutoff, l, n, radialFunctionType, radialFunctionOrder, U, accuracy, stencil);


	//char stencilFilename[128];
	//snprintf(stencilFilename, 128, "stencil_%f_%d_%s_Type%d_%d.csv", rCutoff, l, n, radialFunctionType, radialFunctionOrder);
	//writeMatToFile(stencilFilename, stencil, stencilDimX, stencilDimY, stencilDimZ);

	printf("end stencil calculation \n");
	//time(&end_stencil_t); 
	end_stencil_t  = MPI_Wtime();


	convolve5(image, stencil, imageDimX, imageDimY, imageDimZ, stencilDimX, stencilDimY, stencilDimZ, convolveResult);
	char convolveResultFilename[128];
	snprintf(convolveResultFilename, 128, "%s_%d.csv", "convolve_result", stencilIndex);

	writeMatToFile(convolveResultFilename, convolveResult, imageDimX, imageDimY, imageDimZ);	
	// printf("end convolve1 \n");
	// time(&end_convolve_t);
	end_convolve_t  = MPI_Wtime();

	printf("\n r: %f \t l: %d \t n: %s \t total_time: %f \t stencil: %f \t convolve: %f \n",rCutoff, l, n,end_convolve_t - start_t, end_stencil_t - start_t, end_convolve_t - end_stencil_t);

	free(stencil);
}


void calcStencilAndConvolveAndSave(const double *image, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
								   const double rCutoff, const int l, const char *n, const int radialFunctionType, const int radialFunctionOrder, 
								   const double *U, const int accuracy, const int stencilIndex)
{
	int pixelEvalArrSize = accuracy * accuracy * accuracy;

	// double dv = calcDv(hx, hy, hz, accuracy,U);

	int stencilDimX, stencilDimY, stencilDimZ;
	GetDimensionsPlane(hx, hy, hz, rCutoff, U, &stencilDimX, &stencilDimY, &stencilDimZ);

	double *stencil = calloc( stencilDimX * stencilDimY * stencilDimZ, sizeof(double));
	calculateStencil(stencilDimX, stencilDimY, stencilDimZ, hx, hy, hz, 
					 rCutoff, l, n, radialFunctionType, radialFunctionOrder, U, accuracy, stencil);

	

	char stencilFilename[128];

	snprintf(stencilFilename, 128, "%s_%d.csv", "stencil", stencilIndex);
	printf("%s",stencilFilename);

	writeMatToFile(stencilFilename, stencil, stencilDimX, stencilDimY, stencilDimZ);


	int imageSize = imageDimX * imageDimY * imageDimZ;
	double *convolveResult = calloc( imageDimX * imageDimY * imageDimZ, sizeof(double));


	convolve3(image, stencil, imageDimX, imageDimY, imageDimZ, stencilDimX, stencilDimY, stencilDimZ, convolveResult);

	char convolveResultFilename[128];
	snprintf(convolveResultFilename, 128, "%s_%d.csv", "convolve_result", stencilIndex);

	writeMatToFile(convolveResultFilename, convolveResult, imageDimX, imageDimY, imageDimZ);

	free(convolveResult);
}


void calcStencilAndSave(const double hx, const double hy, const double hz, const double rCutoff, const int l, const char *n, const double *U, 
						const int radialFunctionType, const int radialFunctionOrder, const int accuracy, const int stencilIndex)
{	
	int pixelEvalArrSize = accuracy * accuracy * accuracy;

	// double dv = calcDv(hx, hy, hz, accuracy,U);

	int stencilDimX, stencilDimY, stencilDimZ;
	GetDimensionsPlane(hx, hy, hz, rCutoff, U, &stencilDimX, &stencilDimY, &stencilDimZ);

	double *stencil = calloc( stencilDimX * stencilDimY * stencilDimZ, sizeof(double));
	calculateStencil(stencilDimX, stencilDimY, stencilDimZ, hx, hy, hz, 
					 rCutoff, l, n, radialFunctionType, radialFunctionOrder, U, accuracy, stencil);


	// printf("after stencil calculation");
	char stencilFilename[128];
	snprintf(stencilFilename, 128, "%s_%d.csv", "test_output", stencilIndex);
	writeMatToFile(stencilFilename, stencil, stencilDimX, stencilDimY, stencilDimZ);

	free(stencil);
}







void prepareMCSHFeatureAndSave(const double *image, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
							   const double rCutoff, const int l, const int radialFunctionType, const int radialFunctionOrder, const double *U, const int accuracy)
{	
	int imageSize = imageDimX * imageDimY * imageDimZ;
	double *featureVector = calloc( imageSize, sizeof(double));

	int i, j, k;
	switch (l) 
	{
		case 0:
			calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "000", radialFunctionType, radialFunctionOrder, U, accuracy, featureVector, 0);
			break;

		case 1:
			{
				double *component1 = calloc( imageSize, sizeof(double));
				double *component2 = calloc( imageSize, sizeof(double));
				double *component3 = calloc( imageSize, sizeof(double));

				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "100", radialFunctionType, radialFunctionOrder, U, accuracy, component1, 10);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "010", radialFunctionType, radialFunctionOrder, U, accuracy, component2, 11);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "001", radialFunctionType, radialFunctionOrder, U, accuracy, component3, 12);
	
				powVector(component1, 2, component1, imageSize);
				powVector(component2, 2, component2, imageSize);
				powVector(component3, 2, component3, imageSize);

				addVector(component1, component2, featureVector, imageSize);
				addVector(featureVector, component3, featureVector, imageSize);

				sqrtVector(featureVector, featureVector, imageSize);

				free(component1);
				free(component2);
				free(component3);
			}
			break;

		case 2:
			{
				double *component1 = calloc( imageSize, sizeof(double));
				double *component2 = calloc( imageSize, sizeof(double));
				double *component3 = calloc( imageSize, sizeof(double));
			
				double *component4 = calloc( imageSize, sizeof(double));
				double *component5 = calloc( imageSize, sizeof(double));
				double *component6 = calloc( imageSize, sizeof(double));

				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "200", radialFunctionType, radialFunctionOrder, U, accuracy, component1, 21);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "020", radialFunctionType, radialFunctionOrder, U, accuracy, component2, 22);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "002", radialFunctionType, radialFunctionOrder, U, accuracy, component3, 23);
			
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "110", radialFunctionType, radialFunctionOrder, U, accuracy, component4, 24);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "101", radialFunctionType, radialFunctionOrder, U, accuracy, component5, 25);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "011", radialFunctionType, radialFunctionOrder, U, accuracy, component6, 26);

				powVector(component1, 2, component1, imageSize);
				powVector(component2, 2, component2, imageSize);
				powVector(component3, 2, component3, imageSize);
			
				powVector(component4, 2, component4, imageSize);
				powVector(component5, 2, component5, imageSize);
				powVector(component6, 2, component6, imageSize);	
			
				multiplyScalarVector(component4, 2.0, component4, imageSize);
				multiplyScalarVector(component5, 2.0, component5, imageSize);
				multiplyScalarVector(component6, 2.0, component6, imageSize);

				addVector(component1, component2, featureVector, imageSize);
				addVector(featureVector, component3, featureVector, imageSize);	
				addVector(featureVector, component4, featureVector, imageSize);	
				addVector(featureVector, component5, featureVector, imageSize);	
				addVector(featureVector, component6, featureVector, imageSize);	

				sqrtVector(featureVector, featureVector, imageSize);

				free(component1);
				free(component2);
				free(component3);
				free(component4);
				free(component5);
				free(component6);
			}
			break;

		case 3:
			{
				double *component1 = calloc( imageSize, sizeof(double));
				double *component2 = calloc( imageSize, sizeof(double));
				double *component3 = calloc( imageSize, sizeof(double));

				double *component4 = calloc( imageSize, sizeof(double));
				double *component5 = calloc( imageSize, sizeof(double));
				double *component6 = calloc( imageSize, sizeof(double));
				double *component7 = calloc( imageSize, sizeof(double));
				double *component8 = calloc( imageSize, sizeof(double));
				double *component9 = calloc( imageSize, sizeof(double));

				double *component10 = calloc( imageSize, sizeof(double));

				// group 1
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "300", radialFunctionType, radialFunctionOrder, U, accuracy, component1, 31);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "030", radialFunctionType, radialFunctionOrder, U, accuracy, component2, 32);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "003", radialFunctionType, radialFunctionOrder, U, accuracy, component3, 33);	

				powVector(component1, 2, component1, imageSize);
				powVector(component2, 2, component2, imageSize);
				powVector(component3, 2, component3, imageSize);

				// group 2
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "210", radialFunctionType, radialFunctionOrder, U, accuracy, component4, 34);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "201", radialFunctionType, radialFunctionOrder, U, accuracy, component5, 35);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "021", radialFunctionType, radialFunctionOrder, U, accuracy, component6, 36);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "120", radialFunctionType, radialFunctionOrder, U, accuracy, component7, 37);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "102", radialFunctionType, radialFunctionOrder, U, accuracy, component8, 38);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "012", radialFunctionType, radialFunctionOrder, U, accuracy, component9, 39);	

				powVector(component4, 2, component4, imageSize);
				powVector(component5, 2, component5, imageSize);
				powVector(component6, 2, component6, imageSize);
				powVector(component7, 2, component7, imageSize);
				powVector(component8, 2, component8, imageSize);
				powVector(component9, 2, component9, imageSize);
	
				multiplyScalarVector(component4, 3.0, component4, imageSize);	
				multiplyScalarVector(component5, 3.0, component5, imageSize);
				multiplyScalarVector(component6, 3.0, component6, imageSize);
				multiplyScalarVector(component7, 3.0, component7, imageSize);
				multiplyScalarVector(component8, 3.0, component8, imageSize);
				multiplyScalarVector(component9, 3.0, component9, imageSize);

			// group 3
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "111", radialFunctionType, radialFunctionOrder, U, accuracy, component10, 310);

				powVector(component10, 2, component10, imageSize);

				multiplyScalarVector(component10, 6, component10, imageSize);	

				addVector(component1, component2, featureVector, imageSize);
				addVector(featureVector, component3, featureVector, imageSize);	
				addVector(featureVector, component4, featureVector, imageSize);	
				addVector(featureVector, component5, featureVector, imageSize);	
				addVector(featureVector, component6, featureVector, imageSize);	
				addVector(featureVector, component7, featureVector, imageSize);
				addVector(featureVector, component8, featureVector, imageSize);	
				addVector(featureVector, component9, featureVector, imageSize);	
				addVector(featureVector, component10, featureVector, imageSize);		

				sqrtVector(featureVector, featureVector, imageSize);

				free(component1);
				free(component2);
				free(component3);
				free(component4);
				free(component5);
				free(component6);
				free(component7);
				free(component8);
				free(component9);
				free(component10);
			}
			break;

		case 4:
			{
				double *component1 = calloc( imageSize, sizeof(double));
				double *component2 = calloc( imageSize, sizeof(double));
				double *component3 = calloc( imageSize, sizeof(double));	

				double *component4 = calloc( imageSize, sizeof(double));
				double *component5 = calloc( imageSize, sizeof(double));
				double *component6 = calloc( imageSize, sizeof(double));
				double *component7 = calloc( imageSize, sizeof(double));
				double *component8 = calloc( imageSize, sizeof(double));
				double *component9 = calloc( imageSize, sizeof(double));	

				double *component10 = calloc( imageSize, sizeof(double));
				double *component11 = calloc( imageSize, sizeof(double));
				double *component12 = calloc( imageSize, sizeof(double));

				double *component13 = calloc( imageSize, sizeof(double));
				double *component14 = calloc( imageSize, sizeof(double));
				double *component15 = calloc( imageSize, sizeof(double));	

				// group 1
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "400", radialFunctionType, radialFunctionOrder, U, accuracy, component1, 41);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "040", radialFunctionType, radialFunctionOrder, U, accuracy, component2, 42);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "004", radialFunctionType, radialFunctionOrder, U, accuracy, component3, 43);	

				powVector(component1, 2, component1, imageSize);
				powVector(component2, 2, component2, imageSize);
				powVector(component3, 2, component3, imageSize);

				// group 2
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "310", radialFunctionType, radialFunctionOrder, U, accuracy, component4, 44);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "301", radialFunctionType, radialFunctionOrder, U, accuracy, component5, 45);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "031", radialFunctionType, radialFunctionOrder, U, accuracy, component6, 46);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "130", radialFunctionType, radialFunctionOrder, U, accuracy, component7, 47);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "103", radialFunctionType, radialFunctionOrder, U, accuracy, component8, 48);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "013", radialFunctionType, radialFunctionOrder, U, accuracy, component9, 49);	

				powVector(component4, 2, component4, imageSize);
				powVector(component5, 2, component5, imageSize);
				powVector(component6, 2, component6, imageSize);
				powVector(component7, 2, component7, imageSize);
				powVector(component8, 2, component8, imageSize);
				powVector(component9, 2, component9, imageSize);

				multiplyScalarVector(component4, 4.0, component4, imageSize);	
				multiplyScalarVector(component5, 4.0, component5, imageSize);
				multiplyScalarVector(component6, 4.0, component6, imageSize);
				multiplyScalarVector(component7, 4.0, component7, imageSize);
				multiplyScalarVector(component8, 4.0, component8, imageSize);
				multiplyScalarVector(component9, 4.0, component9, imageSize);	

				// group 3
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "220", radialFunctionType, radialFunctionOrder, U, accuracy, component10, 410);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "202", radialFunctionType, radialFunctionOrder, U, accuracy, component11, 411);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "022", radialFunctionType, radialFunctionOrder, U, accuracy, component12, 412);	

				powVector(component10, 2, component10, imageSize);
				powVector(component11, 2, component11, imageSize);
				powVector(component12, 2, component12, imageSize);	

				multiplyScalarVector(component10, 6.0, component10, imageSize);
				multiplyScalarVector(component11, 6.0, component11, imageSize);
				multiplyScalarVector(component12, 6.0, component12, imageSize);	

				// group 4

				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "211", radialFunctionType, radialFunctionOrder, U, accuracy, component13, 413);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "121", radialFunctionType, radialFunctionOrder, U, accuracy, component14, 414);
				calcStencilAndConvolveAndAddResult(image, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, l, "112", radialFunctionType, radialFunctionOrder, U, accuracy, component15, 415);

				powVector(component13, 2, component13, imageSize);
				powVector(component14, 2, component14, imageSize);
				powVector(component15, 2, component15, imageSize);

				multiplyScalarVector(component13, 12.0, component13, imageSize);
				multiplyScalarVector(component14, 12.0, component14, imageSize);
				multiplyScalarVector(component15, 12.0, component15, imageSize);	

				addVector(component1, component2, featureVector, imageSize);
				addVector(featureVector, component3, featureVector, imageSize);	
				addVector(featureVector, component4, featureVector, imageSize);	
				addVector(featureVector, component5, featureVector, imageSize);	
				addVector(featureVector, component6, featureVector, imageSize);	
				addVector(featureVector, component7, featureVector, imageSize);
				addVector(featureVector, component8, featureVector, imageSize);	
				addVector(featureVector, component9, featureVector, imageSize);	
				addVector(featureVector, component10, featureVector, imageSize);
				addVector(featureVector, component11, featureVector, imageSize);	
				addVector(featureVector, component12, featureVector, imageSize);
				addVector(featureVector, component13, featureVector, imageSize);	
				addVector(featureVector, component14, featureVector, imageSize);	
				addVector(featureVector, component15, featureVector, imageSize);		

				sqrtVector(featureVector, featureVector, imageSize);		
			}
			break;
			

		default:
			die("\nERROR: l is not valid\n");
			break;
	}

	char convolveResultFilename[128];

	if (radialFunctionType == 1)
	{
		snprintf(convolveResultFilename, 128, "%s_%d_%f.csv", "MCSH_feature", l, rCutoff);// change this
	}
	else if (radialFunctionType == 2)
	{
		snprintf(convolveResultFilename, 128, "%s_%d_%f_Legendre_%d.csv", "MCSH_feature", l, rCutoff, radialFunctionOrder);// change this
	}
	
	// printf(DensFilename);
	writeMatToFile(convolveResultFilename, featureVector, imageDimX, imageDimY, imageDimZ);

}


double scoreTask(const double rCutoff, const int l)//, const int group)
{
	double result;
	switch (l) 
	{
		case 0:
			result = 1;
			break;

		case 1:
			result = 3;
			break;

		case 2:
			result = 6;
			break;


		case 3:
			result = 10;
			break;
			
		case 4:
			result = 15;
			break;


		default:
			printf("\n***** WARNING: l not valid *****\n");
			result = 1;
			break;

	}

	result *= rCutoff * rCutoff * rCutoff;
	return result;
}



void LegendrePolynomial(const double *x, const double *y, const double *z, const int polynomialOrder, const double rCutoff, double *result, const int size)
{
	// (2*r_array-r)/r
	double *r = calloc( size, sizeof(double));
	getRArray(x, y, z, r, size);

	int i;
	for (i = 0; i < size; i++)
	{
		r[i] = (2.0 * r[i] - rCutoff) / rCutoff;
	}


	if (polynomialOrder == 0){
		// 1
		for ( i = 0; i < size; i++)
		{
			result[i] = 1.0;
		}
	} else if (polynomialOrder == 1){
		// x
		for ( i = 0; i < size; i++)
		{
			result[i] = r[i];
		}
	} else if (polynomialOrder == 2){
		// 0.5 * (3*x*x - 1)
		double *temp1 = malloc( size * sizeof(double));

		polyArray(r, 2, 3.0, temp1, size);

		addScalarVector(temp1, -1.0, result, size);
		multiplyScalarVector(result, 0.5, result, size);

		free(temp1);
	} else if (polynomialOrder == 3){
		// 0.5 * (5*x*x*x - 3x)
		double *temp1 = malloc( size * sizeof(double));
		double *temp2 = malloc( size * sizeof(double));

		polyArray(r, 3, 5.0, temp1, size);

		multiplyScalarVector(r, -3.0, temp2, size);

		addVector(temp1, temp2, result, size);
		multiplyScalarVector(result, 0.5, result, size);

		free(temp1);
		free(temp2);
	} else if (polynomialOrder == 4){
		// (1/8) * (35*x*x*x*x - 30*x*x +3)
		double *temp1 = malloc( size * sizeof(double));
		double *temp2 = malloc( size * sizeof(double));

		polyArray(r, 4, 35.0, temp1, size);
		polyArray(r, 2, -30.0, temp2, size);

		addVector(temp1, temp2, result, size);
		addScalarVector(result, 3.0, result, size);
		multiplyScalarVector(result, (1.0/8.0), result, size);

		free(temp1);
		free(temp2);
	} else if (polynomialOrder == 5){
		// (1/8) * (63*x*x*x*x*x - 70*x*x*x + 15*x)
		double *temp1 = malloc( size * sizeof(double));
		double *temp2 = malloc( size * sizeof(double));
		double *temp3 = malloc( size * sizeof(double));

		polyArray(r, 5, 63.0, temp1, size);
		polyArray(r, 3, -70.0, temp2, size);
		multiplyScalarVector(r, 15.0, temp3, size);

		addVector(temp1, temp2, result, size);
		addVector(result, temp3, result, size);
		multiplyScalarVector(result, (1.0/8.0), result, size);

		free(temp1);
		free(temp2);
		free(temp3);
	} else if (polynomialOrder == 6){
		// (1/16) * (231*x*x*x*x*x*x - 315*x*x*x*x + 105*x*x -5)
		double *temp1 = malloc( size * sizeof(double));
		double *temp2 = malloc( size * sizeof(double));
		double *temp3 = malloc( size * sizeof(double));

		polyArray(r, 6, 231.0, temp1, size);
		polyArray(r, 4, -315.0, temp2, size);
		polyArray(r, 2, 105.0, temp3, size);

		addVector(temp1, temp2, result, size);
		addVector(result, temp3, result, size);
		addScalarVector(result, -5.0, result, size);
		multiplyScalarVector(result, (1.0/16.0), result, size);

		free(temp1);
		free(temp2);
		free(temp3);
	} else if (polynomialOrder == 7){
		// (1/16) * (429*x*x*x*x*x*x*x - 693*x*x*x*x*x + 315*x*x*x - 35*x)
		double *temp1 = malloc( size * sizeof(double));
		double *temp2 = malloc( size * sizeof(double));
		double *temp3 = malloc( size * sizeof(double));
		double *temp4 = malloc( size * sizeof(double));

		polyArray(r, 7, 231.0, temp1, size);
		polyArray(r, 5, -315.0, temp2, size);
		polyArray(r, 3, 105.0, temp3, size);
		multiplyScalarVector(r, 15.0, temp4, size);

		addVector(temp1, temp2, result, size);
		addVector(result, temp3, result, size);
		addVector(result, temp4, result, size);
		multiplyScalarVector(result, (1.0/16.0), result, size);

		free(temp1);
		free(temp2);
		free(temp3);
		free(temp4);
	} else if (polynomialOrder == 8){
		// (1/128) * (6435*x**8 - 12012*x**6 + 6930*x**4 - 1260*x*2 + 35)
		double *temp1 = malloc( size * sizeof(double));
		double *temp2 = malloc( size * sizeof(double));
		double *temp3 = malloc( size * sizeof(double));
		double *temp4 = malloc( size * sizeof(double));

		polyArray(r, 8, 6435.0, temp1, size);
		polyArray(r, 6, -12012.0, temp2, size);
		polyArray(r, 4, 6930.0, temp3, size);
		polyArray(r, 2, -1260.0, temp4, size);

		addVector(temp1, temp2, result, size);
		addVector(result, temp3, result, size);
		addVector(result, temp4, result, size);
		addScalarVector(result, 35.0, result, size);
		multiplyScalarVector(result, (1.0/128.0), result, size);

		free(temp1);
		free(temp2);
		free(temp3);
		free(temp4);
	} else if (polynomialOrder == 9){
		// (1/128) * (12155*x**9 - 25740*x**7 + 18018*x**5 - 4620*x**3 + 315x)
		double *temp1 = malloc( size * sizeof(double));
		double *temp2 = malloc( size * sizeof(double));
		double *temp3 = malloc( size * sizeof(double));
		double *temp4 = malloc( size * sizeof(double));
		double *temp5 = malloc( size * sizeof(double));

		polyArray(r, 9, 12155.0, temp1, size);
		polyArray(r, 7, -25740.0, temp2, size);
		polyArray(r, 5, 18018.0, temp3, size);
		polyArray(r, 3, -4620.0, temp4, size);
		multiplyScalarVector(r, 315.0, temp5, size);

		addVector(temp1, temp2, result, size);
		addVector(result, temp3, result, size);
		addVector(result, temp4, result, size);
		addVector(result, temp5, result, size);
		multiplyScalarVector(result, (1.0/128.0), result, size);

		free(temp1);
		free(temp2);
		free(temp3);
		free(temp4);
		free(temp5);
	} else if (polynomialOrder == 10){
		// (1/256) * (46189*x**10 - 109395*x**8 + 90090*x**6 - 30030*x**4 + 3465*x**2 -63)
		double *temp1 = malloc( size * sizeof(double));
		double *temp2 = malloc( size * sizeof(double));
		double *temp3 = malloc( size * sizeof(double));
		double *temp4 = malloc( size * sizeof(double));
		double *temp5 = malloc( size * sizeof(double));

		polyArray(r, 10, 46189.0, temp1, size);
		polyArray(r, 8, -109395.0, temp2, size);
		polyArray(r, 6, 90090.0, temp3, size);
		polyArray(r, 4, 30030.0, temp4, size);
		polyArray(r, 2, 3465.0, temp5, size);

		addVector(temp1, temp2, result, size);
		addVector(result, temp3, result, size);
		addVector(result, temp4, result, size);
		addVector(result, temp5, result, size);
		addScalarVector(result, -63.0, result, size);
		multiplyScalarVector(result, (1.0/256.0), result, size);

		free(temp1);
		free(temp2);
		free(temp3);
		free(temp4);
		free(temp5);
	} else {
		printf("\nERROR: Legendre Order Not Valid\n");
	}

	for (i = 0; i < size; i++)
	{	
		if (r[i] > rCutoff)
		{
			result[i] = 0.0;
		}
	}

	free(r);

}







