# include <stdio.h>
# include <stdlib.h>
# include <string.h>
#include <assert.h>
#include <time.h> 
#include <mpi.h>
# include <math.h>
/* BLAS, LAPACK, LAPACKE routines */
#ifdef USE_MKL
    // #define MKL_Complex16 double complex
    #include <mkl.h>
#else
    #include <cblas.h>
#endif

#include "MCSHHelper.h"



void die(const char *message)
{
	printf( "ERROR: %s\n", message );
	exit(1);
}

void GetPlaneEquation(const Point p1, const Point p2, const Point p3, double *a, double *b, double *c, double *d)
{


	double a1 = p2.x - p1.x;
	double b1 = p2.y - p1.y;
	double c1 = p2.z - p1.z;
	double a2 = p3.x - p1.x;
	double b2 = p3.y - p1.y;
	double c2 = p3.z - p1.z;

	double tempA = b1 * c2 - b2 * c1;
	double tempB = a2 * c1 - a1 * c2;
	double tempC = a1 * b2 - b1 * a2;
	double tempD = - tempA * p1.x - tempB * p1.y - tempC * p1.z;

	*a = tempA;
	*b = tempB;
	*c = tempC;
	*d = tempD;
}


int CheckPlaneIntersectWithSphere(const Point p1, const Point p2, const Point p3, const double rCutoff, const Point origin)
{
	double xs = origin.x;
	double ys = origin.y;
	double zs = origin.z;
	double a, b, c, d; 
	int result;
	GetPlaneEquation(p1, p2, p3, &a, &b, &c, &d);

	// printf("a: %22.15f \t b: %22.15f \t c: %22.15f \t d: %22.15f \n", a,b,c,d);
	// printf("xs: %22.15f \t ys: %22.15f \t zs: %22.15f \n", xs,ys,zs);

	double numerator = fabs(a * xs + b * ys + c * zs + d);
	double denominator = sqrt(a*a + b*b + c*c);
	double test_d = numerator / denominator;
	// printf("numerator: %22.15f \t denominator: %22.15f\n", numerator, denominator);

	// double test_d = abs(a * xs + b * ys + c * zs + d) / sqrt(a*a + b*b + c*c);
	// printf("numerator: %22.15f \t denominator: %22.15f\n", abs(a * xs + b * ys + c * zs + d), sqrt(a*a + b*b + c*c));

	// double test_d = abs(a * xs + b * ys + c * zs + d) / sqrt(pow(a,2.0) + pow(b,2.0) + pow(c,2.0));
	// printf("numerator: %22.15f \t denominator: %22.15f\n", abs(a * xs + b * ys + c * zs + d), sqrt(pow(a,2.0) + pow(b,2.0) + pow(c,2.0)));
	// printf("testd: %22.15f \t rCutoff: %22.15f\n", test_d, rCutoff);

	if (rCutoff > test_d)
	{
	    result = 1;
	}
	else
	{
	    result = 0;
	}

	// rCutoff > test_d ? result = 1: result = 0;
	
	return result;
}


Point UTransform(const double x, const double y, const double z, const double *U)
{
	double tempX = x * U[0] + y * U[1] + z * U[2];
	double tempY = x * U[3] + y * U[4] + z * U[5];
	double tempZ = x * U[6] + y * U[7] + z * U[8];

	Point p = {tempX,tempY,tempZ};

	return p;
}


void GetDimensionsPlane(const double hx, const double hy, const double hz, const double rCutoff, const double *U, int *dimXResult, int *dimYResult, int *dimZResult)
{
	int dimX = 2 * ceil(rCutoff / hx) + 1;
	int dimY = 2 * ceil(rCutoff / hy) + 1;
	int dimZ = 2 * ceil(rCutoff / hz) + 1;

	Point origin = {0,0,0};

	double ref_x_min, ref_x_max, ref_y_min, ref_y_max, ref_z_min, ref_z_max;
	Point p1x_1, p2x_1, p3x_1, p1x_2, p2x_2, p3x_2;  

	while (1)
	{
		ref_x_min = - hx * dimX * 0.5;
		ref_x_max = hx * dimX * 0.5;
		ref_y_min = - hy * dimY * 0.5;
		ref_y_max = hy * dimY * 0.5;
		ref_z_min = - hz * dimZ * 0.5;
		ref_z_max = hz * dimZ * 0.5;

		p1x_1 = UTransform(ref_x_min,ref_y_min,ref_z_min,U);
		p2x_1 = UTransform(ref_x_min,ref_y_max,ref_z_min,U);
		p3x_1 = UTransform(ref_x_min,ref_y_min,ref_z_max,U);

		p1x_2 = UTransform(ref_x_max,ref_y_min,ref_z_min,U);
		p2x_2 = UTransform(ref_x_max,ref_y_max,ref_z_min,U);
		p3x_2 = UTransform(ref_x_max,ref_y_min,ref_z_max,U);

		// printf("%22.15f \t %22.15f \t %22.15f \t %22.15f \t %22.15f \t %22.15f\n",ref_x_min,ref_x_max,ref_y_min,ref_y_max,ref_z_min,ref_z_max);
		// printf("%22.15f \t %22.15f \t %22.15f\n",p1x_1.x, p1x_1.y,p1x_1.z);
		// printf("%22.15f \t %22.15f \t %22.15f\n",p2x_1.x, p2x_1.y,p2x_1.z);
		// printf("%22.15f \t %22.15f \t %22.15f\n",p3x_1.x, p3x_1.y,p3x_1.z);
		// printf("%22.15f \t %22.15f \t %22.15f\n",p1x_2.x, p1x_2.y,p1x_2.z);
		// printf("%22.15f \t %22.15f \t %22.15f\n",p2x_2.x, p2x_2.y,p2x_2.z);
		// printf("%22.15f \t %22.15f \t %22.15f\n",p3x_2.x, p3x_2.y,p3x_2.z);

		if (CheckPlaneIntersectWithSphere(p1x_1, p2x_1, p3x_1, rCutoff, origin) || CheckPlaneIntersectWithSphere(p1x_2, p2x_2, p3x_2, rCutoff, origin))
		{
			dimX += 2;
			// printf("added dimension to X\n");
		}
		else
		{
			break;
		}
	}

	Point p1y_1, p2y_1, p3y_1, p1y_2, p2y_2, p3y_2;  

	while (1)
	{
		ref_x_min = - hx * dimX * 0.5;
		ref_x_max = hx * dimX * 0.5;
		ref_y_min = - hy * dimY * 0.5;
		ref_y_max = hy * dimY * 0.5;
		ref_z_min = - hz * dimZ * 0.5;
		ref_z_max = hz * dimZ * 0.5;

		p1y_1 = UTransform(ref_x_min,ref_y_min,ref_z_min,U);
		p2y_1 = UTransform(ref_x_max,ref_y_min,ref_z_min,U);
		p3y_1 = UTransform(ref_x_min,ref_y_min,ref_z_max,U);

		p1y_2 = UTransform(ref_x_min,ref_y_max,ref_z_min,U);
		p2y_2 = UTransform(ref_x_max,ref_y_max,ref_z_min,U);
		p3y_2 = UTransform(ref_x_min,ref_y_max,ref_z_max,U);

		if (CheckPlaneIntersectWithSphere(p1y_1, p2y_1, p3y_1, rCutoff, origin) || CheckPlaneIntersectWithSphere(p1y_2, p2y_2, p3y_2, rCutoff, origin))
		{
			dimY += 2;
			// printf("added dimension to Y\n");
		}
		else
		{
			break;
		}
	}

	Point p1z_1, p2z_1, p3z_1, p1z_2, p2z_2, p3z_2;  

	while (1)
	{
		ref_x_min = - hx * dimX * 0.5;
		ref_x_max = hx * dimX * 0.5;
		ref_y_min = - hy * dimY * 0.5;
		ref_y_max = hy * dimY * 0.5;
		ref_z_min = - hz * dimZ * 0.5;
		ref_z_max = hz * dimZ * 0.5;

		p1z_1 = UTransform(ref_x_min,ref_y_min,ref_z_min,U);
		p2z_1 = UTransform(ref_x_max,ref_y_min,ref_z_min,U);
		p3z_1 = UTransform(ref_x_min,ref_y_max,ref_z_min,U);

		p1z_2 = UTransform(ref_x_min,ref_y_min,ref_z_max,U);
		p2z_2 = UTransform(ref_x_max,ref_y_min,ref_z_max,U);
		p3z_2 = UTransform(ref_x_min,ref_y_max,ref_z_max,U);

		if (CheckPlaneIntersectWithSphere(p1z_1, p2z_1, p3z_1, rCutoff, origin) || CheckPlaneIntersectWithSphere(p1z_2, p2z_2, p3z_2, rCutoff, origin))
		{
			dimZ += 2;
			// printf("added dimension to Z\n");
		}
		else
		{
			break;
		}
	}

	*dimXResult = dimX;
	*dimYResult = dimY;
	*dimZResult = dimZ;
}


void printArr(double *arr, int size)
{
	int i;
	for ( i = 0; i < size; i++ )
	{
		printf("%f\n", arr[i]);
	}
}

void sqrtVector(const double *arr, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		result[i] = sqrt(arr[i]);
	}
}


void powVector(const double *arr, const int power, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		result[i] = pow(arr[i], power);
	}
}

void addVector(const double *arr1, const double *arr2, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		result[i] = arr1[i] + arr2[i];
	}
}

void subtractVector(const double *arr1, const double *arr2, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		result[i] = arr1[i] - arr2[i];
	}
}

void divideVector(const double *arr1, const double *arr2, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		result[i] = arr1[i] / arr2[i];
	}
}

void addScalarVector(const double *arr, const double a, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		result[i] = arr[i] + a;
	}
}

void multiplyScalarVector(const double *arr1, const double a, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		result[i] = arr1[i] * a;
	}
}

void getRArray(const double *x, const double *y, const double *z, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		double r = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
		result[i] = r;
	}
}

void polyArray(const double *x, const int powX, const double a, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		result[i] = a * pow(x[i], powX);
	}
}

void polyXYZArray(const double *x, const double *y, const double *z, const int powX, const int powY, const int powZ, const double a, double *result, const int size)
{
	int i = 0;
	for (i = 0; i < size; i++)
	{
		result[i] = a * pow(x[i], powX) * pow(y[i], powY) * pow(z[i], powZ) ;
	}
}

void applyU(double *X, double *Y, double*Z, const double *U, const int size)
{
	double *combinedArr = malloc( 3 * size * sizeof(double));
	double *matMutResult = malloc( 3 * size * sizeof(double));

	int i;
	for ( i = 0; i < size; i++ )
	{
		combinedArr[ i * 3 ] = X[i];
		combinedArr[ i * 3 + 1 ] = Y[i];
		combinedArr[ i * 3 + 2 ] = Z[i];
	}
	// printf("----------combined-----------\n");
	// printArr(combinedArr, size*3);

	int overallSize = size * 3;
	// cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, overallSize, 3, 3, 
 //                1.0, combinedArr, overallSize, U, 3, 0.0, matMutResult, overallSize);


	for ( i = 0; i < size; i++ )
	{
		X[i] = matMutResult[ i * 3 ];
		Y[i] = matMutResult[ i * 3 + 1 ];
		Z[i] = matMutResult[ i * 3 + 2 ];
	}

	free(combinedArr);
	free(matMutResult);

	//void cblas_dgemm ( Layout, transa,  transb,  m,  n, k, 
	//			alpha, *a,  lda,  *b, ldb, beta, *c, ldc);

}

void applyU2(double *X, double *Y, double*Z, const double *U, const int size)
{
	// double *resultX = malloc( size * sizeof(double));
	// double *resultY = malloc( size * sizeof(double));
	// double *resultZ = malloc( size * sizeof(double));
	int i;
	double tempX, tempY, tempZ;
	for (i = 0; i < size; i++)
	{
		tempX = U[0] * X[i] + U[1] * Y[i] + U[2] * Z[i];
		tempY = U[3] * X[i] + U[4] * Y[i] + U[5] * Z[i];
		tempZ = U[6] * X[i] + U[7] * Y[i] + U[8] * Z[i];

		X[i] = tempX;
		Y[i] = tempY;
		Z[i] = tempZ;

		// resultX[i] = U[0] * X[i] + U[1] * Y[i] + U[2] * Z[i];
		// resultY[i] = U[3] * X[i] + U[4] * Y[i] + U[5] * Z[i];
		// resultZ[i] = U[6] * X[i] + U[7] * Y[i] + U[8] * Z[i];
	}

	// for (i = 0; i < size; i++)
	// {
	// 	X[i] = resultX[i];
	// 	Y[i] = resultX[i];
	// 	Z[i] = resultX[i];
	// }

	// free(resultX);
	// free(resultY);
	// free(resultZ);
}



double calcDv(const double hx, const double hy, const double hz, const int accuracy, const double *U)
{
	double hxAcc = hx / accuracy;
	double hyAcc = hy / accuracy;
	double hzAcc = hz / accuracy;

	Point l1 = {U[0] * hxAcc, U[1] * hxAcc, U[2] * hxAcc};
	Point l2 = {U[3] * hyAcc, U[4] * hyAcc, U[5] * hyAcc};
	Point l3 = {U[6] * hzAcc, U[7] * hzAcc, U[8] * hzAcc};

	Point crossL2L3 = {l2.y * l3.z - l2.z * l3.y, l2.z * l3.x - l2.x * l3.z, l2.x * l3.z - l2.z * l3.x};

	double result = fabs(l1.x * crossL2L3.x + l1.y * crossL2L3.y + l1.z * crossL2L3.z);

	return result;
}


double sumArr(const double *arr, const int size)
{	
	double result = 0;
	int i;
	for (i = 0; i < size; i++)
	{
		result += arr[i];
	}
	return result;

}



void linspace(double start, double end, double *result, int num)
{	
	double stepsize = (end - start) / (double)(num-1);
	double current = start;

	int i;
	for (i = 0; i < num; i++)
	{
		result[i] = current;
		current += stepsize;
	}
}

void meshgrid3D(const double *x, const double *y, const double *z, const int sizex, const int sizey, const int sizez, double *X, double *Y, double *Z)
{
	int i,j,k;
	
	for (k = 0; k < sizez; k++){
		for ( j = 0; j < sizey; j++) {
			for ( i = 0; i < sizey; i++) {
				X[ k * sizex * sizey + j * sizex + i ] = x[i];
				Y[ k * sizex * sizey + j * sizex + i ] = y[j];
				Z[ k * sizex * sizey + j * sizex + i ] = z[k];
				// printf("%f \t %f \t %f \n", x[i], y[j],z[k]);
			}
		}
	}
}

void getCentralCoords(const double hx, const double hy, const double hz, const int accuracy, double *refX, double *refY, double *refZ)
{
	double hxAcc = hx / accuracy;
	double hyAcc = hy / accuracy;
	double hzAcc = hz / accuracy;

	double *ref_x_li = calloc( accuracy, sizeof(double));
	double *ref_y_li = calloc( accuracy, sizeof(double));
	double *ref_z_li = calloc( accuracy, sizeof(double));

	double xStart = -((hx / 2.0) - (hxAcc / 2.0));
	double xEnd = (hx / 2.0) - (hxAcc / 2.0);
	linspace(xStart, xEnd, ref_x_li, accuracy);

	double yStart = -((hy / 2.0) - (hyAcc / 2.0));
	double yEnd = (hy / 2.0) - (hyAcc / 2.0);
	linspace(yStart, yEnd, ref_y_li, accuracy);

	double zStart = -((hz / 2.0) - (hzAcc / 2.0));
	double zEnd = (hz / 2.0) - (hzAcc / 2.0);
	linspace(zStart, zEnd, ref_z_li, accuracy);

	meshgrid3D(ref_x_li, ref_y_li, ref_z_li, accuracy, accuracy, accuracy, refX, refY, refZ);

	free(ref_x_li);
	free(ref_y_li);
	free(ref_z_li);

}

double calcNorm3(double x1, double x2, double x3)
{
	double result = sqrt(x1*x1 + x2*x2 + x3*x3);
	return result;
}

void normalizeU(double *U, double *normalizedU)
{
	// double norm1 = calcNorm3(U[0], U[1], U[2]);
	// double norm2 = calcNorm3(U[3], U[4], U[5]);
	// double norm3 = calcNorm3(U[6], U[7], U[8]);

	// normalizedU[0] = U[0] / norm1;
	// normalizedU[1] = U[1] / norm1;
	// normalizedU[2] = U[2] / norm1;

	// normalizedU[3] = U[3] / norm2;
	// normalizedU[4] = U[4] / norm2;
	// normalizedU[5] = U[5] / norm2;

	// normalizedU[6] = U[6] / norm3;
	// normalizedU[7] = U[7] / norm3;
	// normalizedU[8] = U[8] / norm3;


	double norm1 = calcNorm3(U[0], U[3], U[6]);
	double norm2 = calcNorm3(U[1], U[4], U[7]);
	double norm3 = calcNorm3(U[2], U[5], U[8]);

	normalizedU[0] = U[0] / norm1;
	normalizedU[1] = U[1] / norm2;
	normalizedU[2] = U[2] / norm3;

	normalizedU[3] = U[3] / norm1;
	normalizedU[4] = U[4] / norm2;
	normalizedU[5] = U[5] / norm3;

	normalizedU[6] = U[6] / norm1;
	normalizedU[7] = U[7] / norm2;
	normalizedU[8] = U[8] / norm3;
}

int mod(int a, int b)
{
	if (a<0)
	{
		return a + b;
	}
	else if (a >= b)
	{
		return a - b;
	}
	else
	{
		return a;
	}
    // int r = a % b;
    // return r < 0 ? r + b : r;
}

void calcSingleConvolveStep(const double *image, const double stencilCoeff, const int shiftX, const int shiftY, const int shiftZ, double *result, const int imageSize, const int imageDimX, const int imageDimY, const int imageDimZ)
{
	double *tempResult = calloc( imageSize, sizeof(double));
	multiplyScalarVector(image, stencilCoeff, tempResult, imageSize);
	
	// printf("image %10f \t coeff %10f \t result %10f \t temp result: %10f \t", image[100], stencilCoeff, result[100], tempResult[100]);

	// int i = 0, index;
	// for (i = 0; i < imageSize; i++)
	// {
	// 	newIndex = i;
	// 	result[newIndex] += tempResult[i];
	// }
	// free(tempResult);
	// printf("end multiply scalar\n");
	// printf("image size %d  %d  %d \n", imageDimX, imageDimY, imageDimZ);
	int i, j, k, new_i, new_j, new_k, original_index, new_index;
	// double tempResult;
	for (k = 0; k < imageDimZ; k++){
		for ( j = 0; j < imageDimY; j++) {
			for ( i = 0; i < imageDimX; i++) {

				// printf("start add %d  %d  %d\n", i, j, k);
				original_index = k * imageDimX * imageDimY + j * imageDimX + i;
				new_i = mod ((i - shiftX), imageDimX);
				new_j = mod ((j - shiftY), imageDimY);
				new_k = mod ((k - shiftZ), imageDimZ);
				new_index = new_k * imageDimX * imageDimY + new_j * imageDimX + new_i;


				// tempResult = stencilCoeff * image[new_index] + result[original_index];
				// //printf("new index: %d \t %d \t %d \t\t shift: %d \t %d \t %d \t\t old index: %d \t %d \t %d \t\t image: %10f \t %10f \t %10f\n", new_i, new_j, new_k, shiftX, shiftY, shiftZ, i, j, k, image[new_index], tempResult, result[original_index]);
				// result[original_index] = tempResult;


				// printf("end calc new index %d  %d  %d     %d\n", new_i, new_j, new_k, new_index);
				result[original_index] += tempResult[new_index];
				// printf("end add %d  %d  %d \t%10f \t%10f\n", i, j, k, result[original_index], tempResult[new_index]);
				// printf("end add %d  %d  %d\n", i, j, k);
			}
		}
	}
	// printf("added result %10f \t", result[100]);

	//free(tempResult);
}


void convolve(const double *image, const double *stencil, const int imageDimX, const int imageDimY, const int imageDimZ, const int stencilDimX, const int stencilDimY, const int stencilDimZ, double *result)
{
	int imageSize = imageDimX * imageDimY * imageDimZ;
	int stencilSize = stencilDimX * stencilDimY * stencilDimZ;

	int *xShiftList = malloc( stencilSize * sizeof(int));
	int *yShiftList = malloc( stencilSize * sizeof(int));
	int *zShiftList = malloc( stencilSize * sizeof(int));
	double *coeffList = malloc( stencilSize * sizeof(double));

	

	int i, j, k, index = 0;
	for (k = 0; k < stencilDimZ; k++){
		for ( j = 0; j < stencilDimY; j++) {
			for ( i = 0; i < stencilDimX; i++) {
				// index = k * stencilDimX * stencilDimY + j * stencilDimX + i;
				// ((nx+1)/2)-i, ((ny+1)/2)-j, ((nz+1)/2)-k
				xShiftList[index] = ((stencilDimX - 1) / 2) - i;
				yShiftList[index] = ((stencilDimY - 1) / 2) - j;
				zShiftList[index] = ((stencilDimZ - 1) / 2) - k;
				coeffList[index] = stencil[index];
				index++;

				// fprintf(output_fp,"%d,%d,%d,%22f\n",i,j,k,stencil[index]);
				// printf("%10.8f\t",stencil[index]);
			}
		}
	} 

	// printf("end ordering \n");


	for (i = 0; i < stencilSize; i++)
	{
		// printf("start convolve step %d \n", i);
		//const double *image, const double kernelCoeff, const int shiftX, const int shiftY, const int shiftZ, double *result, const int imageSize, const int kernelDimX, const int kernelDimY, const int kernelDimZ
		calcSingleConvolveStep(image, coeffList[i], xShiftList[i], yShiftList[i], zShiftList[i], result, imageSize, imageDimX, imageDimY, imageDimZ);
		//printf("after result %10f \n", result[100]);
		// printf("end convolve step %d \n", i);
	}

	// printf("\nafter everything %10f \n\n", result[100]);


	free(xShiftList);
	free(yShiftList);
	free(zShiftList);
	free(coeffList);

}

void convolve2(const double *image, const double *stencil, const int imageDimX, const int imageDimY, const int imageDimZ, const int stencilDimX, const int stencilDimY, const int stencilDimZ, double *result)
{
	int imageSize = imageDimX * imageDimY * imageDimZ;
	int stencilSize = stencilDimX * stencilDimY * stencilDimZ;

	int i, j, k, l, m, n;
	int stencilCenterX = (stencilDimX - 1) / 2;
	int stencilCenterY = (stencilDimY - 1) / 2;
	int stencilCenterZ = (stencilDimZ - 1) / 2;
	int xShift,yShift,zShift;
	int outIndex, outI, outJ, outK;
	int imageIndex = 0, stencilIndex = 0;

	for (k = 0; k < imageDimZ; k++){
		for ( j = 0; j < imageDimY; j++) {
			for ( i = 0; i < imageDimX; i++) {

				stencilIndex = 0;
				for (n = 0; n < stencilDimZ; n++){
					for ( m = 0; m < stencilDimY; m++) {
						for ( l = 0; l < stencilDimX; l++) {
							// xShift = stencilCenterX - l;
							// yShift = stencilCenterY - m;
							// zShift = stencilCenterZ - n;

							xShift = l - stencilCenterX;
							yShift = m - stencilCenterY;
							zShift = n - stencilCenterZ;

							outI = mod ((i - xShift), imageDimX);
							outJ = mod ((j - yShift), imageDimY);
							outK = mod ((k - zShift), imageDimZ);
							// printf("%d \t %d \t %d \t\t %d \t %d \t %d \t\t %d \t %d \t %d \t\t %d \t %d \t %d \n",i,j,k,l,m,n,xShift, yShift,zShift,outI,outJ,outK);
							
							outIndex = outK * imageDimX * imageDimY + outJ * imageDimX + outI;
							result[outIndex] += stencil[stencilIndex]* image[imageIndex];
							//result[outIndex] = fma(stencil[stencilIndex], image[imageIndex], result[outIndex]);
							stencilIndex++;

						}
					}
				} 

				imageIndex ++;
			}
		}
	} 
}


void convolve3(const double *image, const double *stencil, const int imageDimX, const int imageDimY, const int imageDimZ, const int stencilDimX, const int stencilDimY, const int stencilDimZ, double *result)
{
	int imageSize = imageDimX * imageDimY * imageDimZ;
	int stencilSize = stencilDimX * stencilDimY * stencilDimZ;

	int i, j, k, l, m, n;
	int stencilCenterX = (stencilDimX - 1) / 2;
	int stencilCenterY = (stencilDimY - 1) / 2;
	int stencilCenterZ = (stencilDimZ - 1) / 2;
	int xShift,yShift,zShift;
	int outIndex, outI, outJ, outK;
	int imageIndex = 0, stencilIndex = 0;

	for (k = 0; k < imageDimZ; k++){
		for ( j = 0; j < imageDimY; j++) {
			for ( i = 0; i < imageDimX; i++) {

				stencilIndex = 0;
				for (n = 0, zShift = - stencilCenterZ; n < stencilDimZ; n++, zShift++){
					outK = mod ((k - zShift), imageDimZ);
					for ( m = 0, yShift = - stencilCenterY; m < stencilDimY; m++, yShift++) {
						outJ = mod ((j - yShift), imageDimY);
						for ( l = 0, xShift = - stencilCenterX; l < stencilDimX; l++, xShift++) {
							outI = mod ((i - xShift), imageDimX);
							// printf("%d \t %d \t %d \t\t %d \t %d \t %d \t\t %d \t %d \t %d \t\t %d \t %d \t %d \n",i,j,k,l,m,n,xShift, yShift,zShift,outI,outJ,outK);
							
							outIndex = outK * imageDimX * imageDimY + outJ * imageDimX + outI;
							result[outIndex] += stencil[stencilIndex]* image[imageIndex];
							//result[outIndex] = fma(stencil[stencilIndex], image[imageIndex], result[outIndex]);
							stencilIndex++;
						}
					}
				}
				imageIndex ++;
			}
		}
	} 
}


void writeMatToFile(const char *filename, const double *data, const int dimX, const int dimY, const int dimZ)
{
	
	FILE *output_fp = fopen(filename,"w");
	if (output_fp == NULL) {
		printf("\nCannot open file \"%s\"\n",filename);
		die("cant open file");
	}

	int i,j,k, index = 0;
	for (k = 0; k < dimZ; k++){
		for ( j = 0; j < dimY; j++) {
			for ( i = 0; i < dimX; i++) {
				fprintf(output_fp,"%d,%d,%d,%.15f\n",i,j,k,data[index]);
				index ++;
			}
		}
	}


	fclose(output_fp);
}












void getMainParameter(const double rStepsize, const double rMaxCutoff, const int maxOrder, const int length, double* rCutoffList, int* lList, int* groupList)
{
	// length
	// int length = getDescriptorListLength(rStepsize, rMaxCutoff, maxOrder)
	int numGroup = getNumGroup(maxOrder);
	int numRCutoff = getNumRCutoff(rStepsize, rMaxCutoff);

	int i, j, index = 0;

	for (i = 0; i < numGroup; i++)
	{
		for (j = 0; j < numRCutoff; j++)
		{
			rCutoffList[index] = (j+1) * rStepsize;
			lList[index] = getCurrentLNumber(i);
			groupList[index] = getCurrentGroupNumber(i);
			index++;
		}
	}


}

int getNumRCutoff(const double rStepsize, const double rMaxCutoff)
{
	return (int) ceil(rMaxCutoff / rStepsize);
}

int getNumGroup(const int maxOrder)
{
	int numGroupList[5] = {1,1,2,3,4};
	if (maxOrder > 4) die("\n Error: Maximum Order Not Implemented \n");

	int numGroup=0;
	int i;
	for (i = 0; i < maxOrder + 1; i++)
	{
		numGroup += numGroupList[i];
	}
	
	return numGroup;
}

int getDescriptorListLength(const double rStepsize, const double rMaxCutoff, const int maxOrder)
{

	int numGroup = getNumGroup(maxOrder);
	// printf("\nnumber groups:%d \n", numGroup);

	int numRCutoff = getNumRCutoff(rStepsize, rMaxCutoff);
	// printf("\nnumber r cutoff:%d \n", numRCutoff);

	return numGroup * numRCutoff;

}

int getCurrentGroupNumber(const int currentIndex)
{	
	// 1 --> 1
	// 2 --> 1
	// 3 --> 1
	// 4 --> 2
	// 5 --> 1
	// 6 --> 2
	// 7 --> 3
	// 8 --> 1
	// 9 --> 2
	// 10 --> 3
	// 11 --> 4


	int resultList[11] = {1,1,1,2,1,2,3,1,2,3,4};
	return resultList[currentIndex];




	// int numGroupList[5] = {1,1,2,3,4};
	// int L = 0, groupNumber = currentIndex, currentTotalGroups = 0, indexCursor = currentIndex;

	// int i;
	// for (i = 0; i < 5; i++)
	// {
	// 	indexCursor -= numGroupList[i];
	// 	if (indexCursor >= 0) 
	// 		groupNumber -= numGroupList[i];
	// }

	// return groupNumber+1;
}


int getCurrentLNumber(const int currentIndex)
{	
	// 1 --> 0
	// 2 --> 1
	// 3 --> 2
	// 4 --> 2
	// 5 --> 3
	// 6 --> 3
	// 7 --> 3
	// 8 --> 4
	// 9 --> 4
	// 10 --> 4
	// 11 --> 4

	int resultList[11] = {0,1,2,2,3,3,3,4,4,4,4};
	return resultList[currentIndex];


	// int numGroupList[5] = {1,1,2,3,4};
	// int L = 0, indexCursor = currentIndex;

	// int i;
	// for (i = 0; i < 5; i++)
	// {
	// 	indexCursor -= numGroupList[i];
	// 	if (indexCursor > 0) 
	// 		L++;
		
	// }

	// return L;
}
