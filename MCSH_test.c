# include <stdio.h>
# include <stdlib.h>
# include <string.h>
#include <assert.h>
#include <mpi.h>
# include <math.h>
/* BLAS, LAPACK, LAPACKE routines */
//#ifdef USE_MKL
    // #define MKL_Complex16 double complex
   // #include <mkl.h>
//#else
//    #include <cblas.h>
//#endif

#include "MCSHHelper.h"
#include "MCSH.h"
#include "MCSHDescriptorMain.h"
//#include "test_data.h"
//#include "pseudo_dens_CO_60.h"
//#include "pseudo_dens_CO_90.h"
//#include "pseudo_dens_CO_30.h"
#include "/home/ssahoo41/Documents/MCSH_codes/generate_images/pseudo_dens_O2_10.h"

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
	int world_rank, world_size;

	double start = MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);//rank or process ID
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);//total number of processors
	//printf("I am process %i out of %i processes\n", world_rank, world_size);

	int numProcPerRow = 1;
	int color = world_rank  / numProcPerRow; // Determine color based on row

	// Split the communicator based on the color and use the
	// original rank for ordering
	MPI_Comm row_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);

	int numParallelComm = world_size / numProcPerRow;

	//int imageDimX = 101, imageDimY = 101, imageDimZ = 101;
	int imageDimX = 11, imageDimY = 11, imageDimZ = 11;
	//double *rho = malloc( imageDimX * imageDimY * imageDimZ * sizeof(double));
	//double *rho = Pt_test;
	//double *rho = CO_test_60;
	//double *rho = CO_test_90;
	//double *rho = CO_test_30;
	double *rho = O_test;
	int ii, imageSize = imageDimX * imageDimY * imageDimZ;
	//double current = 0.0;
	//for (ii = 0; ii < imageSize; ii++)
	//{
	//	rho[ii] = current;
	//	current += 1.0;
	//}

	double hx = 1.0, hy = 1.0, hz = 1.0;
	//double Uvec[9] = {0.0, 0.707106781186548, 0.707106781186548, 0.707106781186548, 0.0, 0.707106781186548, 0.707106781186548, 0.707106781186548, 0.0};
	double Uvec[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	double *U = Uvec;
	int accuracy = 2;

	int MCSHMaxOrder = 0;
	double MCSHMaxR = 2.0;
	//double MCSHRStepsize = 0.1;


	// printf("\nstart MCSH\n");
	// MCSHDescriptorMain_RadialRStep(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, U, accuracy, MCSHMaxOrder, MCSHMaxR, MCSHRStepsize, color, numParallelComm, row_comm);
	MCSHDescriptorMain_RadialLegendre(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, U, accuracy, MCSHMaxR, MCSHMaxOrder, 5, color, numParallelComm, row_comm);

	double end = MPI_Wtime();
	printf("\n\n END execution time %10f\n\n", end - start);

	MPI_Finalize();
	// printf("after MPI finalize\n");

	return 0;
}
/*
int main_old(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
	int world_rank, world_size;

	double start = MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int numProcPerRow = 1;
	int color = world_rank  / numProcPerRow; // Determine color based on row

	// Split the communicator based on the color and use the
	// original rank for ordering
	MPI_Comm row_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);

	int numParallelComm = world_size / numProcPerRow;

	int imageDimX = 60, imageDimY = 60, imageDimZ = 60;
	double *rho = malloc( imageDimX * imageDimY * imageDimZ * sizeof(double));
	int ii, imageSize = imageDimX * imageDimY * imageDimZ;
	double current = 0.0;
	for (ii = 0; ii < imageSize; ii++)
	{
		rho[ii] = current;
		current += 1.0;
	}

	double hx = 0.1, hy = 0.1, hz = 0.1;
	double Uvec[9] = {1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0};
	double *U = Uvec;
	int accuracy = 6;

	int MCSHMaxOrder = 3;
	double MCSHMaxR = 0.6;
	double MCSHRStepsize = 0.1;

	// printf("\nstart MCSH\n");
	// MCSHDescriptorMain_RadialRStep(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, U, accuracy, MCSHMaxOrder, MCSHMaxR, MCSHRStepsize, color, numParallelComm, row_comm);
	MCSHDescriptorMain_RadialLegendre(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, U, accuracy, 0.5, MCSHMaxOrder, 5, color, numParallelComm, row_comm);



	double end = MPI_Wtime();
	printf("\n\n END execution time %10f\n\n", end - start);

	MPI_Finalize();
	// printf("after MPI finalize\n");

	return 0;
}


int main_set(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
	int world_rank, world_size;

	double start = MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int numProcPerRow = 1;
	int color = world_rank  / numProcPerRow; // Determine color based on row

	// Split the communicator based on the color and use the
	// original rank for ordering
	MPI_Comm row_comm;
	MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);

	int numParallelComm = world_size / numProcPerRow;

	int imageDimX = 60, imageDimY = 60, imageDimZ = 60;
	double *rho = malloc( imageDimX * imageDimY * imageDimZ * sizeof(double));
	int ii, imageSize = imageDimX * imageDimY * imageDimZ;
	double current = 0.0;
	for (ii = 0; ii < imageSize; ii++)
	{
		rho[ii] = current;
		current += 1.0;
	}

	double hx = 0.1, hy = 0.1, hz = 0.1;
	double Uvec[9] = {1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0};
	double *U = Uvec;
	int accuracy = 6;

	MCSHDescriptorMainFixed_RadialRStep(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, U, accuracy, color, numParallelComm, row_comm);



	double end = MPI_Wtime();
	printf("\n\n END execution time %10f\n\n", end - start);

	MPI_Finalize();
	// printf("after MPI finalize\n");

	return 0;
}



int main_single(int argc, char *argv[])
{
	//printf("start\n");
	//double Uvec[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	//double *U = Uvec;
	//calcStencil(0.1, 0.1, 0.1, 0.5, 1, "100", U, 6, 0);


	MPI_Init(&argc,&argv);
	int rank, numProc;

	MPI_Comm_size(MPI_COMM_WORLD, &numProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	printf("rank %d out of %d\n", rank, numProc);


	int imageDimX = 100, imageDimY = 100, imageDimZ = 100;
	double *rho = malloc( imageDimX * imageDimY * imageDimZ * sizeof(double));

	int ii, imageSize = imageDimX * imageDimY * imageDimZ;
	double current = 0.0;
	for (ii = 0; ii < imageSize; ii++)
	{
		rho[ii] = current;
		current += 1.0;
	}


	double hx = 0.1, hy = 0.1, hz = 0.1;
	// double Uvec[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	// double Uvec[9] = {1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0};
	double Uvec[9] = {1.0, -0.292371704722737, 0.17364817766693, 0.0, 0.956304755963035, -0.110492654830881, 0.0, 0.0, 0.978589640054184};
	double *U = Uvec;
	double normalizedU[9];
	normalizeU(U, normalizedU);
	int accuracy = 6;

	printf("\n hx: %f \t hy: %f \t hz: %f \nU: \n %f \t %f \t %f \n %f \t %f \t %f \n %f \t %f \t %f \n",
			 hx, hy, hz, normalizedU[0], normalizedU[3], normalizedU[6], normalizedU[1], normalizedU[4], normalizedU[7], normalizedU[2], normalizedU[5], normalizedU[8] );
	// int length = 10;
	// double rCutoffList[10] = {0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5};
	// int lList[10] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
	// char *nList[] = {"000", "100", "000", "100", "000", "100", "000", "100", "000", "100"};

	// int length = 7;
	// double rCutoffList[7] = {0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5};
	// int lList[7] = {0, 1,2,2,3,3,3};
	// int groupList[7] = {1,1,1,2,1,2,3};


	int length = 2;
	double rCutoffList[2] = {0.5, 0.5};
	int lList[2] = {0, 1};
	char *nList[] = {"000", "100"};
	int groupList[2] = {1,1};

	// char **nListPointer = nList;

	int i;
	int rankCount;
	double *featureVector = malloc( imageSize * sizeof(double));
	for ( i = 0; i < imageSize; i++)
	{
		featureVector[i] = 0.0;
	}
	for (i = 0; i < length; i++){
		if (rank == i % numProc){
			printf("rank %d, index %d\n", rank, i);
			//calcStencil(const double hx, const double hx, const double hz, const double rCutoff, const int l, const char *n, const double *U, const int accuracy);
			// calcStencilAndSave(hx, hy, hz, rCutoffList[i], lList[i], nList[i], U, accuracy,i);

			// calcStencilAndConvolveAndSave(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoffList[i], lList[i], nList[i], U, accuracy,i);
			prepareMCSHFeatureAndSave(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoffList[i], lList[i], groupList[i], 1, 0, normalizedU, accuracy);
			// calcStencilAndConvolveAndAddResult(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoffList[i], lList[i], nList[i], normalizedU, accuracy, featureVector);

		}
	}

	MPI_Finalize();

	return 0;

}



*/




// int main(int argc, char *argv[])
// {


// 	int imageDimX = 100, imageDimY = 100, imageDimZ = 100;
// 	double *rho = malloc( imageDimX * imageDimY * imageDimZ * sizeof(double));

// 	int ii, imageSize = imageDimX * imageDimY * imageDimZ;
// 	double current = 0.0;
// 	for (ii = 0; ii < imageSize; ii++)
// 	{
// 		rho[ii] = current;
// 		current += 1.0;
// 	}


// 	double hx = 0.1, hy = 0.1, hz = 0.1;
// 	double Uvec[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
// 	double *U = Uvec;
// 	int accuracy = 6;

// 	int length = 10;
// 	double rCutoffList[10] = {0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5};
// 	int lList[10] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
// 	char *nList[] = {"000", "100", "000", "100", "000", "100", "000", "100", "000", "100"};
// 	// char **nListPointer = nList;

// 	int i;
// 	int rankCount;
// 	for (i = 0; i < length; i++){
// 		calcStencilAndConvolve(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoffList[i], lList[i], nList[i], U, accuracy,i);
// 	}


// }
