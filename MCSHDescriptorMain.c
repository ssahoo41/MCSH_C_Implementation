# include <stdio.h>
# include <stdlib.h>
# include <string.h>
#include <assert.h>
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
#include "MCSH.h"
#include "MCSHDescriptorMain.h"




void taskPartition_RadialLegendre(const int length, const double rCutoff, const int *lList, const int *groupList, const int numParallelComm, int *taskAssignmentList)
{	
	// printf("before calc total score\n");
	double totalScore = 0;
	int i;
	for (i = 0; i < length; i++)
	{	
		// printf("\nrCutoff: %f \t l: %d \t group: %d",rCutoffList[i], lList[i], groupList[i]);
		totalScore += scoreTask(rCutoff, lList[i], groupList[i]);
		// averageGroupLoad = (totalScore / numParallelComm) + 1;
	}

	// printf("after calc total score, total score: %f\n", totalScore);

	double currentTotal = 0, currentRemaindingScore = totalScore;
	int currentRemaindingComm = numParallelComm, currentComm = 0;
	double currentTargetScore = currentRemaindingScore/currentRemaindingComm;

	// printf("step -1, currentTotal: %f, currentRemaindingScore: %f, currentRemaindingComm: %d, currentComm: %d, currentTargetScore: %f \n", currentTotal, currentRemaindingScore, currentRemaindingComm, currentComm, currentTargetScore);	
	
	for (i = 0; i < length; i++)
	{
		taskAssignmentList[i] = currentComm;
		currentTotal += scoreTask(rCutoff, lList[i], groupList[i]);
		// printf("current total: %f\n", currentTotal);
		if (currentTotal >= currentTargetScore )
		{
			currentComm++;
			currentRemaindingScore -= currentTotal; 
			currentRemaindingComm--; 
			currentTotal = 0;
			if (currentRemaindingScore != 0)
			{
				currentTargetScore = currentRemaindingScore/currentRemaindingComm;
			}
			
		}
		// printf("step %d, currentTotal: %f, currentRemaindingScore: %f, currentRemaindingComm: %d, currentComm: %d, currentTargetScore: %f \n", i, currentTotal, currentRemaindingScore, currentRemaindingComm, currentComm, currentTargetScore);
	}

	// printf("assigning task score\n");

	return;
}




void taskPartition_RadialRStep(const int length, const double *rCutoffList, const int *lList, const int *groupList, const int numParallelComm, int *taskAssignmentList)
{	
	// printf("before calc total score\n");
	double totalScore = 0;
	int i;
	for (i = 0; i < length; i++)
	{	
		// printf("\nrCutoff: %f \t l: %d \t group: %d",rCutoffList[i], lList[i], groupList[i]);
		totalScore += scoreTask(rCutoffList[i], lList[i], groupList[i]);
		// averageGroupLoad = (totalScore / numParallelComm) + 1;
	}

	// printf("after calc total score, total score: %f\n", totalScore);

	double currentTotal = 0, currentRemaindingScore = totalScore;
	int currentRemaindingComm = numParallelComm, currentComm = 0;
	double currentTargetScore = currentRemaindingScore/currentRemaindingComm;

	// printf("step -1, currentTotal: %f, currentRemaindingScore: %f, currentRemaindingComm: %d, currentComm: %d, currentTargetScore: %f \n", currentTotal, currentRemaindingScore, currentRemaindingComm, currentComm, currentTargetScore);	
	
	for (i = 0; i < length; i++)
	{
		taskAssignmentList[i] = currentComm;
		currentTotal += scoreTask(rCutoffList[i], lList[i], groupList[i]);
		// printf("current total: %f\n", currentTotal);
		if (currentTotal >= currentTargetScore )
		{
			currentComm++;
			currentRemaindingScore -= currentTotal; 
			currentRemaindingComm--; 
			currentTotal = 0;
			if (currentRemaindingScore != 0)
			{
				currentTargetScore = currentRemaindingScore/currentRemaindingComm;
			}
			
		}
		// printf("step %d, currentTotal: %f, currentRemaindingScore: %f, currentRemaindingComm: %d, currentComm: %d, currentTargetScore: %f \n", i, currentTotal, currentRemaindingScore, currentRemaindingComm, currentComm, currentTargetScore);
	}

	// printf("assigning task score\n");

	return;
}

void MCSHDescriptorMain_RadialLegendre(const double *rho, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
										double *U, const int accuracy, const double rCutoff, const int maxMCSHOrder, const int maxLegendreOrder,  
										const int commIndex, const int numParallelComm, const MPI_Comm communicator)
{
	int rank, numProc;

	MPI_Comm_size(communicator, &numProc);
	MPI_Comm_rank(communicator, &rank);
	// printf("start length");

	int length = getDescriptorListLength_RadialLegendre(maxLegendreOrder, maxMCSHOrder);
	// printf("\nlength: %d\n", length);
	int *LegendreOrderList = calloc( length, sizeof(int));
	int *lList = calloc( length, sizeof(int));
	int *groupList = calloc( length, sizeof(int));

	// int LegendreOrderList[70] = {0,1,2,3,4,5,6,7,8,9,
	// 								0,1,2,3,4,5,6,7,8,9,
	// 								0,1,2,3,4,5,6,7,8,9,
	// 								0,1,2,3,4,5,6,7,8,9,
	// 								0,1,2,3,4,5,6,7,8,9,
	// 								0,1,2,3,4,5,6,7,8,9,
	// 								0,1,2,3,4,5,6,7,8,9};

	// int lList[70] = {0,0,0,0,0,0,0,0,0,0,
	// 				 1,1,1,1,1,1,1,1,1,1,
	// 				 2,2,2,2,2,2,2,2,2,2,
	// 				 2,2,2,2,2,2,2,2,2,2,
	// 				 3,3,3,3,3,3,3,3,3,3,
	// 				 3,3,3,3,3,3,3,3,3,3,
	// 				 3,3,3,3,3,3,3,3,3,3};

	// int groupList[70] = {1,1,1,1,1,1,1,1,1,1,
	// 					 1,1,1,1,1,1,1,1,1,1,
	// 					 1,1,1,1,1,1,1,1,1,1,
	// 					 2,2,2,2,2,2,2,2,2,2,
	// 					 1,1,1,1,1,1,1,1,1,1,
	// 					 2,2,2,2,2,2,2,2,2,2,
	// 					 3,3,3,3,3,3,3,3,3,3};
	//void getMainParameter_RadialLegendre(const int maxMCSHOrder, const int maxLegendreOrder, const int length, int* LegendreOrderList, int* lList, int* groupList)
	getMainParameter_RadialLegendre(maxMCSHOrder, maxLegendreOrder, length, LegendreOrderList, lList, groupList);

	// int counter;
	// for (counter = 0; counter < length; counter ++)
	// {
	// 	printf("\ni:%d \t order:%d \t group:%d \t r:%.10f", counter, lList[counter], groupList[counter], rCutoffList[counter]);
	// }


	if (rank == 0)
	{
		double normalizedU[9];
		normalizeU(U, normalizedU);

		int *taskAssignmentList = calloc( length, sizeof(int));
		taskPartition_RadialLegendre(length, rCutoff, lList, groupList, numParallelComm, taskAssignmentList);

		int i;
		for (i = 0; i < length; i++)
		{
			// printf("\n %.5f \t %d \t %d \t %d \t %d\n ",rCutoffList[i], lList[i], groupList[i], taskAssignmentList[i], commIndex);
			if (commIndex == taskAssignmentList[i])
			{	
				//radial type: 2 for Legendre Polynomial
				prepareMCSHFeatureAndSave(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoff, lList[i], groupList[i], 2, LegendreOrderList[i], normalizedU, accuracy);
			}
		}

		free(taskAssignmentList);
	}

}


void MCSHDescriptorMain_RadialRStep(const double *rho, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
									double *U, const int accuracy, const int maxMCSHOrder, const double rMaxCutoff, const double rStepsize, 
									const int commIndex, const int numParallelComm, const MPI_Comm communicator)
{
	int rank, numProc;

	MPI_Comm_size(communicator, &numProc);
	MPI_Comm_rank(communicator, &rank);
	// printf("start length");
	int length = getDescriptorListLength_RadialRStep(rStepsize, rMaxCutoff, maxMCSHOrder);
	// printf("\nlength: %d\n", length);
	double *rCutoffList = calloc( length, sizeof(double));
	int *lList = calloc( length, sizeof(int));
	int *groupList = calloc( length, sizeof(int));

	getMainParameter_RadialRStep(rStepsize, rMaxCutoff, maxMCSHOrder, length, rCutoffList, lList, groupList);

	// int counter;
	// for (counter = 0; counter < length; counter ++)
	// {
	// 	printf("\ni:%d \t order:%d \t group:%d \t r:%.10f", counter, lList[counter], groupList[counter], rCutoffList[counter]);
	// }


	if (rank == 0)
	{
		double normalizedU[9];
		normalizeU(U, normalizedU);

		int *taskAssignmentList = calloc( length, sizeof(int));
		taskPartition_RadialRStep(length, rCutoffList, lList, groupList, numParallelComm, taskAssignmentList);

		int i;
		for (i = 0; i < length; i++)
		{
			// printf("\n %.5f \t %d \t %d \t %d \t %d\n ",rCutoffList[i], lList[i], groupList[i], taskAssignmentList[i], commIndex);
			if (commIndex == taskAssignmentList[i])
			{
				// radial type 1: r step
				// radial order: 0 (default, not relevant)
				prepareMCSHFeatureAndSave(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoffList[i], lList[i], groupList[i], 1, 0, normalizedU, accuracy);
			}
		}

		free(taskAssignmentList);
	}

}



void MCSHDescriptorMainFixed_RadialRStep(const double *rho, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
						double *U, const int accuracy, const int commIndex, const int numParallelComm, const MPI_Comm communicator)
{
	int rank, numProc;

	MPI_Comm_size(communicator, &numProc);
	MPI_Comm_rank(communicator, &rank);

	// int length = 10;
	// double rCutoffList[10] = {0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5};
	// int lList[10] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
	// int groupList[10] = {1,1,1,1,1,1,1,1,1,1};

	int length = 70;

	double rCutoffList[70] = {	0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
								0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
								0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
								0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
								0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
								0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
								0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};

	int lList[70] = {0,0,0,0,0,0,0,0,0,0,
					 1,1,1,1,1,1,1,1,1,1,
					 2,2,2,2,2,2,2,2,2,2,
					 2,2,2,2,2,2,2,2,2,2,
					 3,3,3,3,3,3,3,3,3,3,
					 3,3,3,3,3,3,3,3,3,3,
					 3,3,3,3,3,3,3,3,3,3};

	int groupList[70] = {1,1,1,1,1,1,1,1,1,1,
						 1,1,1,1,1,1,1,1,1,1,
						 1,1,1,1,1,1,1,1,1,1,
						 2,2,2,2,2,2,2,2,2,2,
						 1,1,1,1,1,1,1,1,1,1,
						 2,2,2,2,2,2,2,2,2,2,
						 3,3,3,3,3,3,3,3,3,3};

	if (rank == 0)
	{
		double normalizedU[9];
		normalizeU(U, normalizedU);

		// int length = 2;
		// double rCutoffList[2] = {0.3, 0.3};
		// int lList[2] = {0, 1};
		// int groupList[2] = {1,1};

		

		// printf("before assigning task, num processors: %d\n", numParallelComm);

		int *taskAssignmentList = malloc( length * sizeof(int));
		taskPartition_RadialRStep(length, rCutoffList, lList, groupList, numParallelComm, taskAssignmentList);

		// printf("after assigning task\n");

		int i;
		for (i = 0; i < length; i++)
		{
			if (rank == 0 && commIndex == taskAssignmentList[i])
			{
				prepareMCSHFeatureAndSave(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoffList[i], lList[i], groupList[i], 1, 0, normalizedU, accuracy);
			}
		}

		free(taskAssignmentList);
	}

}

void calcAndSaveCoords(const double *rho, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, double *U)
{
	int imageSize = imageDimX * imageDimY * imageDimZ;

	double *refX = calloc( imageDimX, sizeof(double));
	double *refY = calloc( imageDimY, sizeof(double));
	double *refZ = calloc( imageDimZ, sizeof(double));

	linspace(0.5 * hx, (imageDimX - 0.5) * hx, refX, imageDimX);
	linspace(0.5 * hy, (imageDimY - 0.5) * hy, refY, imageDimY);
	linspace(0.5 * hz, (imageDimZ - 0.5) * hz, refZ, imageDimZ);

	// int ii;
	// printf("\n Nx: %d \t Ny: %d \t Nz: %d \t \t hx: %.10f \t hy: %.10f \t hz: %.10f \n", imageDimX, imageDimY, imageDimZ, hx, hy, hz);
	// printf("\nX linspace: \t");
	// for (ii = 0; ii < imageDimX; ii++)
	// {
	// 	printf("%.10f\t", refX[ii]);
	// }

	// printf("\nY linspace: \t");
	// for (ii = 0; ii < imageDimY; ii++)
	// {
	// 	printf("%.10f\t", refY[ii]);
	// }

	// printf("\nZ linspace: \t");
	// for (ii = 0; ii < imageDimZ; ii++)
	// {
	// 	printf("%.10f\t", refZ[ii]);
	// }

	double *X = calloc( imageSize, sizeof(double));
	double *Y = calloc( imageSize, sizeof(double));
	double *Z = calloc( imageSize, sizeof(double));

	meshgrid3D(refX, refY, refZ, imageDimX, imageDimY, imageDimZ, X, Y, Z);
	applyU2(X, Y, Z, U, imageSize);

	char filename[128] = "Mesh_Coordinates.csv";
	//snprintf(stencilFilename, 128, "%s_%d.csv", "test_output", stencilIndex);

	FILE *output_fp = fopen(filename,"w");
	if (output_fp == NULL) {
		printf("\nCannot open file \"%s\"\n",filename);
		die("cant open file");
	}

	int i,j,k, index = 0;
	for (k = 0; k < imageDimZ; k++){
		for ( j = 0; j < imageDimY; j++) {
			for ( i = 0; i < imageDimX; i++) {
				fprintf(output_fp,"%d,%d,%d,%.15f,%.15f,%.15f,%.15f\n",i,j,k,X[index],Y[index],Z[index],rho[index]);
				index ++;
			}
		}
	}


	fclose(output_fp);

	free(refX);
	free(refY);
	free(refZ);

	free(X);
	free(Y);
	free(Z);
}