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









void taskPartition(const int length, const double *rCutoffList, const int *lList, const int *groupList, const int numParallelComm, int *taskAssignmentList)
{	
	// printf("before calc total score\n");
	double totalScore = 0;
	int i;
	for (i = 0; i < length; i++)
	{
		totalScore += scoreTask(rCutoffList[i], lList[i], groupList[i]);
		// averageGroupLoad = (totalScore / numParallelComm) + 1;
	}

	// printf("after calc total score, total score: %d\n", totalScore);

	double currentTotal = 0, currentRemaindingScore = totalScore;
	int currentRemaindingComm = numParallelComm, currentComm = 0;
	double currentTargetScore = currentRemaindingScore/currentRemaindingComm;

	// printf("step -1, currentTotal: %d, currentRemaindingScore: %d, currentRemaindingComm: %d, currentComm: %d, currentTargetScore: %d \n", currentTotal, currentRemaindingScore, currentRemaindingComm, currentComm, currentTargetScore);	
	
	for (i = 0; i < length; i++)
	{
		taskAssignmentList[i] = currentComm;
		currentTotal += scoreTask(rCutoffList[i], lList[i], groupList[i]);
		// printf("current total: %d\n", currentTotal);
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
		// printf("step %d, currentTotal: %d, currentRemaindingScore: %d, currentRemaindingComm: %d, currentComm: %d, currentTargetScore: %d \n", i, currentTotal, currentRemaindingScore, currentRemaindingComm, currentComm, currentTargetScore);
	}

	// printf("assigning task score\n");

	return;
}

void MCSHDescriptorMain(const double *rho, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
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
		taskPartition(length, rCutoffList, lList, groupList, numParallelComm, taskAssignmentList);

		// printf("after assigning task\n");

		int i;
		for (i = 0; i < length; i++)
		{
			if (rank == 0 && commIndex == taskAssignmentList[i])
			{
				prepareMCSHFeatureAndSave(rho, imageDimX, imageDimY, imageDimZ, hx, hy, hz, rCutoffList[i], lList[i], groupList[i], normalizedU, accuracy);
			}
		}

		free(taskAssignmentList);
	}

}
