void taskPartition_RadialLegendre(const int length, const double rCutoff, const int *lList, const int numParallelComm, int *taskAssignmentList);

void taskPartition_RadialRStep(const int length, const double *rCutoffList, const int *lList, const int numParallelComm, int *taskAssignmentList);

void MCSHDescriptorMainFixed_RadialRStep(const double *rho, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
						double *U, const int accuracy, const int commIndex, const int numParallelComm, const MPI_Comm communicator);

void MCSHDescriptorMain_RadialLegendre(const double *rho, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
										double *U, const int accuracy, const double rCutoff, const int maxMCSHOrder, const int maxLegendreOrder,  
										const int commIndex, const int numParallelComm, const MPI_Comm communicator);

void MCSHDescriptorMain_RadialRStep(const double *rho, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
									double *U, const int accuracy, const int maxMCSHOrder, const double rMaxCutoff, const double rStepsize, 
									const int commIndex, const int numParallelComm, const MPI_Comm communicator);

void calcAndSaveCoords(const double *rho, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, double *U);