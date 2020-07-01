void MaxwellCartesianSphericalHarmonics(const double *x, const double *y, const double *z, const int l, const char *n, const double rCutoff, double *result, const int size);

void prepareMCSHFeatureAndSave(const double *image, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
							   const double rCutoff, const int l, const int group, const double *U, const int accuracy);

double scoreTask(const double rCutoff, const int l, const int group);

void calculateStencil(const int stencilDimX, const int stencilDimY, const int stencilDimZ, const double hx, const double hy, const double hz, 
					  const double rCutoff, const int l, const char *n, const double *U, const int accuracy, double *stencil);

void calcStencilAndSave(const double hx, const double hy, const double hz, const double rCutoff, const int l, const char *n, const double *U, const int accuracy, const int stencilIndex);

void calcStencilAndConvolveAndSave(const double *image, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
								   const double rCutoff, const int l, const char *n, const double *U, const int accuracy, const int stencilIndex);

void calcStencilAndConvolveAndAddResult(const double *image, const int imageDimX, const int imageDimY, const int imageDimZ, const double hx, const double hy, const double hz, 
										const double rCutoff, const int l, const char *n, const double *U, const int accuracy, double *convolveResult);