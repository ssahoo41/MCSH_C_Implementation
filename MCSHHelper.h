typedef struct Point {
   double x;
   double y;
   double z;
} Point;

void die(const char *message);

void GetPlaneEquation(const Point p1, const Point p2, const Point p3, double *a, double *b, double *c, double *d);

int CheckPlaneIntersectWithSphere(const Point p1, const Point p2, const Point p3, const double rCutoff, const Point origin);

Point UTransform(const double x, const double y, const double z, const double *U);

void GetDimensionsPlane(const double hx, const double hy, const double hz, const double rCutoff, const double *U, int *dimXResult, int *dimYResult, int *dimZResult);

void printArr(double *arr, int size);

void sqrtVector(const double *arr, double *result, const int size);

void powVector(const double *arr, const int power, double *result, const int size);

void addVector(const double *arr1, const double *arr2, double *result, const int size);

void subtractVector(const double *arr1, const double *arr2, double *result, const int size);

void multiplyVector(const double *arr1, const double *arr2, double *result, const int size);

void divideVector(const double *arr1, const double *arr2, double *result, const int size);

void addScalarVector(const double *arr, const double a, double *result, const int size);

void multiplyScalarVector(const double *arr1, const double a, double *result, const int size);

void getRArray(const double *x, const double *y, const double *z, double *result, const int size);

void polyXYZArray(const double *x, const double *y, const double *z, const int powX, const int powY, const int powZ, const double a, double *result, const int size);

void polyArray(const double *x, const int powX, const double a, double *result, const int size);

void applyU(double *X, double *Y, double*Z, const double *U, const int size);

void applyU2(double *X, double *Y, double*Z, const double *U, const int size);

double calcDv(const double hx, const double hy, const double hz, const int accuracy, const double *U);

double sumArr(const double *arr, const int size);

double sumAbsArr(const double *arr, const int size);

void linspace(double start, double end, double *result, int num);

void meshgrid3D(const double *x, const double *y, const double *z, const int sizex, const int sizey, const int sizez, double *X, double *Y, double *Z);

void getCentralCoords(const double hx, const double hy, const double hz, const int accuracy, double *refX, double *refY, double *refZ);

double calcNorm3(double x1, double x2, double x3);

void normalizeU(double *U, double *normalizedU);

void calcSingleConvolveStep(const double *image, const double stencilCoeff, const int shiftX, const int shiftY, const int shiftZ, double *result, const int imageSize, const int imageDimX, const int imageDimY, const int imageDimZ);

void convolve(const double *image, const double *stencil, const int imageDimX, const int imageDimY, const int imageDimZ, const int stencilDimX, const int stencilDimY, const int stencilDimZ, double *result);

void convolve2(const double *image, const double *stencil, const int imageDimX, const int imageDimY, const int imageDimZ, const int stencilDimX, const int stencilDimY, const int stencilDimZ, double *result);

void convolve3(const double *image, const double *stencil, const int imageDimX, const int imageDimY, const int imageDimZ, const int stencilDimX, const int stencilDimY, const int stencilDimZ, double *result);

void writeMatToFile(const char *filename, const double *data, const int dimX, const int dimY, const int dimZ);

int mod(int a, int b);

void getMainParameter_RadialLegendre(const int maxMCSHOrder, const int maxLegendreOrder, const int length, int* LegendreOrderList, int* lList, int* groupList);

int getDescriptorListLength_RadialLegendre(const int maxLegendreOrder, const int maxMCSHOrder);

void getMainParameter_RadialRStep(const double rStepsize, const double rMaxCutoff, const int maxOrder, const int length, double* rCutoffList, int* lList, int* groupList);

int getDescriptorListLength_RadialRStep(const double rStepsize, const double rMaxCutoff, const int maxOrder);

int getNumRCutoff(const double rStepsize, const double rMaxCutoff);

int getNumGroup(const int maxOrder);

int getCurrentGroupNumber(const int currentIndex);

int getCurrentLNumber(const int currentIndex);
