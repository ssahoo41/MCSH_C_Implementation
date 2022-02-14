# include <stdio.h>
# include <stdlib.h>
# include <string.h>
#define RHO_FILE "image_0.dat"

int main(){
	int imageDimX = 101, imageDimY = 101, imageDimZ = 101;
	double *rho = malloc( imageDimX * imageDimY * imageDimZ * sizeof(double));
	int count = 0;
	int i = 0;
	char *eptr;
	FILE *file;

	file = fopen(RHO_FILE, "r");
	if (!file){
		perror("Error opening file");
        	return -1;
	}
	memset(rho, 0, sizeof(rho));
	while (!feof(file) && (count < imageDimX * imageDimY * imageDimZ))
	{
		
		fscanf(file, "%s", &(rho[count++]));
	}
	fclose(file);
	printf("Prices count: %d\n", count);
    	for(i = 0; i < count; i++)
    	{
        	printf("rho[%d] = %lf\n", i, rho[i]);
        }
	return 0;
}

