#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include "tools_mlff.h"
#include "spherical_harmonics.h"
#include "soap_descriptor.h"
#include "mlff_types.h"
#include "ddbp_tools.h"
#include "tools.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

/*
Prints descriptor for a given atom
*/
void print_descriptors(double *X2, double *X3, double **dX2_dX, double **dX2_dY, double **dX2_dZ, double **dX2_dF, double **dX3_dF, double **dX3_dX, double **dX3_dY, double **dX3_dZ, int size_X2, int size_X3, int neighs){
	FILE *fptr;
	fptr = fopen("descriptors.txt","w");
	if(fptr == NULL)
   {
      printf("Error opening the descriptor.txt file!");   
      exit(1);             
   }
   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"%s\n","              X2:            ");
   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"\n");
   for (int i = 0; i < size_X2; i++){
   	fprintf(fptr, "%10.9f", X2[i]);
   	fprintf(fptr, "\n");
   }

   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"%s\n","              X2 derivative with atom positions:            ");
   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"\n");
   for (int i = 0; i < size_X2; i++){
   	for (int j = 0; j < 1+neighs; j++){
   		fprintf(fptr, "%10.9f %10.9f %10.9f ", dX2_dX[j][i], dX2_dY[j][i], dX2_dZ[j][i]);
   	}
   	fprintf(fptr, "\n");
   }

   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"%s\n","              X2 derivative with deformation gradient:            ");
   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"\n");
   for (int i = 0; i < size_X2; i++){
   	for (int j = 0; j < 6; j++){
   		fprintf(fptr, "%10.9f ", dX2_dF[j][i]);
   	}
   	fprintf(fptr, "\n");
   }


   fprintf(fptr,"\n");
   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"%s\n","              X3:            ");
   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"\n");
   for (int i = 0; i < size_X3; i++){
   	fprintf(fptr, "%10.9f ", X3[i]);
   	fprintf(fptr, "\n");
   }
   fclose(fptr);

   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"%s\n","              X3 derivative with atom positions:            ");
   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"\n");
   for (int i = 0; i < size_X3; i++){
   	for (int j = 0; j < 1+neighs; j++){
   		fprintf(fptr, "%10.9f %10.9f %10.9f ", dX3_dX[j][i], dX3_dY[j][i], dX3_dZ[j][i]);
   	}
   	fprintf(fptr, "\n");
   }

   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"%s\n","              X3 derivative with deformation gradient:            ");
   fprintf(fptr,"%s\n","********************************************************");
   fprintf(fptr,"\n");
   for (int i = 0; i < size_X3; i++){
   	for (int j = 0; j < 6; j++){
   		fprintf(fptr, "%10.9f ", dX3_dF[j][i]);
   	}
   	fprintf(fptr, "\n");
   }
}

void print_cnlm(double complex **cnlm, double complex ***dcnlm_dX, double complex ***dcnlm_dY, double complex ***dcnlm_dZ, double complex ***dcnlm_dF, int size_cnlm, int* neighs, int nelem){
	FILE *fptr;
	fptr = fopen("cnlm.txt","w");
	if(fptr == NULL)
   {
      printf("Error opening the cnlm.txt file!");   
      exit(1);             
   }
   fprintf(fptr,"cnlm:\n");

   

   for (int j = 1; j < nelem; j++){
   	for (int i = 0; i < size_cnlm; i++){
   		// fprintf(fptr, "%10.9f+%10.9fi  %10.9f+%10.9fi", creal(cnlm[j][i]),cimag(cnlm[j][i]),
   		// 		creal(dcnlm_dX[j][0][i]), cimag(dcnlm_dX[j][0][i]));
   		fprintf(fptr, " %10.9f+%10.9fi",
   				creal(dcnlm_dY[j][0][i]), cimag(dcnlm_dY[j][0][i]));

   		// for (int k = 0; k < 1+neighs[j]; k++){
   		// 	fprintf(fptr, "%10.9f+%10.9fi %10.9f+%10.9fi %10.9f+%10.9fi ", creal(dcnlm_dX[j][k][i]), cimag(dcnlm_dX[j][k][i]),
   		// 	 creal(dcnlm_dY[j][k][i]), cimag(dcnlm_dY[j][k][i]), creal(dcnlm_dZ[j][k][i]), cimag(dcnlm_dZ[j][k][i]));
   		// }
   		// for (int k = 0; k < 6; k++){
   		// 	fprintf(fptr, "%10.9f+%10.9fi %10.9f+%10.9fi %10.9f+%10.9fi ", creal(dcnlm_dF[j][k][i]),cimag(dcnlm_dF[j][k][i]));
   		// }
   		fprintf(fptr, "\n");
   	}
   }
   fclose(fptr);

}


// *****************Later might replace this function with C version of hnl generator. *************************************
void read_h_nl(const int N, const int L, double *rgrid, double *h_nl, double *dh_nl){
	int i, j, info;
	FILE *fp;
	char line[512];
	char a1[512], a2[512], a3[512], a4[512];
	int count1=0, count2=0;

	fp = fopen("hnl.txt","r");

	fgets(line, sizeof (line), fp);
	sscanf(line, "%s%s%s%s", a1, a2, a3, a4);
	int N_r = atoi(a4);

	for (i = 0; i < N_r; i++){
		info = fscanf(fp, "%s", a1);
		rgrid[i] = atof(a1);

	}

	for (i =0; i < N*(L+1); i++){
		info = fscanf(fp, "%s %s", a1, a2);

		for (j = 0; j < N_r; j++){
			info = fscanf(fp, "%s", a1);
			h_nl[count1] = atof(a1);
			count1++;	
		}
		for (j=0; j<N_r; j++){
			info = fscanf(fp, "%s", a1);
			dh_nl[count2] = atof(a1);
			count2++;
		}
	}


	fclose(fp);
}
// *****************Later might replace this function with C version of hnl generator. *************************************




/*
initialize_nlist function initializes the a objects in the NeighList structure and also allocates memory for the dunamic arrays.

[Input]
1. natom: Number of atoms in the system
2. rcut: Cutoff distance (bohr)
3. nelem: Number of element species 
[Output]
1. nlist: pointer to Neighlist structure to be initialized
*/


void initialize_nlist(NeighList* nlist, const int natom, const double rcut, const int nelem ){
	int i, j;
 	nlist->natom = natom;
 	nlist->nelem = nelem;
	nlist->rcut = rcut;
	nlist->natom_elem = (int *) malloc(nelem*sizeof(int));
	for (i=0; i<nelem; i++){
		nlist->natom_elem[i] = 0;
	}
	nlist->Nneighbors = (int *) malloc(natom*sizeof(int));
	nlist->unique_Nneighbors = (int *) malloc(natom*sizeof(int));
	nlist->Nneighbors_elemWise = (int **) malloc(natom*sizeof(int*));
	nlist->unique_Nneighbors_elemWise = (int **) malloc(natom*sizeof(int*));
	for (i=0; i<natom; i++){
		nlist->Nneighbors[i] = 0;
		nlist->unique_Nneighbors[i] = 0;
		nlist->Nneighbors_elemWise[i] = (int *) malloc(nelem*sizeof(int));
		nlist->unique_Nneighbors_elemWise[i] = (int *) malloc(nelem*sizeof(int));
		for (j = 0; j<nelem; j++){
			nlist->Nneighbors_elemWise[i][j] = 0;
			nlist->unique_Nneighbors_elemWise[i][j] = 0;
		}
	}
	
	nlist->neighborList = (dyArray *) malloc(sizeof(dyArray)*natom);
	nlist->unique_neighborList = (dyArray *) malloc(sizeof(dyArray)*natom);
	nlist->neighborList_elemWise = (dyArray **) malloc(sizeof(dyArray*)*natom);
	nlist->unique_neighborList_elemWise = (dyArray **) malloc(sizeof(dyArray*)*natom);
	nlist->neighborAtmTyp = (dyArray *) malloc(sizeof(dyArray)*natom);
	nlist->neighborList_imgX = (dyArray *) malloc(sizeof(dyArray)*natom);
	nlist->neighborList_imgY = (dyArray *) malloc(sizeof(dyArray)*natom);
	nlist->neighborList_imgZ = (dyArray *) malloc(sizeof(dyArray)*natom);
	nlist->neighborList_elemWise_imgX = (dyArray **) malloc(sizeof(dyArray*)*natom);
	nlist->neighborList_elemWise_imgY = (dyArray **) malloc(sizeof(dyArray*)*natom);
	nlist->neighborList_elemWise_imgZ = (dyArray **) malloc(sizeof(dyArray*)*natom);
	nlist->localID_neighbours = (dyArray *) malloc(sizeof(dyArray)*natom);
	nlist->localID_neighbours_elem = (dyArray **) malloc(sizeof(dyArray*)*natom);


	
	for (i =0; i < natom; i++){
		init_dyarray(nlist->neighborList + i);
		init_dyarray(nlist->unique_neighborList + i);
		init_dyarray(nlist->neighborAtmTyp + i);
		init_dyarray(nlist->neighborList_imgX + i);
		init_dyarray(nlist->neighborList_imgY + i);
		init_dyarray(nlist->neighborList_imgZ + i);
		init_dyarray(nlist->localID_neighbours + i);


		nlist->neighborList_elemWise[i] = (dyArray *) malloc(sizeof(dyArray)*nelem);
		nlist->unique_neighborList_elemWise[i] = (dyArray *) malloc(sizeof(dyArray)*nelem);
		nlist->neighborList_elemWise_imgX[i] = (dyArray *) malloc(sizeof(dyArray)*nelem);
		nlist->neighborList_elemWise_imgY[i] = (dyArray *) malloc(sizeof(dyArray)*nelem);
		nlist->neighborList_elemWise_imgZ[i] = (dyArray *) malloc(sizeof(dyArray)*nelem);
		nlist->localID_neighbours_elem[i] = (dyArray *) malloc(sizeof(dyArray)*nelem);

		for (j=0; j<nelem; j++){
			init_dyarray(&(nlist->neighborList_elemWise[i][j]));
			init_dyarray(&(nlist->unique_neighborList_elemWise[i][j]));
			init_dyarray(&(nlist->neighborList_elemWise_imgX[i][j]));
			init_dyarray(&(nlist->neighborList_elemWise_imgY[i][j]));
			init_dyarray(&(nlist->neighborList_elemWise_imgZ[i][j]));
			init_dyarray(&(nlist->localID_neighbours_elem[i][j]));

		}
	}
}

/*
lin_search function searches for the first occurence of an element in the array 

[Input]
1. arr: pointer to array
2. n: length of array
3. x: element to search

[Output]
1. i: ID of the element in the array
*/

int lin_search(int *arr, int n, int x)
{
    int i;
    for (i = 0; i < n; i++)
        if (arr[i] == x)
            return i;
    return -1;
}

/*
build_nlist function calculate the list of neighbours for each atoms in the structure.

[Input]
1. rcut: Cutoff distance (bohr) 
2. nelem: Number of element species 
3. natom: Number of atoms in the system
4. atompos: Pointer to the array containing the atom positions of the atoms in the fundamental cell (stored RowMajor format)
5. atomtyp: Pointer to the array containing the element type of the atoms in the fundamental cell (stored RowMajor format)
6. BC[]: Boundary condition (0 for Dirchelet, 1 for Periodic) [Not used currently, add later]
7. cell: 3*1 array to store the length of fundamental cell. [Currently only orthogonal cells are implemented]

[Output]
1. nlist: pointer to Neighlist structure
*/

void build_nlist(const double rcut, const int nelem, const int natom, const double * const atompos,
				 int * atomtyp, int BC[], double* cell, NeighList* nlist) {

	int img_px, img_nx, img_py, img_ny, img_pz, img_nz, i, j, img_x, img_y, img_z, k, ID, idx_neigh;
	double xi, yi, zi, xj, yj, zj, dx, dy, dz, dr, L1, L2, L3;
	nlist->cell[0] = cell[0];
	nlist->cell[1] = cell[1];
	nlist->cell[2] = cell[2];
	L1 = cell[0];
	L2 = cell[1];
	L3 = cell[2];
    //Why does this get called here and in sparc_interface.c?
	initialize_nlist(nlist, natom, rcut, nelem);

	for (int i=0; i<natom; i++){
		nlist->natom_elem[atomtyp[i]] = nlist->natom_elem[atomtyp[i]] +1;
	}

	for (i = 0; i < natom; i++){
		xi = atompos[3*i];
		yi = atompos[3*i+1];
		zi = atompos[3*i+2];
		img_px = max(0,ceil((rcut - L1 + xi) /L1));
		img_nx = max(0,ceil((rcut - xi) /L1));
		img_py = max(0,ceil((rcut - L2 + yi) /L2));
		img_ny = max(0,ceil((rcut - yi) /L2));
		img_pz = max(0,ceil((rcut - L3 + zi) /L3));
		img_nz = max(0,ceil((rcut - zi) /L3));
		// count_neighs = 0;
		for (j = 0; j < natom; j++){
			int count_unique=0;
			xj = atompos[3*j];
			yj = atompos[3*j+1];
			zj = atompos[3*j+2];
			for (img_x = -img_nx; img_x <= img_px; img_x++){
				if (abs(xj + img_x*L1 - xi) >= rcut) continue;		
				for (img_y = -img_ny; img_y <= img_py; img_y++){
					if (abs(yj + img_y*L2 - yi) >= rcut) continue;		
					for (img_z = -img_nz; img_z <= img_pz; img_z++){
						if (abs(zj + img_z*L3 - zi) >= rcut) continue;
						if ((i==j) && (img_x==0) && (img_y==0) && (img_z==0)) continue;
						
						dx = xj + img_x*L1 - xi;
						dy = yj + img_y*L2 - yi;
						dz = zj + img_z*L3 - zi;
						dr = sqrt(dx*dx + dy*dy + dz*dz);
						if (dr >= rcut) continue;

						nlist->Nneighbors[i] += 1;
						nlist->Nneighbors_elemWise[i][atomtyp[j]] += 1;

						
						if (count_unique==0) {
							nlist->unique_Nneighbors[i] += 1;
							nlist->unique_Nneighbors_elemWise[i][atomtyp[j]] += 1;
							append_dyarray(nlist->unique_neighborList + i, j);
							append_dyarray(&(nlist->unique_neighborList_elemWise[i][atomtyp[j]]), j);	
							// count_neighs++;			
						}
						// append_dyarray(nlist->localID_neighbours + i, count_neighs-1);
	

						append_dyarray(&(nlist->neighborList_elemWise_imgX[i][atomtyp[j]]), img_x);
						append_dyarray(&(nlist->neighborList_elemWise_imgY[i][atomtyp[j]]), img_y);
						append_dyarray(&(nlist->neighborList_elemWise_imgZ[i][atomtyp[j]]), img_z);

						append_dyarray(nlist->neighborList + i, j);
						append_dyarray(&(nlist->neighborList_elemWise[i][atomtyp[j]]), j);



						append_dyarray(nlist->neighborList_imgX + i, img_x);
						append_dyarray(nlist->neighborList_imgY + i, img_y);
						append_dyarray(nlist->neighborList_imgZ + i, img_z);

						append_dyarray(nlist->neighborAtmTyp+i, atomtyp[j]);

						count_unique++;
					}
				}
			} 	
		}
	}
	
	for (i = 0; i < natom; i++) {
		for (j = 0; j < nelem; j++){
			for (k = 0; k < nlist->Nneighbors_elemWise[i][j]; k++){
				idx_neigh =  nlist->neighborList_elemWise[i][j].array[k];
				ID = lin_search((nlist->unique_neighborList_elemWise[i][j]).array, nlist->unique_Nneighbors_elemWise[i][j], idx_neigh);
				append_dyarray(&(nlist->localID_neighbours[i]), ID);
			}
		}
	}

	for (i = 0; i < natom; i++) {
		for (j = 0; j < nlist->unique_Nneighbors[i]; j++) {
			idx_neigh =  nlist->unique_neighborList[i].array[j];
			for (k = 0; k < nelem; k++){
				ID = lin_search((nlist->unique_neighborList_elemWise[i][k]).array, nlist->unique_Nneighbors_elemWise[i][k], idx_neigh);
				append_dyarray(&(nlist->localID_neighbours_elem[i][k]), ID);
			}
		}
	}

}

/*
clear_nlist function frees the memory dynamically allocated for Neighlist 

[Input]
1. nlist: pointer to Neighlist structure

[Output]
1. nlist: pointer to Neighlist structure
*/

void clear_nlist(NeighList* nlist) {
	int i, j;
	free(nlist->Nneighbors);
	free(nlist->unique_Nneighbors);
	for (i=0; i < nlist->natom; i++){
		free(nlist->Nneighbors_elemWise[i]);
		free(nlist->unique_Nneighbors_elemWise[i]);
		delete_dyarray(nlist->neighborList+i);
		delete_dyarray(nlist->neighborList+i);
		delete_dyarray(nlist->localID_neighbours+i);
		for (j=0; j < nlist->nelem; j++){
			delete_dyarray(&(nlist->neighborList_elemWise[i][j]));
			delete_dyarray(&(nlist->unique_neighborList_elemWise[i][j]));
			delete_dyarray(&(nlist->neighborList_elemWise_imgX[i][j]));
			delete_dyarray(&(nlist->neighborList_elemWise_imgY[i][j]));
			delete_dyarray(&(nlist->neighborList_elemWise_imgZ[i][j]));
			delete_dyarray(&(nlist->localID_neighbours_elem[i][j]));
		}
		free(nlist->neighborList_elemWise[i]);
		free(nlist->neighborList_elemWise_imgX[i]);
		free(nlist->neighborList_elemWise_imgY[i]);
		free(nlist->neighborList_elemWise_imgZ[i]);
		free(nlist->localID_neighbours_elem[i]);
	}
	free(nlist->localID_neighbours);
	free(nlist->localID_neighbours_elem);
	free(nlist->neighborList_elemWise);
	free(nlist->unique_neighborList_elemWise);
	free(nlist->neighborList);
	free(nlist->unique_neighborList);
	free(nlist->neighborAtmTyp);
	free(nlist->neighborList_imgX);
	free(nlist->neighborList_imgY);
	free(nlist->neighborList_imgZ);
	free(nlist->neighborList_elemWise_imgX);
	free(nlist->neighborList_elemWise_imgY);
	free(nlist->neighborList_elemWise_imgZ);
}


/*
uniqueEle function finds number of unique entries in an array

[Input]
1. a: pointer to the array
2. n: size of array

[Output]
1. count: Number of unique entries
*/

int uniqueEle(int* a, int n)      //Function Definition
{
   int i, j, count = 1;
   //Traverse the array
   for (i = 1; i < n; i++)      //hold an array element
   {
      for (j = 0; j < i; j++)   
      {
         if (a[i] == a[j])    //Check for duplicate elements
         {
            break;             //If duplicate elements found then break
         }
      }
      if (i == j)
      {
         count++;     //increment the number of distinct elements
      }
   }
   return count;      //Return the number of distinct elements
}


/*
initialize_soapObj function initializes the objects in SoapObj structure and also allocates memory of dynamic arrays

[Input]
1. nlist: pointer to Neighlist structure
2. Lmax: maximum angular quantum number of the spherical harmonics considered in SOAP
3. Nmax: number of radial basis functions considered in SOAP (spherical Bessel function is used in this)
4. beta_3: Weight to the three body term in SOAP
5. xi_3: exponent in the dot product kernel of SOAP
[Output]
1. soap_str: pointer to SoapObj structure
*/

void initialize_soapObj(SoapObj *soap_str, NeighList *nlist, int Lmax, int Nmax, double beta_3, double xi_3){

	int size_X2, size_X3, Nmod, i, j, k;
	int nelem = nlist->nelem, natom = nlist->natom;
	Nmod = nelem * Nmax;
	// size_cnlm = (Lmax+1)*(Lmax+1)*Nmax;
	size_X3 = ((Nmod*(1 + Nmod))/2)*(Lmax+1);
	size_X2 = Nmax*nelem;
	soap_str->size_X2 = size_X2;
	soap_str->size_X3 = size_X3;
	soap_str->beta_2 = 1 - beta_3;
	soap_str->beta_3 = beta_3;
	soap_str->xi_3= xi_3;

	soap_str->natom = natom;
	soap_str->rcut = nlist->rcut;
	soap_str->Lmax = Lmax;
	soap_str->Nmax = Nmax;
	soap_str->nelem = nelem;
	soap_str->cell[0] = nlist->cell[0];
	soap_str->cell[1] = nlist->cell[1];
	soap_str->cell[2] = nlist->cell[2];
	soap_str->N_rgrid = 72;

	soap_str->Nneighbors = (int *) malloc(sizeof(int)*natom);
	for (i=0; i < natom; i++){
		soap_str->Nneighbors[i] = nlist->Nneighbors[i];
	}

	soap_str->unique_Nneighbors = (int *) malloc(sizeof(int)*natom);
	for (i=0; i < natom; i++){
		soap_str->unique_Nneighbors[i] = nlist->unique_Nneighbors[i];
	}

	soap_str->unique_Nneighbors_elemWise = (int **) malloc(natom*sizeof(int*));
	for (i=0; i < natom; i++){
		soap_str->unique_Nneighbors_elemWise[i] = (int *) malloc(nelem*sizeof(int));
	}
	for (i=0; i < natom; i++){
		for(j=0; j < nelem; j++){
			soap_str->unique_Nneighbors_elemWise[i][j] = nlist->unique_Nneighbors_elemWise[i][j];
		}
	}

	soap_str->unique_neighborList_elemWise = (dyArray **) malloc(sizeof(dyArray*)*natom);
	for (i =0; i < natom; i++){
		soap_str->unique_neighborList_elemWise[i] = (dyArray *) malloc(sizeof(dyArray)*nelem);
	}

	for (i =0; i < natom; i++){
		for(j=0; j < nelem; j++){
			soap_str->unique_neighborList_elemWise[i][j].len = nlist->unique_neighborList_elemWise[i][j].len;
			soap_str->unique_neighborList_elemWise[i][j].capacity = nlist->unique_neighborList_elemWise[i][j].capacity;
			soap_str->unique_neighborList_elemWise[i][j].array = (int *)malloc(sizeof(int)*nlist->unique_neighborList_elemWise[i][j].len);
			for (k =0; k < nlist->unique_neighborList_elemWise[i][j].len; k++){
				soap_str->unique_neighborList_elemWise[i][j].array[k] = nlist->unique_neighborList_elemWise[i][j].array[k];
			}
		}
	}

	soap_str->neighborList = (dyArray *) malloc(sizeof(dyArray)*natom);
	for (i =0; i < natom; i++){
		soap_str->neighborList[i].len = nlist->neighborList[i].len;
		soap_str->neighborList[i].capacity = nlist->neighborList[i].capacity;
		soap_str->neighborList[i].array = (int *)malloc(sizeof(int)*nlist->neighborList[i].len);
		for (k =0; k<nlist->neighborList[i].len; k++){
			soap_str->neighborList[i].array[k] = nlist->neighborList[i].array[k];
		}
	}

	soap_str->unique_neighborList = (dyArray *) malloc(sizeof(dyArray)*natom);
	for (i =0; i < natom; i++){
		soap_str->unique_neighborList[i].len = nlist->unique_neighborList[i].len;
		soap_str->unique_neighborList[i].capacity = nlist->unique_neighborList[i].capacity;
		soap_str->unique_neighborList[i].array = (int *)malloc(sizeof(int)*nlist->unique_neighborList[i].len);
		for (k =0; k<nlist->unique_neighborList[i].len; k++){
			soap_str->unique_neighborList[i].array[k] = nlist->unique_neighborList[i].array[k];
		}
	}


	soap_str->natom_elem = (int *) malloc(sizeof(int)*nelem);
	for (i=0; i<nelem; i++){
		soap_str->natom_elem[i] = nlist->natom_elem[i];
	}

	soap_str->X2 = (double **) malloc(natom * sizeof(double*));
	soap_str->X3 = (double **) malloc(natom * sizeof(double*));
	soap_str->dX2_dF = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX2_dX = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX2_dY = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX2_dZ = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX3_dF = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX3_dX = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX3_dY = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX3_dZ = (double ***) malloc(natom * sizeof(double**));

	for (i=0; i < natom; i++){
		soap_str->X2[i] = (double *) malloc(size_X2 * sizeof(double));
		soap_str->X3[i] = (double *) malloc(size_X3 * sizeof(double));
		for (int sz=0; sz < size_X2; sz++){
			soap_str->X2[i][sz] = 0.0;
		}
		for (int sz=0; sz < size_X3; sz++){
			soap_str->X3[i][sz] = 0.0;
		}

		int uniq_natms = uniqueEle((nlist->neighborList[i]).array, nlist->Nneighbors[i]);
		soap_str->dX2_dX[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));
		soap_str->dX2_dY[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));
		soap_str->dX2_dZ[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));

		soap_str->dX3_dX[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));
		soap_str->dX3_dY[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));
		soap_str->dX3_dZ[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));



		soap_str->dX2_dF[i] = (double **) malloc(6 * sizeof(double*));
		soap_str->dX3_dF[i] = (double **) malloc(6 * sizeof(double*));
		for (j=0; j < 1+uniq_natms; j++){
			soap_str->dX2_dX[i][j] = (double *) malloc(size_X2 * sizeof(double));
			soap_str->dX2_dY[i][j] = (double *) malloc(size_X2 * sizeof(double));
			soap_str->dX2_dZ[i][j] = (double *) malloc(size_X2 * sizeof(double));
			soap_str->dX3_dX[i][j] = (double *) malloc(size_X3 * sizeof(double));
			soap_str->dX3_dY[i][j] = (double *) malloc(size_X3 * sizeof(double));
			soap_str->dX3_dZ[i][j] = (double *) malloc(size_X3 * sizeof(double));
			for (int sz=0; sz < size_X2; sz++){
				soap_str->dX2_dX[i][j][sz] = 0.0;
				soap_str->dX2_dY[i][j][sz] = 0.0;
				soap_str->dX2_dZ[i][j][sz] = 0.0;
			}
			for (int sz=0; sz < size_X3; sz++){
				soap_str->dX3_dX[i][j][sz] = 0.0;
				soap_str->dX3_dY[i][j][sz] = 0.0;
				soap_str->dX3_dZ[i][j][sz] = 0.0;
			}
		}
		for (j=0; j < 6; j++){
			soap_str->dX2_dF[i][j] = (double *) malloc(size_X2 * sizeof(double));
			soap_str->dX3_dF[i][j] = (double *) malloc(size_X3 * sizeof(double));
			for (int sz=0; sz < size_X2; sz++){
				soap_str->dX2_dF[i][j][sz] = 0.0;
			}
			for (int sz=0; sz < size_X3; sz++){
				soap_str->dX3_dF[i][j][sz] = 0.0;
			}
		}
	}

}



/*
delete_soapObj function frees the memory dynamically allocated for SoapObj 

[Input]
1. soap_str: pointer to SoapObj structure

[Output]
1. soap_str: pointer to SoapObj structure
*/

void delete_soapObj(SoapObj *soap_str){
	int size_cnlm, size_X2, size_X3, Nmod, nelem, i, j;
	Nmod = soap_str->nelem * soap_str->Nmax;
	nelem = soap_str->nelem;
	// size_cnlm = (soap_str->Lmax+1)*(soap_str->Lmax+1)*soap_str->Nmax;
	size_X3 = ( ((Nmod)*(1 + Nmod))/2)*(soap_str->Lmax+1);
	size_X2 = soap_str->Nmax * soap_str->nelem;

	for (i=0; i<soap_str->natom; i++){
		free(soap_str->X2[i]);
		free(soap_str->X3[i]);
		free(soap_str->unique_Nneighbors_elemWise[i]);

		for (j=0; j < soap_str->nelem; j++){
			delete_dyarray(&(soap_str->unique_neighborList_elemWise[i][j]));
		}

		free(soap_str->unique_neighborList_elemWise[i]);
		int uniq_natms = uniqueEle((soap_str->neighborList[i]).array, soap_str->Nneighbors[i]);
		delete_dyarray(&(soap_str->neighborList[i]));
		delete_dyarray(&(soap_str->unique_neighborList[i]));
		for (j=0; j<1+uniq_natms; j++){
			free(soap_str->dX2_dX[i][j]);
			free(soap_str->dX2_dY[i][j]);
			free(soap_str->dX2_dZ[i][j]);
			free(soap_str->dX3_dX[i][j]);
			free(soap_str->dX3_dY[i][j]);
			free(soap_str->dX3_dZ[i][j]);
		}
		for (j=0; j<6; j++){
			free(soap_str->dX2_dF[i][j]);
			free(soap_str->dX3_dF[i][j]);
		}
		free(soap_str->dX2_dX[i]);
		free(soap_str->dX2_dY[i]);
		free(soap_str->dX2_dZ[i]);
		free(soap_str->dX3_dX[i]);
		free(soap_str->dX3_dY[i]);
		free(soap_str->dX3_dZ[i]);
		free(soap_str->dX2_dF[i]);
		free(soap_str->dX3_dF[i]);
	}

	free(soap_str->unique_Nneighbors);
	free(soap_str->Nneighbors);
	free(soap_str->natom_elem);
	free(soap_str->neighborList);
	free(soap_str->unique_neighborList);
	free(soap_str->unique_neighborList_elemWise);

	free(soap_str->X2);
	free(soap_str->X3);
	free(soap_str->dX2_dX);
	free(soap_str->dX2_dY);
	free(soap_str->dX2_dZ);
	free(soap_str->dX3_dX);
	free(soap_str->dX3_dY);
	free(soap_str->dX3_dZ);
	free(soap_str->dX2_dF);
	free(soap_str->dX3_dF);
}




/*
build_soapObj function calculates the SOAP descriptors and their derivatives w.r.t the atom positions and the deformation gradient

[Input]
1. nlist: pointer to the NeighList structure
2. rgrid: pointer to the rgrid for the spline interpolation
3. h_nl: pointer to the h_nl function for spline interpolation
4. dh_nl: pointer to the derivative of h_nl for spline interpolation
5. atompos: pointer to the atom positions [stored in ColMajor]

[Output]
1. soap_str: pointer to the SoapObj structure
*/

void build_soapObj(SoapObj *soap_str, NeighList *nlist, double* rgrid, double* h_nl, double* dh_nl, double *atompos, int Nmax, int Lmax, double beta_3, double xi_3) {


	const double PI = 3.141592653589793;
	int i, j, k, N_r;
	double xi, yi, zi, xj, yj, zj;
	double dx,dy,dz,dr,dtheta,dphi;
	double L1, L2, L3;
	double complex *Y, *dY_theta, *dY_phi;
	double *YD_hnl, *YD_dhnl;
	double cth, sth, cph, sph, dx2dy2z, dx2dy2;
	double dtheta_dx, dtheta_dy, dtheta_dz, dphi_dx, dphi_dy, dphi_dz;
	double dr_dF11, dr_dF22, dr_dF33, dr_dF12, dr_dF23, dr_dF13;
	double dth_dF11, dth_dF22, dth_dF33, dth_dF12, dth_dF23, dth_dF13;
	double dphi_dF11, dphi_dF22, dphi_dF33, dphi_dF12, dphi_dF23, dphi_dF13, dxdr, dydr, dzdr;
	int idx1, idx2, idx3, idx_neigh, elem_typ, n, l, m;
	double *hnl_temp, *dhnl_temp, *rtemp;
	int *ntemp, *imgx_temp, *imgy_temp, *imgz_temp, *elemtyp_temp;
	double complex ***cnlm, ****dcnlm_dX, ****dcnlm_dY, ****dcnlm_dZ, ****dcnlm_dF;
	int natom= nlist->natom, nelem = nlist->nelem, uniq_el_atms, local_index;
	int size_cnlm =  (Lmax + 1) * (Lmax + 2) * Nmax;
	size_cnlm = size_cnlm/2;

	initialize_soapObj(soap_str, nlist, Lmax, Nmax, beta_3, xi_3);
	cnlm = (double complex***) malloc(natom*sizeof(double complex**));
	dcnlm_dX = (double complex****) malloc(natom*sizeof(double complex***));
	dcnlm_dY = (double complex****) malloc(natom*sizeof(double complex***));
	dcnlm_dZ = (double complex****) malloc(natom*sizeof(double complex***));
	dcnlm_dF = (double complex****) malloc(natom*sizeof(double complex***));

	for (i = 0; i < natom; i++){
		cnlm[i] = (double complex**) malloc(nelem*sizeof(double complex*));
		for (j = 0; j < nelem; j++){
			cnlm[i][j] = (double complex*) malloc(size_cnlm*sizeof(double complex));
			for (k = 0; k < size_cnlm; k++){
				cnlm[i][j][k] = 0.0 + 0.0*I;
			}
		}
		dcnlm_dX[i] = (double complex***) malloc(nelem*sizeof(double complex**));
		dcnlm_dY[i] = (double complex***) malloc(nelem*sizeof(double complex**));
		dcnlm_dZ[i] = (double complex***) malloc(nelem*sizeof(double complex**));
		dcnlm_dF[i] = (double complex***) malloc(nelem*sizeof(double complex**));

		for (j = 0; j < nelem; j++){
			uniq_el_atms = nlist->unique_Nneighbors_elemWise[i][j];
			dcnlm_dX[i][j] = (double complex**) malloc((1+uniq_el_atms)*sizeof(double complex*));
			dcnlm_dY[i][j] = (double complex**) malloc((1+uniq_el_atms)*sizeof(double complex*));
			dcnlm_dZ[i][j] = (double complex**) malloc((1+uniq_el_atms)*sizeof(double complex*));
			dcnlm_dF[i][j] = (double complex**) malloc(6*sizeof(double complex*));
			for (k = 0; k < 1+uniq_el_atms; k++){
				dcnlm_dX[i][j][k] = (double complex*) malloc(size_cnlm*sizeof(double complex));
				dcnlm_dY[i][j][k] = (double complex*) malloc(size_cnlm*sizeof(double complex));
				dcnlm_dZ[i][j][k] = (double complex*) malloc(size_cnlm*sizeof(double complex));
				for (l = 0; l < size_cnlm; l++){
					dcnlm_dX[i][j][k][l] = 0.0 + 0.0*I;
					dcnlm_dY[i][j][k][l] = 0.0 + 0.0*I;
					dcnlm_dZ[i][j][k][l] = 0.0 + 0.0*I;
				}
			}
			for (k = 0; k < 6; k++){
				dcnlm_dF[i][j][k] = (double complex*) malloc(size_cnlm*sizeof(double complex));
				for (l = 0; l < size_cnlm; l++){
					dcnlm_dF[i][j][k][l] = 0.0 + 0.0*I;
				}
			}
		}
	}


	N_r = soap_str->N_rgrid;
	L1 = soap_str->cell[0];
	L2 = soap_str->cell[1];
	L3 = soap_str->cell[2];
	Lmax = soap_str->Lmax;
	Nmax = soap_str->Nmax;



	Y = (double complex *) malloc((Lmax+1)*(Lmax+1)*sizeof(double complex));
	dY_theta = (double complex *) malloc((Lmax+1)*(Lmax+1)*sizeof(double complex));
	dY_phi = (double complex *) malloc((Lmax+1)*(Lmax+1)*sizeof(double complex));
	hnl_temp = (double *) malloc(sizeof(double)*Nmax*(Lmax+1));
	dhnl_temp = (double *) malloc(sizeof(double)*Nmax*(Lmax+1));
	rtemp = (double *) malloc(sizeof(double)*1);

	YD_hnl = (double *) malloc(N_r*(Nmax)*(Lmax+1)*sizeof(double));
	YD_dhnl = (double *) malloc(N_r*(Nmax)*(Lmax+1)*sizeof(double));

	// derivatives of tabulated h_nl, dh_nl required for spline interpolation (used later in the inner loop).
	for (i =0; i < Nmax; i++){
		for (j =0; j < Lmax+1; j++){
			getYD_gen(rgrid, h_nl + (i*(Lmax+1)+j)*N_r, YD_hnl + (i*(Lmax+1)+j)*N_r, N_r);
			getYD_gen(rgrid, dh_nl + (i*(Lmax+1)+j)*N_r, YD_dhnl + (i*(Lmax+1)+j)*N_r, N_r);
		}
	}





	for (i = 0; i < soap_str->natom; i++){
		xi = atompos[3*i];
		yi = atompos[3*i+1];
		zi = atompos[3*i+2];
		ntemp = (nlist->neighborList +i)->array;
		imgx_temp=(nlist->neighborList_imgX + i)->array;
		imgy_temp=(nlist->neighborList_imgY + i)->array;
		imgz_temp=(nlist->neighborList_imgZ + i)->array;
		elemtyp_temp=(nlist->neighborAtmTyp +i)->array;

		for (j = 0; j < nlist->Nneighbors[i]; j++){
			idx_neigh = ntemp[j];
			elem_typ = elemtyp_temp[j];
			xj = atompos[3*idx_neigh] + L1 * imgx_temp[j];
			yj = atompos[3*idx_neigh+1] + L2 * imgy_temp[j];
			zj = atompos[3*idx_neigh+2] + L3 * imgz_temp[j];

			dx = xj - xi;
			dy = yj - yi;
			dz = zj - zi;

			dr = sqrt(dx*dx + dy*dy + dz*dz);

			dtheta = atan(sqrt(dx*dx + dy*dy)/dz);
			if (dz==0) dtheta = PI/2;
			if (dtheta<0) dtheta = PI + dtheta;

			dphi = atan(dy/dx);
			if (dx>0 && dy>0) dphi = dphi;  // 1st Quandrant
			if (dx<0 && dy>0) dphi = PI+dphi;  // 2nd Quandrant
			if (dx<0 && dy<0) dphi = PI+dphi;  // 3rd Quadrant
			if (dx>0 && dy<0) dphi = 2*PI+dphi;  // 4th Quadrant
			if (dx==0 && dy >0) dphi = PI/2;
			if (dx==0 && dy <0) dphi = 3*PI/2;
			if (dx==0 && dy ==0) dphi = 0;


			cth = cos(dtheta);
			sth = sin (dtheta);

			dx2dy2z = dz*sqrt(dx*dx+dy*dy);
			dx2dy2 = sqrt(dx*dx+dy*dy);

			
			if (dx2dy2z == 0){
				dtheta_dx = 0.0;
				dtheta_dy = 0.0;
			} else {
				dtheta_dx = (dx*cth*cth)/(dx2dy2z);
				dtheta_dy = (dy*cth*cth)/(dx2dy2z);
			}

			if (dz == 0){
				dtheta_dz = 0.0;
			} else {
				dtheta_dz = -1*sth*cth / dz;
			}


			cph = cos(dphi);
			sph = sin (dphi);

			dphi_dz = 0;
			if (dx == 0){
				dphi_dx = 0.0;
				dphi_dy = 0.0;
			} else {
				dphi_dx = (-1*sph*cph)/dx;
				dphi_dy = (cph*cph)/dx;
			}

			dr_dF11 = dx*dx/dr;
			dr_dF22 = dy*dy/dr;
			dr_dF33 = dz*dz/dr;
			dr_dF12 = 2*dx*dy/dr;
			dr_dF23 = 2*dy*dz/dr;
			dr_dF13 = 2*dz*dx/dr;

			
			if (dx2dy2z==0){
				dth_dF11=0.0;
				dth_dF22=0.0;
				dth_dF12=0.0;
				dth_dF23=0.0;
				dth_dF13=0.0;

			} else {
				dth_dF11 = (dx*dx*cth*cth)/(dx2dy2z);
				dth_dF22 = (dy*dy*cth*cth)/(dx2dy2z);
				dth_dF12 = (2*dx*dy*cth*cth)/(dx2dy2z);
				dth_dF23 = (dy*cth*cth)/dx2dy2 - dy*cth*sth/dz;
				dth_dF13 = (dx*cth*cth)/dx2dy2 - dx*cth*sth/dz;
			}
			dth_dF33 = -1*sth*cth;
			

			dphi_dF11 = -1*sph*cph;
			if (dx==0){
				dphi_dF22=0.0;
				dphi_dF12=0.0;
				dphi_dF23=0.0;
				dphi_dF13=0.0;

			} else {
				dphi_dF22 = (dy*cph*cph)/dx;
				dphi_dF12 = (dx*cph*cph - sph*cph*dy)/dx;
				dphi_dF23 = cph*cph*dz/dx;
				dphi_dF13 = -1*sph*cph*dz/dx;
			}
			
			dphi_dF33 = 0;

			dxdr = dx/dr, dydr = dy/dr, dzdr = dz/dr;

			sph_harmonics(dtheta, dphi, soap_str->Lmax, Y, dY_theta, dY_phi);

			rtemp[0] = dr;
			// Calculate h_nl and dh_nl from dpline interpoation for all combinations of n and l
			for (n=0; n < Nmax; n++){
				for (l=0; l < Lmax+1; l++){
					SplineInterp(rgrid, h_nl+(n*(Lmax+1)+l)*N_r,
						N_r, rtemp, hnl_temp + n*(Lmax+1)+l, 1, YD_hnl+(n*(Lmax+1)+l)*N_r);

					SplineInterp(rgrid, dh_nl+(n*(Lmax+1)+l)*N_r,
						N_r, rtemp, dhnl_temp + n*(Lmax+1)+l, 1, YD_dhnl+(n*(Lmax+1)+l)*N_r);
				}
			}



			// int local_index = lin_search((nlist->unique_neighborList_elemWise[i][elem_typ]).array,
			 							// nlist->unique_Nneighbors_elemWise[i][elem_typ], idx_neigh);
			local_index = 1+nlist->localID_neighbours[i].array[j];

			for (n=0; n < Nmax; n++){
				for (l=0; l < Lmax+1; l++){
					for (m=0; m < l+1; m++){
						idx1 = (n*(Lmax+1)*(Lmax+2))/2 + (l*(l+1))/2 + m;
						idx2 = n*(Lmax+1)+l;
						idx3 = l*l + l + m;

						cnlm[i][elem_typ][idx1] +=  hnl_temp[idx2] * conj(Y[idx3]);


						if (local_index != -1){

							dcnlm_dX[i][elem_typ][local_index][idx1] += dhnl_temp[idx2] * (dxdr)*conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dx + dY_phi[idx3]*dphi_dx);

							dcnlm_dY[i][elem_typ][local_index][idx1] += dhnl_temp[idx2] * (dydr)* conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dy + dY_phi[idx3]*dphi_dy);

							dcnlm_dZ[i][elem_typ][local_index][idx1] += dhnl_temp[idx2] * (dzdr)* conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dz + dY_phi[idx3]*dphi_dz);


							dcnlm_dX[i][elem_typ][0][idx1] -= dhnl_temp[idx2] * (dxdr)*conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dx + dY_phi[idx3]*dphi_dx);

							dcnlm_dY[i][elem_typ][0][idx1] -= dhnl_temp[idx2] * (dydr)* conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dy + dY_phi[idx3]*dphi_dy);

							dcnlm_dZ[i][elem_typ][0][idx1] -= dhnl_temp[idx2] * (dzdr)* conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dz + dY_phi[idx3]*dphi_dz);

							
						}
						
						
						dcnlm_dF[i][elem_typ][0][idx1] += dhnl_temp[idx2]*dr_dF11 * conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF11 + dY_phi[idx3]*dphi_dF11);

						dcnlm_dF[i][elem_typ][1][idx1] += dhnl_temp[idx2]*dr_dF22 * conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF22 + dY_phi[idx3]*dphi_dF22);

						dcnlm_dF[i][elem_typ][2][idx1] += dhnl_temp[idx2]*dr_dF33 * conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF33 + dY_phi[idx3]*dphi_dF33);

						dcnlm_dF[i][elem_typ][3][idx1] += dhnl_temp[idx2]*dr_dF12 * conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF12 + dY_phi[idx3]*dphi_dF12);

						dcnlm_dF[i][elem_typ][4][idx1] += dhnl_temp[idx2]*dr_dF23 * conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF23 + dY_phi[idx3]*dphi_dF23);

						dcnlm_dF[i][elem_typ][5][idx1] += dhnl_temp[idx2]*dr_dF13 * conj(Y[idx3]) +
									 hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF13 + dY_phi[idx3]*dphi_dF13);

					}
				}
			}
			
		}
	}

	double complex el1_xder, el1_yder, el1_zder, el2_xder, el2_yder, el2_zder,el_xder,el_yder,el_zder;
	int na, el, count_X2, neigh, count_X3;
	
	for (na = 0; na < natom; na++){
		 count_X2=0;
		for (el = 0; el < nelem; el++){
			for (n = 0; n < soap_str->Nmax; n++){
				idx1 = (n*(Lmax+1)*(Lmax+2))/2 ;

				soap_str->X2[na][count_X2] += (double) (cnlm[na][el][idx1]);
														 
				soap_str->dX2_dF[na][0][count_X2] += (double) (dcnlm_dF[na][el][0][idx1]);
														 
				soap_str->dX2_dF[na][1][count_X2] += (double) (dcnlm_dF[na][el][1][idx1]);
														 
				soap_str->dX2_dF[na][2][count_X2] += (double) (dcnlm_dF[na][el][2][idx1]);
														 
				soap_str->dX2_dF[na][3][count_X2] += (double) (dcnlm_dF[na][el][3][idx1]);
														 
				soap_str->dX2_dF[na][4][count_X2] += (double) (dcnlm_dF[na][el][4][idx1]);
														 
				soap_str->dX2_dF[na][5][count_X2] += (double) (dcnlm_dF[na][el][5][idx1]);
														 

				for (neigh = 0; neigh < nlist->unique_Nneighbors[na]; neigh++){

					// local_index  = lin_search((nlist->unique_neighborList_elemWise[na][el]).array,
 					// nlist->unique_Nneighbors_elemWise[na][el], (nlist->unique_neighborList[na]).array[neigh]);
					local_index = 1+nlist->localID_neighbours_elem[na][el].array[neigh];
			
					if (local_index == 0){
						el_xder = 0.0 + 0.0*I;
						el_yder = 0.0 + 0.0*I;
						el_zder = 0.0 + 0.0*I;
					} else{
						el_xder = dcnlm_dX[na][el][local_index][idx1];
						el_yder = dcnlm_dY[na][el][local_index][idx1];
						el_zder = dcnlm_dZ[na][el][local_index][idx1];
					}

					soap_str->dX2_dX[na][neigh+1][count_X2]  += (double)(el_xder);

					soap_str->dX2_dY[na][neigh+1][count_X2]  += (double)(el_yder);

					soap_str->dX2_dZ[na][neigh+1][count_X2]  +=  (double)(el_zder);
				}

				soap_str->dX2_dX[na][0][count_X2] +=  (double) dcnlm_dX[na][el][0][idx1];
				soap_str->dX2_dY[na][0][count_X2] +=  (double) dcnlm_dY[na][el][0][idx1];
				soap_str->dX2_dZ[na][0][count_X2] +=  (double) dcnlm_dZ[na][el][0][idx1];

				count_X2++;
			}
		}
	}


	
	int el1, n1, n2_el2, n2, el2, local_index2, local_index1;
	double const_X3, temp1, temp2, temp3;

	for (na = 0; na < natom; na++){
		count_X3=0;
		for (el1 = 0; el1 < nelem; el1++){
			for (n1 = 0; n1 < soap_str->Nmax; n1++){
				temp1 = (n1*(Lmax+1)*(Lmax+2))/2;
				for (n2_el2 = el1*soap_str->Nmax + n1; n2_el2 < nelem*soap_str->Nmax; n2_el2++){
					n2 = n2_el2 % soap_str->Nmax;
					el2 = n2_el2 / soap_str->Nmax;
					temp2 = (n2*(Lmax+1)*(Lmax+2))/2;
					for (l = 0; l < soap_str->Lmax+1; l++){
						const_X3 = sqrt(8*PI*PI/(2*l+1));
						temp3 = (l*(l+1))/2;
						for (m = 0; m < l+1; m++){
							idx1 = temp1 + temp3 + m;
							idx2 = temp2 + temp3 + m;
							if (m==0){
								soap_str->X3[na][count_X3] += const_X3 * (double) (cnlm[na][el1][idx1]*conj(cnlm[na][el2][idx2]));
							} else {
								soap_str->X3[na][count_X3] += const_X3 * (double) (2* creal(cnlm[na][el1][idx1]*conj(cnlm[na][el2][idx2])));
							}
							// The derivative w.r.t to itself
							el1_xder = dcnlm_dX[na][el1][0][idx1];
							el1_yder = dcnlm_dY[na][el1][0][idx1];
							el1_zder = dcnlm_dZ[na][el1][0][idx1];
							el2_xder = dcnlm_dX[na][el2][0][idx2];
							el2_yder = dcnlm_dY[na][el2][0][idx2];
							el2_zder = dcnlm_dZ[na][el2][0][idx2];
							if (m==0){
								soap_str->dX3_dX[na][0][count_X3]  += const_X3 * (double)( (el1_xder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_xder)));

								soap_str->dX3_dY[na][0][count_X3]  += const_X3 * (double)( (el1_yder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_yder)));

								soap_str->dX3_dZ[na][0][count_X3]  += const_X3 * (double)( (el1_zder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_zder)));
							} else {
								soap_str->dX3_dX[na][0][count_X3]  += const_X3 * 2*(double)creal( (el1_xder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_xder)));

								soap_str->dX3_dY[na][0][count_X3]  += const_X3 * 2*(double)creal( (el1_yder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_yder)));

								soap_str->dX3_dZ[na][0][count_X3]  += const_X3 * 2*(double)creal( (el1_zder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_zder)));
							}
							// The derivative w.r.t to itself end
							for (neigh = 0; neigh < nlist->unique_Nneighbors[na]; neigh++){
								local_index1 = 1+nlist->localID_neighbours_elem[na][el1].array[neigh];
								local_index2 = 1+nlist->localID_neighbours_elem[na][el2].array[neigh];
								if (local_index1 == 0){
									el1_xder = 0.0 + 0.0*I;
									el1_yder = 0.0 + 0.0*I;
									el1_zder = 0.0 + 0.0*I;
								} else{
									el1_xder = dcnlm_dX[na][el1][local_index1][idx1];
									el1_yder = dcnlm_dY[na][el1][local_index1][idx1];
									el1_zder = dcnlm_dZ[na][el1][local_index1][idx1];
								}
								if (local_index2 == 0){
									el2_xder = 0.0 + 0.0*I;
									el2_yder = 0.0 + 0.0*I;
									el2_zder = 0.0 + 0.0*I;
								} else{
									el2_xder = dcnlm_dX[na][el2][local_index2][idx2];
									el2_yder = dcnlm_dY[na][el2][local_index2][idx2];
									el2_zder = dcnlm_dZ[na][el2][local_index2][idx2];
								}

 
								if (m==0){
									soap_str->dX3_dX[na][1+neigh][count_X3]  += const_X3 * (double)( (el1_xder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_xder)));

									soap_str->dX3_dY[na][1+neigh][count_X3]  += const_X3 * (double)( (el1_yder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_yder)));

									soap_str->dX3_dZ[na][1+neigh][count_X3]  += const_X3 * (double)( (el1_zder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_zder)));
								} else {
									soap_str->dX3_dX[na][1+neigh][count_X3]  += const_X3 * 2*(double)creal( (el1_xder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_xder)));

									soap_str->dX3_dY[na][1+neigh][count_X3]  += const_X3 * 2*(double)creal( (el1_yder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_yder)));

									soap_str->dX3_dZ[na][1+neigh][count_X3]  += const_X3 * 2*(double)creal( (el1_zder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_zder)));
								}
							}
							if (m==0){
								soap_str->dX3_dF[na][0][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][0][idx1]*conj(cnlm[na][el2][idx2])) +
															 (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][0][idx2])));

								soap_str->dX3_dF[na][1][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][1][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][1][idx2])));

								soap_str->dX3_dF[na][2][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][2][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][2][idx2])));

								soap_str->dX3_dF[na][3][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][3][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][3][idx2])));

								soap_str->dX3_dF[na][4][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][4][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][4][idx2])));

								soap_str->dX3_dF[na][5][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][5][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][5][idx2])));
							} else {
								soap_str->dX3_dF[na][0][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][0][idx1]*conj(cnlm[na][el2][idx2])) +
															 (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][0][idx2])));

								soap_str->dX3_dF[na][1][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][1][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][1][idx2])));

								soap_str->dX3_dF[na][2][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][2][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][2][idx2])));

								soap_str->dX3_dF[na][3][count_X3] += const_X3 * 2* (double)creal( (dcnlm_dF[na][el1][3][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][3][idx2])));

								soap_str->dX3_dF[na][4][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][4][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][4][idx2])));

								soap_str->dX3_dF[na][5][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][5][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][5][idx2])));
							}
						}
						count_X3++;
					}
				}
			}
			
		}

	}
	// if (1)
	// 	print_cnlm(cnlm[0], dcnlm_dX[0], dcnlm_dY[0], dcnlm_dZ[0], dcnlm_dF[0], size_cnlm, soap_str->unique_Nneighbors_elemWise[0], nelem);
	// free the memory for temperory Y, dY, hnl etc.
	free(Y);
	free(dY_theta);
	free(dY_phi);
	free(YD_hnl);
	free(YD_dhnl);
	free(hnl_temp);
	free(dhnl_temp);
	free(rtemp);

	// free the memory for cnlm
	for (i = 0; i < natom; i++){
		for (j = 0; j < nelem; j++){
			free(cnlm[i][j]);
		}
		free(cnlm[i]);

		for (j = 0; j < nelem; j++){
			uniq_el_atms = nlist->unique_Nneighbors_elemWise[i][j];
			for (k = 0; k < 1+uniq_el_atms; k++){
				free(dcnlm_dX[i][j][k]); free(dcnlm_dY[i][j][k]); free(dcnlm_dZ[i][j][k]);
			}
			for (k = 0; k < 6; k++){
				free(dcnlm_dF[i][j][k]);
			}
			free(dcnlm_dX[i][j]); free(dcnlm_dY[i][j]); free(dcnlm_dZ[i][j]); free(dcnlm_dF[i][j]); 
		}
		free(dcnlm_dX[i]); free(dcnlm_dY[i]); free(dcnlm_dZ[i]); free(dcnlm_dF[i]);
	}
	free(dcnlm_dX); free(dcnlm_dY); free(dcnlm_dZ); free(dcnlm_dF); free(cnlm);

}


/*
initialize_soapObj_wZ function initializes the objects in SoapObj structure and also allocates memory of dynamic arrays

[Input]
1. nlist: pointer to Neighlist structure
2. Lmax: maximum angular quantum number of the spherical harmonics considered in SOAP
3. Nmax: number of radial basis functions considered in SOAP (spherical Bessel function is used in this)
4. beta_3: Weight to the three body term in SOAP
5. xi_3: exponent in the dot product kernel of SOAP
[Output]
1. soap_str: pointer to SoapObj structure
*/

void initialize_soapObj_wZ(SoapObj *soap_str, NeighList *nlist, int Lmax, int Nmax, double beta_3, double xi_3){

	int size_X2, size_X3, Nmod, i, j, k;
	int nelem, natom = nlist->natom;
	nelem = 1;
	Nmod = nelem * Nmax;
	// size_cnlm = (Lmax+1)*(Lmax+1)*Nmax;
	size_X3 = ((Nmod*(1 + Nmod))/2)*(Lmax+1);
	size_X2 = Nmax*nelem;
	soap_str->size_X2 = size_X2;
	soap_str->size_X3 = size_X3;
	soap_str->beta_2 = 1 - beta_3;
	soap_str->beta_3 = beta_3;
	soap_str->xi_3= xi_3;

	soap_str->natom = natom;
	soap_str->rcut = nlist->rcut;
	soap_str->Lmax = Lmax;
	soap_str->Nmax = Nmax;
	soap_str->nelem = nlist->nelem;
	soap_str->cell[0] = nlist->cell[0];
	soap_str->cell[1] = nlist->cell[1];
	soap_str->cell[2] = nlist->cell[2];
	soap_str->N_rgrid = 72;

	soap_str->Nneighbors = (int *) malloc(sizeof(int)*natom);
	for (i=0; i < natom; i++){
		soap_str->Nneighbors[i] = nlist->Nneighbors[i];
	}

	soap_str->unique_Nneighbors = (int *) malloc(sizeof(int)*natom);
	for (i=0; i < natom; i++){
		soap_str->unique_Nneighbors[i] = nlist->unique_Nneighbors[i];
	}

	soap_str->unique_Nneighbors_elemWise = (int **) malloc(natom*sizeof(int*));
	for (i=0; i < natom; i++){
		soap_str->unique_Nneighbors_elemWise[i] = (int *) malloc(nlist->nelem*sizeof(int));
	}
	for (i=0; i < natom; i++){
		for(j=0; j < nlist->nelem; j++){
			soap_str->unique_Nneighbors_elemWise[i][j] = nlist->unique_Nneighbors_elemWise[i][j];
		}
	}

	soap_str->unique_neighborList_elemWise = (dyArray **) malloc(sizeof(dyArray*)*natom);
	for (i =0; i < natom; i++){
		soap_str->unique_neighborList_elemWise[i] = (dyArray *) malloc(sizeof(dyArray)*nlist->nelem);
	}

	for (i =0; i < natom; i++){
		for(j=0; j < nlist->nelem; j++){
			soap_str->unique_neighborList_elemWise[i][j].len = nlist->unique_neighborList_elemWise[i][j].len;
			soap_str->unique_neighborList_elemWise[i][j].capacity = nlist->unique_neighborList_elemWise[i][j].capacity;
			soap_str->unique_neighborList_elemWise[i][j].array = (int *)malloc(sizeof(int)*nlist->unique_neighborList_elemWise[i][j].len);
			for (k =0; k < nlist->unique_neighborList_elemWise[i][j].len; k++){
				soap_str->unique_neighborList_elemWise[i][j].array[k] = nlist->unique_neighborList_elemWise[i][j].array[k];
			}
		}
	}

	soap_str->neighborList = (dyArray *) malloc(sizeof(dyArray)*natom);
	for (i =0; i < natom; i++){
		soap_str->neighborList[i].len = nlist->neighborList[i].len;
		soap_str->neighborList[i].capacity = nlist->neighborList[i].capacity;
		soap_str->neighborList[i].array = (int *)malloc(sizeof(int)*nlist->neighborList[i].len);
		for (k =0; k<nlist->neighborList[i].len; k++){
			soap_str->neighborList[i].array[k] = nlist->neighborList[i].array[k];
		}
	}

	soap_str->unique_neighborList = (dyArray *) malloc(sizeof(dyArray)*natom);
	for (i =0; i < natom; i++){
		soap_str->unique_neighborList[i].len = nlist->unique_neighborList[i].len;
		soap_str->unique_neighborList[i].capacity = nlist->unique_neighborList[i].capacity;
		soap_str->unique_neighborList[i].array = (int *)malloc(sizeof(int)*nlist->unique_neighborList[i].len);
		for (k =0; k<nlist->unique_neighborList[i].len; k++){
			soap_str->unique_neighborList[i].array[k] = nlist->unique_neighborList[i].array[k];
		}
	}


	soap_str->natom_elem = (int *) malloc(sizeof(int)*nlist->nelem);
	for (i=0; i<nlist->nelem; i++){
		soap_str->natom_elem[i] = nlist->natom_elem[i];
	}

	soap_str->X2 = (double **) malloc(natom * sizeof(double*));
	soap_str->X3 = (double **) malloc(natom * sizeof(double*));
	soap_str->dX2_dF = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX2_dX = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX2_dY = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX2_dZ = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX3_dF = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX3_dX = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX3_dY = (double ***) malloc(natom * sizeof(double**));
	soap_str->dX3_dZ = (double ***) malloc(natom * sizeof(double**));

	for (i=0; i < natom; i++){
		soap_str->X2[i] = (double *) malloc(size_X2 * sizeof(double));
		soap_str->X3[i] = (double *) malloc(size_X3 * sizeof(double));
		for (int sz=0; sz < size_X2; sz++){
			soap_str->X2[i][sz] = 0.0;
		}
		for (int sz=0; sz < size_X3; sz++){
			soap_str->X3[i][sz] = 0.0;
		}

		int uniq_natms = uniqueEle((nlist->neighborList[i]).array, nlist->Nneighbors[i]);
		soap_str->dX2_dX[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));
		soap_str->dX2_dY[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));
		soap_str->dX2_dZ[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));

		soap_str->dX3_dX[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));
		soap_str->dX3_dY[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));
		soap_str->dX3_dZ[i] = (double **) malloc((1+uniq_natms) * sizeof(double*));



		soap_str->dX2_dF[i] = (double **) malloc(6 * sizeof(double*));
		soap_str->dX3_dF[i] = (double **) malloc(6 * sizeof(double*));
		for (j=0; j < 1+uniq_natms; j++){
			soap_str->dX2_dX[i][j] = (double *) malloc(size_X2 * sizeof(double));
			soap_str->dX2_dY[i][j] = (double *) malloc(size_X2 * sizeof(double));
			soap_str->dX2_dZ[i][j] = (double *) malloc(size_X2 * sizeof(double));
			soap_str->dX3_dX[i][j] = (double *) malloc(size_X3 * sizeof(double));
			soap_str->dX3_dY[i][j] = (double *) malloc(size_X3 * sizeof(double));
			soap_str->dX3_dZ[i][j] = (double *) malloc(size_X3 * sizeof(double));
			for (int sz=0; sz < size_X2; sz++){
				soap_str->dX2_dX[i][j][sz] = 0.0;
				soap_str->dX2_dY[i][j][sz] = 0.0;
				soap_str->dX2_dZ[i][j][sz] = 0.0;
			}
			for (int sz=0; sz < size_X3; sz++){
				soap_str->dX3_dX[i][j][sz] = 0.0;
				soap_str->dX3_dY[i][j][sz] = 0.0;
				soap_str->dX3_dZ[i][j][sz] = 0.0;
			}
		}
		for (j=0; j < 6; j++){
			soap_str->dX2_dF[i][j] = (double *) malloc(size_X2 * sizeof(double));
			soap_str->dX3_dF[i][j] = (double *) malloc(size_X3 * sizeof(double));
			for (int sz=0; sz < size_X2; sz++){
				soap_str->dX2_dF[i][j][sz] = 0.0;
			}
			for (int sz=0; sz < size_X3; sz++){
				soap_str->dX3_dF[i][j][sz] = 0.0;
			}
		}
	}

}



/*
delete_soapObj_wZ function frees the memory dynamically allocated for SoapObj 

[Input]
1. soap_str: pointer to SoapObj structure

[Output]
1. soap_str: pointer to SoapObj structure
*/

void delete_soapObj_wZ(SoapObj *soap_str){
	int size_cnlm, size_X2, size_X3, Nmod, nelem, i, j;
	
	// nelem = soap_str->nelem;
	nelem=1;
	Nmod = nelem * soap_str->Nmax;
	// size_cnlm = (soap_str->Lmax+1)*(soap_str->Lmax+1)*soap_str->Nmax;
	size_X3 = ( ((Nmod)*(1 + Nmod))/2)*(soap_str->Lmax+1);
	size_X2 = soap_str->Nmax * nelem;

	for (i=0; i<soap_str->natom; i++){
		free(soap_str->X2[i]);
		free(soap_str->X3[i]);
		free(soap_str->unique_Nneighbors_elemWise[i]);

		for (j=0; j < soap_str->nelem; j++){
			delete_dyarray(&(soap_str->unique_neighborList_elemWise[i][j]));
		}

		free(soap_str->unique_neighborList_elemWise[i]);
		int uniq_natms = uniqueEle((soap_str->neighborList[i]).array, soap_str->Nneighbors[i]);
		delete_dyarray(&(soap_str->neighborList[i]));
		delete_dyarray(&(soap_str->unique_neighborList[i]));
		for (j=0; j<1+uniq_natms; j++){
			free(soap_str->dX2_dX[i][j]);
			free(soap_str->dX2_dY[i][j]);
			free(soap_str->dX2_dZ[i][j]);
			free(soap_str->dX3_dX[i][j]);
			free(soap_str->dX3_dY[i][j]);
			free(soap_str->dX3_dZ[i][j]);
		}
		for (j=0; j<6; j++){
			free(soap_str->dX2_dF[i][j]);
			free(soap_str->dX3_dF[i][j]);
		}
		free(soap_str->dX2_dX[i]);
		free(soap_str->dX2_dY[i]);
		free(soap_str->dX2_dZ[i]);
		free(soap_str->dX3_dX[i]);
		free(soap_str->dX3_dY[i]);
		free(soap_str->dX3_dZ[i]);
		free(soap_str->dX2_dF[i]);
		free(soap_str->dX3_dF[i]);
	}

	free(soap_str->unique_Nneighbors);
	free(soap_str->Nneighbors);
	free(soap_str->natom_elem);
	free(soap_str->neighborList);
	free(soap_str->unique_neighborList);
	free(soap_str->unique_neighborList_elemWise);

	free(soap_str->X2);
	free(soap_str->X3);
	free(soap_str->dX2_dX);
	free(soap_str->dX2_dY);
	free(soap_str->dX2_dZ);
	free(soap_str->dX3_dX);
	free(soap_str->dX3_dY);
	free(soap_str->dX3_dZ);
	free(soap_str->dX2_dF);
	free(soap_str->dX3_dF);
}




/*
build_soapObj_wZ function calculates the SOAP descriptors and their derivatives w.r.t the atom positions and the deformation gradient

[Input]
1. nlist: pointer to the NeighList structure
2. rgrid: pointer to the rgrid for the spline interpolation
3. h_nl: pointer to the h_nl function for spline interpolation
4. dh_nl: pointer to the derivative of h_nl for spline interpolation
5. atompos: pointer to the atom positions [stored in ColMajor]

[Output]
1. soap_str: pointer to the SoapObj structure
*/

void build_soapObj_wZ(SoapObj *soap_str, NeighList *nlist, double* rgrid, double* h_nl, double* dh_nl, double *atompos, int Nmax, int Lmax, double beta_3, double xi_3, int *Z) {


	const double PI = 3.141592653589793;
	int i, j, k, N_r, elem_typ1;
	double xi, yi, zi, xj, yj, zj;
	double dx,dy,dz,dr,dtheta,dphi;
	double L1, L2, L3;
	double complex *Y, *dY_theta, *dY_phi;
	double *YD_hnl, *YD_dhnl;
	double cth, sth, cph, sph, dx2dy2z, dx2dy2;
	double dtheta_dx, dtheta_dy, dtheta_dz, dphi_dx, dphi_dy, dphi_dz;
	double dr_dF11, dr_dF22, dr_dF33, dr_dF12, dr_dF23, dr_dF13;
	double dth_dF11, dth_dF22, dth_dF33, dth_dF12, dth_dF23, dth_dF13;
	double dphi_dF11, dphi_dF22, dphi_dF33, dphi_dF12, dphi_dF23, dphi_dF13, dxdr, dydr, dzdr;
	int idx1, idx2, idx3, idx_neigh, elem_typ, n, l, m;
	double *hnl_temp, *dhnl_temp, *rtemp;
	int *ntemp, *imgx_temp, *imgy_temp, *imgz_temp, *elemtyp_temp;
	double complex ***cnlm, ****dcnlm_dX, ****dcnlm_dY, ****dcnlm_dZ, ****dcnlm_dF;
	int natom= nlist->natom, nelem, uniq_el_atms, local_index;
	nelem = 1;
	int size_cnlm =  (Lmax + 1) * (Lmax + 2) * Nmax;
	size_cnlm = size_cnlm/2;




	initialize_soapObj_wZ(soap_str, nlist, Lmax, Nmax, beta_3, xi_3);
	cnlm = (double complex***) malloc(natom*sizeof(double complex**));
	dcnlm_dX = (double complex****) malloc(natom*sizeof(double complex***));
	dcnlm_dY = (double complex****) malloc(natom*sizeof(double complex***));
	dcnlm_dZ = (double complex****) malloc(natom*sizeof(double complex***));
	dcnlm_dF = (double complex****) malloc(natom*sizeof(double complex***));

	for (i = 0; i < natom; i++){
		cnlm[i] = (double complex**) malloc(nelem*sizeof(double complex*));
		for (j = 0; j < nelem; j++){
			cnlm[i][j] = (double complex*) malloc(size_cnlm*sizeof(double complex));
			for (k = 0; k < size_cnlm; k++){
				cnlm[i][j][k] = 0.0 + 0.0*I;
			}
		}
		dcnlm_dX[i] = (double complex***) malloc(nelem*sizeof(double complex**));
		dcnlm_dY[i] = (double complex***) malloc(nelem*sizeof(double complex**));
		dcnlm_dZ[i] = (double complex***) malloc(nelem*sizeof(double complex**));
		dcnlm_dF[i] = (double complex***) malloc(nelem*sizeof(double complex**));

		for (j = 0; j < nelem; j++){
			uniq_el_atms=0;
			for (int jj = 0; jj < soap_str->nelem; jj++){
				uniq_el_atms += nlist->unique_Nneighbors_elemWise[i][jj];
			}
			// uniq_el_atms = nlist->unique_Nneighbors_elemWise[i][j];
			dcnlm_dX[i][j] = (double complex**) malloc((1+uniq_el_atms)*sizeof(double complex*));
			dcnlm_dY[i][j] = (double complex**) malloc((1+uniq_el_atms)*sizeof(double complex*));
			dcnlm_dZ[i][j] = (double complex**) malloc((1+uniq_el_atms)*sizeof(double complex*));
			dcnlm_dF[i][j] = (double complex**) malloc(6*sizeof(double complex*));
			for (k = 0; k < 1+uniq_el_atms; k++){
				dcnlm_dX[i][j][k] = (double complex*) malloc(size_cnlm*sizeof(double complex));
				dcnlm_dY[i][j][k] = (double complex*) malloc(size_cnlm*sizeof(double complex));
				dcnlm_dZ[i][j][k] = (double complex*) malloc(size_cnlm*sizeof(double complex));
				for (l = 0; l < size_cnlm; l++){
					dcnlm_dX[i][j][k][l] = 0.0 + 0.0*I;
					dcnlm_dY[i][j][k][l] = 0.0 + 0.0*I;
					dcnlm_dZ[i][j][k][l] = 0.0 + 0.0*I;
				}
			}
			for (k = 0; k < 6; k++){
				dcnlm_dF[i][j][k] = (double complex*) malloc(size_cnlm*sizeof(double complex));
				for (l = 0; l < size_cnlm; l++){
					dcnlm_dF[i][j][k][l] = 0.0 + 0.0*I;
				}
			}
		}
	}

	N_r = soap_str->N_rgrid;
	L1 = soap_str->cell[0];
	L2 = soap_str->cell[1];
	L3 = soap_str->cell[2];
	Lmax = soap_str->Lmax;
	Nmax = soap_str->Nmax;



	Y = (double complex *) malloc((Lmax+1)*(Lmax+1)*sizeof(double complex));
	dY_theta = (double complex *) malloc((Lmax+1)*(Lmax+1)*sizeof(double complex));
	dY_phi = (double complex *) malloc((Lmax+1)*(Lmax+1)*sizeof(double complex));
	hnl_temp = (double *) malloc(sizeof(double)*Nmax*(Lmax+1));
	dhnl_temp = (double *) malloc(sizeof(double)*Nmax*(Lmax+1));
	rtemp = (double *) malloc(sizeof(double)*1);

	YD_hnl = (double *) malloc(N_r*(Nmax)*(Lmax+1)*sizeof(double));
	YD_dhnl = (double *) malloc(N_r*(Nmax)*(Lmax+1)*sizeof(double));

	// derivatives of tabulated h_nl, dh_nl required for spline interpolation (used later in the inner loop).
	for (i =0; i < Nmax; i++){
		for (j =0; j < Lmax+1; j++){
			getYD_gen(rgrid, h_nl + (i*(Lmax+1)+j)*N_r, YD_hnl + (i*(Lmax+1)+j)*N_r, N_r);
			getYD_gen(rgrid, dh_nl + (i*(Lmax+1)+j)*N_r, YD_dhnl + (i*(Lmax+1)+j)*N_r, N_r);
		}
	}

	for (i = 0; i < soap_str->natom; i++){
		xi = atompos[3*i];
		yi = atompos[3*i+1];
		zi = atompos[3*i+2];
		// ntemp = (nlist->neighborList +i)->array;
		// imgx_temp=(nlist->neighborList_imgX + i)->array;
		// imgy_temp=(nlist->neighborList_imgY + i)->array;
		// imgz_temp=(nlist->neighborList_imgZ + i)->array;
		// elemtyp_temp=(nlist->neighborAtmTyp +i)->array;

		for (j = 0; j < nlist->Nneighbors[i]; j++){
			// idx_neigh = ntemp[j];
			// elem_typ = elemtyp_temp[j];
			// xj = atompos[3*idx_neigh] + L1 * imgx_temp[j];
			// yj = atompos[3*idx_neigh+1] + L2 * imgy_temp[j];
			// zj = atompos[3*idx_neigh+2] + L3 * imgz_temp[j];
			idx_neigh = (nlist->neighborList +i)->array[j];
			elem_typ = (nlist->neighborAtmTyp +i)->array[j];
			xj = atompos[3*idx_neigh] + L1 * (nlist->neighborList_imgX + i)->array[j];
			yj = atompos[3*idx_neigh+1] + L2 * (nlist->neighborList_imgY + i)->array[j];
			zj = atompos[3*idx_neigh+2] + L3 * (nlist->neighborList_imgZ + i)->array[j];

			dx = xj - xi;
			dy = yj - yi;
			dz = zj - zi;

			dr = sqrt(dx*dx + dy*dy + dz*dz);

			dtheta = atan(sqrt(dx*dx + dy*dy)/dz);
			if (dz==0) dtheta = PI/2;
			if (dtheta<0) dtheta = PI + dtheta;

			dphi = atan(dy/dx);
			if (dx>0 && dy>0) dphi = dphi;  // 1st Quandrant
			if (dx<0 && dy>0) dphi = PI+dphi;  // 2nd Quandrant
			if (dx<0 && dy<0) dphi = PI+dphi;  // 3rd Quadrant
			if (dx>0 && dy<0) dphi = 2*PI+dphi;  // 4th Quadrant
			if (dx==0 && dy >0) dphi = PI/2;
			if (dx==0 && dy <0) dphi = 3*PI/2;
			if (dx==0 && dy ==0) dphi = 0;


			cth = cos(dtheta);
			sth = sin (dtheta);

			dx2dy2z = dz*sqrt(dx*dx+dy*dy);
			dx2dy2 = sqrt(dx*dx+dy*dy);

			
			if (dx2dy2z == 0){
				dtheta_dx = 0.0;
				dtheta_dy = 0.0;
			} else {
				dtheta_dx = (dx*cth*cth)/(dx2dy2z);
				dtheta_dy = (dy*cth*cth)/(dx2dy2z);
			}

			if (dz == 0){
				dtheta_dz = 0.0;
			} else {
				dtheta_dz = -1*sth*cth / dz;
			}


			cph = cos(dphi);
			sph = sin (dphi);

			dphi_dz = 0;
			if (dx == 0){
				dphi_dx = 0.0;
				dphi_dy = 0.0;
			} else {
				dphi_dx = (-1*sph*cph)/dx;
				dphi_dy = (cph*cph)/dx;
			}

			dr_dF11 = dx*dx/dr;
			dr_dF22 = dy*dy/dr;
			dr_dF33 = dz*dz/dr;
			dr_dF12 = 2*dx*dy/dr;
			dr_dF23 = 2*dy*dz/dr;
			dr_dF13 = 2*dz*dx/dr;

			
			if (dx2dy2z==0){
				dth_dF11=0.0;
				dth_dF22=0.0;
				dth_dF12=0.0;
				dth_dF23=0.0;
				dth_dF13=0.0;

			} else {
				dth_dF11 = (dx*dx*cth*cth)/(dx2dy2z);
				dth_dF22 = (dy*dy*cth*cth)/(dx2dy2z);
				dth_dF12 = (2*dx*dy*cth*cth)/(dx2dy2z);
				dth_dF23 = (dy*cth*cth)/dx2dy2 - dy*cth*sth/dz;
				dth_dF13 = (dx*cth*cth)/dx2dy2 - dx*cth*sth/dz;
			}
			dth_dF33 = -1*sth*cth;
			

			dphi_dF11 = -1*sph*cph;
			if (dx==0){
				dphi_dF22=0.0;
				dphi_dF12=0.0;
				dphi_dF23=0.0;
				dphi_dF13=0.0;

			} else {
				dphi_dF22 = (dy*cph*cph)/dx;
				dphi_dF12 = (dx*cph*cph - sph*cph*dy)/dx;
				dphi_dF23 = cph*cph*dz/dx;
				dphi_dF13 = -1*sph*cph*dz/dx;
			}
			
			dphi_dF33 = 0;

			dxdr = dx/dr, dydr = dy/dr, dzdr = dz/dr;

			sph_harmonics(dtheta, dphi, soap_str->Lmax, Y, dY_theta, dY_phi);

			rtemp[0] = dr;
			// Calculate h_nl and dh_nl from dpline interpoation for all combinations of n and l
			for (n=0; n < Nmax; n++){
				for (l=0; l < Lmax+1; l++){
					SplineInterp(rgrid, h_nl+(n*(Lmax+1)+l)*N_r,
						N_r, rtemp, hnl_temp + n*(Lmax+1)+l, 1, YD_hnl+(n*(Lmax+1)+l)*N_r);

					SplineInterp(rgrid, dh_nl+(n*(Lmax+1)+l)*N_r,
						N_r, rtemp, dhnl_temp + n*(Lmax+1)+l, 1, YD_dhnl+(n*(Lmax+1)+l)*N_r);
				}
			}

			// int local_index = lin_search((nlist->unique_neighborList_elemWise[i][elem_typ]).array,
			 							// nlist->unique_Nneighbors_elemWise[i][elem_typ], idx_neigh);
			local_index = 1+nlist->localID_neighbours[i].array[j];

			for (n=0; n < Nmax; n++){
				for (l=0; l < Lmax+1; l++){
					for (m=0; m < l+1; m++){
						idx1 = (n*(Lmax+1)*(Lmax+2))/2 + (l*(l+1))/2 + m;
						idx2 = n*(Lmax+1)+l;
						idx3 = l*l + l + m;
						elem_typ1=0;
						cnlm[i][elem_typ1][idx1] +=  Z[elem_typ]*hnl_temp[idx2] * conj(Y[idx3]);

						if (local_index != -1){

							dcnlm_dX[i][elem_typ1][local_index][idx1] += Z[elem_typ]*dhnl_temp[idx2] * (dxdr)*conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dx + dY_phi[idx3]*dphi_dx);

							dcnlm_dY[i][elem_typ1][local_index][idx1] += Z[elem_typ]*dhnl_temp[idx2] * (dydr)* conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dy + dY_phi[idx3]*dphi_dy);

							dcnlm_dZ[i][elem_typ1][local_index][idx1] += Z[elem_typ]*dhnl_temp[idx2] * (dzdr)* conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dz + dY_phi[idx3]*dphi_dz);


							dcnlm_dX[i][elem_typ1][0][idx1] -= Z[elem_typ]*dhnl_temp[idx2] * (dxdr)*conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dx + dY_phi[idx3]*dphi_dx);

							dcnlm_dY[i][elem_typ1][0][idx1] -= Z[elem_typ]*dhnl_temp[idx2] * (dydr)* conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dy + dY_phi[idx3]*dphi_dy);

							dcnlm_dZ[i][elem_typ1][0][idx1] -= Z[elem_typ]*dhnl_temp[idx2] * (dzdr)* conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dtheta_dz + dY_phi[idx3]*dphi_dz);

							
						}
						
						
						dcnlm_dF[i][elem_typ1][0][idx1] += Z[elem_typ]*dhnl_temp[idx2]*dr_dF11 * conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF11 + dY_phi[idx3]*dphi_dF11);

						dcnlm_dF[i][elem_typ1][1][idx1] += Z[elem_typ]*dhnl_temp[idx2]*dr_dF22 * conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF22 + dY_phi[idx3]*dphi_dF22);

						dcnlm_dF[i][elem_typ1][2][idx1] += Z[elem_typ]*dhnl_temp[idx2]*dr_dF33 * conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF33 + dY_phi[idx3]*dphi_dF33);

						dcnlm_dF[i][elem_typ1][3][idx1] += dhnl_temp[idx2]*dr_dF12 * conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF12 + dY_phi[idx3]*dphi_dF12);

						dcnlm_dF[i][elem_typ1][4][idx1] += dhnl_temp[idx2]*dr_dF23 * conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF23 + dY_phi[idx3]*dphi_dF23);

						dcnlm_dF[i][elem_typ1][5][idx1] += dhnl_temp[idx2]*dr_dF13 * conj(Y[idx3]) +
									 Z[elem_typ]*hnl_temp[idx2] * conj(dY_theta[idx3]*dth_dF13 + dY_phi[idx3]*dphi_dF13);

					}
				}
			}
			
		}
	}

	double complex el1_xder, el1_yder, el1_zder, el2_xder, el2_yder, el2_zder,el_xder,el_yder,el_zder;
	int na, el, count_X2, neigh, count_X3;
	
	for (na = 0; na < natom; na++){
		 count_X2=0;
		for (el = 0; el < nelem; el++){
			for (n = 0; n < soap_str->Nmax; n++){
				idx1 = (n*(Lmax+1)*(Lmax+2))/2 ;

				soap_str->X2[na][count_X2] += (double) (cnlm[na][el][idx1]);
														 
				soap_str->dX2_dF[na][0][count_X2] += (double) (dcnlm_dF[na][el][0][idx1]);
														 
				soap_str->dX2_dF[na][1][count_X2] += (double) (dcnlm_dF[na][el][1][idx1]);
														 
				soap_str->dX2_dF[na][2][count_X2] += (double) (dcnlm_dF[na][el][2][idx1]);
														 
				soap_str->dX2_dF[na][3][count_X2] += (double) (dcnlm_dF[na][el][3][idx1]);
														 
				soap_str->dX2_dF[na][4][count_X2] += (double) (dcnlm_dF[na][el][4][idx1]);
														 
				soap_str->dX2_dF[na][5][count_X2] += (double) (dcnlm_dF[na][el][5][idx1]);
														 

				for (neigh = 0; neigh < nlist->unique_Nneighbors[na]; neigh++){

					// local_index  = lin_search((nlist->unique_neighborList_elemWise[na][el]).array,
 					// nlist->unique_Nneighbors_elemWise[na][el], (nlist->unique_neighborList[na]).array[neigh]);
					local_index = 1+nlist->localID_neighbours_elem[na][el].array[neigh];
			
					if (local_index == 0){
						el_xder = 0.0 + 0.0*I;
						el_yder = 0.0 + 0.0*I;
						el_zder = 0.0 + 0.0*I;
					} else{
						el_xder = dcnlm_dX[na][el][local_index][idx1];
						el_yder = dcnlm_dY[na][el][local_index][idx1];
						el_zder = dcnlm_dZ[na][el][local_index][idx1];
					}

					soap_str->dX2_dX[na][neigh+1][count_X2]  += (double)(el_xder);

					soap_str->dX2_dY[na][neigh+1][count_X2]  += (double)(el_yder);

					soap_str->dX2_dZ[na][neigh+1][count_X2]  +=  (double)(el_zder);
				}

				soap_str->dX2_dX[na][0][count_X2] +=  (double) dcnlm_dX[na][el][0][idx1];
				soap_str->dX2_dY[na][0][count_X2] +=  (double) dcnlm_dY[na][el][0][idx1];
				soap_str->dX2_dZ[na][0][count_X2] +=  (double) dcnlm_dZ[na][el][0][idx1];

				count_X2++;
			}
		}
	}


	
	int el1, n1, n2_el2, n2, el2, local_index2, local_index1;
	double const_X3, temp1, temp2, temp3;

	for (na = 0; na < natom; na++){
		count_X3=0;
		for (el1 = 0; el1 < nelem; el1++){
			for (n1 = 0; n1 < soap_str->Nmax; n1++){
				temp1 = (n1*(Lmax+1)*(Lmax+2))/2;
				for (n2_el2 = el1*soap_str->Nmax + n1; n2_el2 < nelem*soap_str->Nmax; n2_el2++){
					n2 = n2_el2 % soap_str->Nmax;
					el2 = n2_el2 / soap_str->Nmax;
					el2=0;
					temp2 = (n2*(Lmax+1)*(Lmax+2))/2;
					for (l = 0; l < soap_str->Lmax+1; l++){
						const_X3 = sqrt(8*PI*PI/(2*l+1));
						temp3 = (l*(l+1))/2;
						for (m = 0; m < l+1; m++){
							idx1 = temp1 + temp3 + m;
							idx2 = temp2 + temp3 + m;
							if (m==0){
								soap_str->X3[na][count_X3] += const_X3 * (double) (cnlm[na][el1][idx1]*conj(cnlm[na][el2][idx2]));
							} else {
								soap_str->X3[na][count_X3] += const_X3 * (double) (2* creal(cnlm[na][el1][idx1]*conj(cnlm[na][el2][idx2])));
							}
							// The derivative w.r.t to itself
							el1_xder = dcnlm_dX[na][el1][0][idx1];
							el1_yder = dcnlm_dY[na][el1][0][idx1];
							el1_zder = dcnlm_dZ[na][el1][0][idx1];
							el2_xder = dcnlm_dX[na][el2][0][idx2];
							el2_yder = dcnlm_dY[na][el2][0][idx2];
							el2_zder = dcnlm_dZ[na][el2][0][idx2];
							if (m==0){
								soap_str->dX3_dX[na][0][count_X3]  += const_X3 * (double)( (el1_xder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_xder)));

								soap_str->dX3_dY[na][0][count_X3]  += const_X3 * (double)( (el1_yder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_yder)));

								soap_str->dX3_dZ[na][0][count_X3]  += const_X3 * (double)( (el1_zder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_zder)));
							} else {
								soap_str->dX3_dX[na][0][count_X3]  += const_X3 * 2*(double)creal( (el1_xder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_xder)));

								soap_str->dX3_dY[na][0][count_X3]  += const_X3 * 2*(double)creal( (el1_yder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_yder)));

								soap_str->dX3_dZ[na][0][count_X3]  += const_X3 * 2*(double)creal( (el1_zder*conj(cnlm[na][el2][idx2])) +
														  (cnlm[na][el1][idx1]*conj(el2_zder)));
							}
							// The derivative w.r.t to itself end
							for (neigh = 0; neigh < nlist->unique_Nneighbors[na]; neigh++){
								local_index1 = 1+nlist->localID_neighbours_elem[na][el1].array[neigh];
								local_index2 = 1+nlist->localID_neighbours_elem[na][el2].array[neigh];
								if (local_index1 == 0){
									el1_xder = 0.0 + 0.0*I;
									el1_yder = 0.0 + 0.0*I;
									el1_zder = 0.0 + 0.0*I;
								} else{
									el1_xder = dcnlm_dX[na][el1][local_index1][idx1];
									el1_yder = dcnlm_dY[na][el1][local_index1][idx1];
									el1_zder = dcnlm_dZ[na][el1][local_index1][idx1];
								}
								if (local_index2 == 0){
									el2_xder = 0.0 + 0.0*I;
									el2_yder = 0.0 + 0.0*I;
									el2_zder = 0.0 + 0.0*I;
								} else{
									el2_xder = dcnlm_dX[na][el2][local_index2][idx2];
									el2_yder = dcnlm_dY[na][el2][local_index2][idx2];
									el2_zder = dcnlm_dZ[na][el2][local_index2][idx2];
								}

 
								if (m==0){
									soap_str->dX3_dX[na][1+neigh][count_X3]  += const_X3 * (double)( (el1_xder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_xder)));

									soap_str->dX3_dY[na][1+neigh][count_X3]  += const_X3 * (double)( (el1_yder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_yder)));

									soap_str->dX3_dZ[na][1+neigh][count_X3]  += const_X3 * (double)( (el1_zder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_zder)));
								} else {
									soap_str->dX3_dX[na][1+neigh][count_X3]  += const_X3 * 2*(double)creal( (el1_xder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_xder)));

									soap_str->dX3_dY[na][1+neigh][count_X3]  += const_X3 * 2*(double)creal( (el1_yder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_yder)));

									soap_str->dX3_dZ[na][1+neigh][count_X3]  += const_X3 * 2*(double)creal( (el1_zder*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(el2_zder)));
								}
							}
							if (m==0){
								soap_str->dX3_dF[na][0][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][0][idx1]*conj(cnlm[na][el2][idx2])) +
															 (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][0][idx2])));

								soap_str->dX3_dF[na][1][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][1][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][1][idx2])));

								soap_str->dX3_dF[na][2][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][2][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][2][idx2])));

								soap_str->dX3_dF[na][3][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][3][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][3][idx2])));

								soap_str->dX3_dF[na][4][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][4][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][4][idx2])));

								soap_str->dX3_dF[na][5][count_X3] += const_X3 * (double)( (dcnlm_dF[na][el1][5][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][5][idx2])));
							} else {
								soap_str->dX3_dF[na][0][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][0][idx1]*conj(cnlm[na][el2][idx2])) +
															 (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][0][idx2])));

								soap_str->dX3_dF[na][1][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][1][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][1][idx2])));

								soap_str->dX3_dF[na][2][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][2][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][2][idx2])));

								soap_str->dX3_dF[na][3][count_X3] += const_X3 * 2* (double)creal( (dcnlm_dF[na][el1][3][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][3][idx2])));

								soap_str->dX3_dF[na][4][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][4][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][4][idx2])));

								soap_str->dX3_dF[na][5][count_X3] += const_X3 * 2*(double)creal( (dcnlm_dF[na][el1][5][idx1]*conj(cnlm[na][el2][idx2])) +
															  (cnlm[na][el1][idx1]*conj(dcnlm_dF[na][el2][5][idx2])));
							}
						}
						count_X3++;
					}
				}
			}
			
		}

	}
	// if (1)
	// 	print_cnlm(cnlm[0], dcnlm_dX[0], dcnlm_dY[0], dcnlm_dZ[0], dcnlm_dF[0], size_cnlm, soap_str->unique_Nneighbors_elemWise[0], nelem);
	// free the memory for temperory Y, dY, hnl etc.
	free(Y);
	free(dY_theta);
	free(dY_phi);
	free(YD_hnl);
	free(YD_dhnl);
	free(hnl_temp);
	free(dhnl_temp);
	free(rtemp);

	// free the memory for cnlm
	for (i = 0; i < natom; i++){
		for (j = 0; j < nelem; j++){
			free(cnlm[i][j]);
		}
		free(cnlm[i]);

		for (j = 0; j < nelem; j++){
			uniq_el_atms = nlist->unique_Nneighbors_elemWise[i][j];
			for (k = 0; k < 1+uniq_el_atms; k++){
				free(dcnlm_dX[i][j][k]); free(dcnlm_dY[i][j][k]); free(dcnlm_dZ[i][j][k]);
			}
			for (k = 0; k < 6; k++){
				free(dcnlm_dF[i][j][k]);
			}
			free(dcnlm_dX[i][j]); free(dcnlm_dY[i][j]); free(dcnlm_dZ[i][j]); free(dcnlm_dF[i][j]); 
		}
		free(dcnlm_dX[i]); free(dcnlm_dY[i]); free(dcnlm_dZ[i]); free(dcnlm_dF[i]);
	}
	free(dcnlm_dX); free(dcnlm_dY); free(dcnlm_dZ); free(dcnlm_dF); free(cnlm);

}
