#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "helper.h"
#include "solid_harmonics.h"
#include "surface_harmonics.h"
#include "gmp_deriv_descriptor.h"
#include "local_type_tools.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

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

	initialize_nlist(nlist, natom, rcut, nelem);
    // Working - gives number of each type of element
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

initialize_GMPObj function initializes the objects in GMPObj structure and also allocates memory of dynamic arrays

[Input]
1. nlist: pointer to Neighlist structure
2. cal_atoms: pointer to list of indexes of atoms requiring descriptor calculations
3. cal_num: number of atoms requiring descriptor calculations
4. integer GMP params: mcsh order, square (1) or square root (0), solid (1) or surface (0) harmonic
5. double GMP params: gaussian probe sigma, 1 (unexplained, never referenced - will change later), sigma manipulation, sigma manipulation, sigma manipulation
6. nmcsh: total number of descriptors/atom
7. atom_gaussian: pointer to primitive gaussian parameters
8. ngaussians: pointer to number of primitives
9. element_index_to_order: param to convert atomic number to order in list of gaussian parameters
[Output]
1. gmp_str: pointer to GMPObj structure

*/

void initialize_gmpObj(GMPObj *gmp_str, NeighList *nlist, int *cal_atoms, int cal_num, int **params_i, double **params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order){
    int size_X, i,j,k;
    int nelem = nlist->nelem, natom = nlist->natom;
    size_X = nmcsh;
    gmp_str->size_X = size_X;
    gmp_str->cal_atoms = cal_atoms;
    gmp_str->cal_num = cal_num;
    gmp_str->params_i = params_i;
    gmp_str->params_d = params_d;
    gmp_str->nmcsh = nmcsh;
    gmp_str->atom_gaussian = atom_gaussian;
    gmp_str->ngaussians = ngaussians;
    gmp_str->element_index_to_order = element_index_to_order;

    gmp_str->natom = natom;
    gmp_str->nelem = nelem;
    gmp_str->rcut = nlist->rcut;
    gmp_str->cell[0] = nlist->cell[0];
    gmp_str->cell[1] = nlist->cell[1];
    gmp_str->cell[2] = nlist->cell[2];

    gmp_str->Nneighbors = (int *) malloc(sizeof(int)*natom);
    for (i=0; i < natom; i++){
        gmp_str->Nneighbors[i] = nlist->Nneighbors[i];

    }
    gmp_str->unique_Nneighbors = (int *) malloc(sizeof(int)*natom);
    for (i=0; i < natom; i ++){
        gmp_str->unique_Nneighbors[i] = nlist->unique_Nneighbors[i];

    }

    gmp_str->unique_Nneighbors_elemWise = (int **) malloc(natom*sizeof(int*));
	for (i=0; i < natom; i++){
		gmp_str->unique_Nneighbors_elemWise[i] = (int *) malloc(nelem*sizeof(int));
	}
	for (i=0; i < natom; i++){
		for(j=0; j < nelem; j++){
			gmp_str->unique_Nneighbors_elemWise[i][j] = nlist->unique_Nneighbors_elemWise[i][j];
		}
	}

	gmp_str->unique_neighborList_elemWise = (dyArray **) malloc(sizeof(dyArray*)*natom);
	for (i =0; i < natom; i++){
		gmp_str->unique_neighborList_elemWise[i] = (dyArray *) malloc(sizeof(dyArray)*nelem);
	}

	for (i =0; i < natom; i++){
		for(j=0; j < nelem; j++){
			gmp_str->unique_neighborList_elemWise[i][j].len = nlist->unique_neighborList_elemWise[i][j].len;
			gmp_str->unique_neighborList_elemWise[i][j].capacity = nlist->unique_neighborList_elemWise[i][j].capacity;
			gmp_str->unique_neighborList_elemWise[i][j].array = (int *)malloc(sizeof(int)*nlist->unique_neighborList_elemWise[i][j].len);
			for (k =0; k < nlist->unique_neighborList_elemWise[i][j].len; k++){
				gmp_str->unique_neighborList_elemWise[i][j].array[k] = nlist->unique_neighborList_elemWise[i][j].array[k];
			}
		}
	}

	gmp_str->neighborList = (dyArray *) malloc(sizeof(dyArray)*natom);
	for (i =0; i < natom; i++){
		gmp_str->neighborList[i].len = nlist->neighborList[i].len;
		gmp_str->neighborList[i].capacity = nlist->neighborList[i].capacity;
		gmp_str->neighborList[i].array = (int *)malloc(sizeof(int)*nlist->neighborList[i].len);
		for (k =0; k<nlist->neighborList[i].len; k++){
			gmp_str->neighborList[i].array[k] = nlist->neighborList[i].array[k];
		}
	}

	gmp_str->unique_neighborList = (dyArray *) malloc(sizeof(dyArray)*natom);
	for (i =0; i < natom; i++){
		gmp_str->unique_neighborList[i].len = nlist->unique_neighborList[i].len;
		gmp_str->unique_neighborList[i].capacity = nlist->unique_neighborList[i].capacity;
		gmp_str->unique_neighborList[i].array = (int *)malloc(sizeof(int)*nlist->unique_neighborList[i].len);
		for (k =0; k<nlist->unique_neighborList[i].len; k++){
			gmp_str->unique_neighborList[i].array[k] = nlist->unique_neighborList[i].array[k];
		}
	}


	gmp_str->natom_elem = (int *) malloc(sizeof(int)*nelem);
	for (i=0; i<nelem; i++){
		gmp_str->natom_elem[i] = nlist->natom_elem[i];
	}

    gmp_str->X = (double **) malloc(natom * sizeof(double*));
    gmp_str->dX_dX = (double ***) malloc(natom * sizeof(double**));
	gmp_str->dX_dY = (double ***) malloc(natom * sizeof(double**));
	gmp_str->dX_dZ = (double ***) malloc(natom * sizeof(double**));

    for (i=0; i < natom; i++){
		gmp_str->X[i] = (double *) malloc(size_X * sizeof(double));
		for (int sz=0; sz < size_X; sz++){
			gmp_str->X[i][sz] = 0.0;
		}

		int uniq_natms = uniqueEle((nlist->neighborList[i]).array, nlist->Nneighbors[i]);
		gmp_str->dX_dX[i] = (double **) malloc((uniq_natms) * sizeof(double*));
		gmp_str->dX_dY[i] = (double **) malloc((uniq_natms) * sizeof(double*));
		gmp_str->dX_dZ[i] = (double **) malloc((uniq_natms) * sizeof(double*));

		for (j=0; j < uniq_natms; j++){
			gmp_str->dX_dX[i][j] = (double *) malloc(size_X * sizeof(double));
			gmp_str->dX_dY[i][j] = (double *) malloc(size_X * sizeof(double));
			gmp_str->dX_dZ[i][j] = (double *) malloc(size_X * sizeof(double));
            for (int sz=0; sz < size_X; sz++){
				gmp_str->dX_dX[i][j][sz] = 0.0;
				gmp_str->dX_dY[i][j][sz] = 0.0;
				gmp_str->dX_dZ[i][j][sz] = 0.0;
			}
        }
    }
}

void build_gmpObj(GMPObj *gmp_str, NeighList *nlist, FeatureScaler *ftr_scale, int nmcsh, double *atompos, int **params_i, double **params_d, double** atom_gaussian, int* ngaussians, int* element_index_to_order, int* atom_type_to_indices, int* atom_indices){
    
    int cal_num = gmp_str->cal_num, train = ftr_scale->train;
    int i, j, k, N_r;
    int *cal_atoms = gmp_str->cal_atoms, *nneigh = gmp_str->Nneighbors;
    double xi, yi, zi, xj, yj, zj;
    int *ntemp, *imgx_temp, *imgy_temp, *imgz_temp, *elemtyp_temp;
    int idx1, idx2, idx3, idx_neigh, elem_typ, n, l, m;
    double L1, L2, L3;
    double x0,y0,z0,r0_sqr,dtheta,dphi;

    L1 = gmp_str->cell[0];
	L2 = gmp_str->cell[1];
	L3 = gmp_str->cell[2];

    for (int ii=0; ii < cal_num; ++ii) {
        i=cal_atoms[ii];       

		xi = atompos[3*i];
		yi = atompos[3*i+1];
		zi = atompos[3*i+2];

		ntemp = (nlist->neighborList +i)->array;
		imgx_temp=(nlist->neighborList_imgX + i)->array;
		imgy_temp=(nlist->neighborList_imgY + i)->array;
		imgz_temp=(nlist->neighborList_imgZ + i)->array;
		elemtyp_temp=(nlist->neighborAtmTyp +i)->array;
	 

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            double A = params_d[m][2], alpha = params_d[m][3], inv_rs = params_d[m][5];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                GMPFunction mcsh_function = get_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                if (mcsh_type == 1){
                    double sum_desc = 0.0;
					double *sum_dmiu_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
					double *sum_dmiu_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
					double *sum_dmiu_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));

					for (int j=0; j<(nneigh[i]+1); j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;
                    }

                    double m_desc[1], deriv[3];

                    for (int j = -1; j < nneigh[i]; ++j) {

                        int neigh_atom_element_order;
                        if (j < 0) {
                            xj = xi;
                            yj = yi;
                            zj = zi;
                            int neigh_atom_element_index = atom_indices[i];
                            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        }
                        else{
                            idx_neigh = ntemp[j];
                            elem_typ = elemtyp_temp[j];
                            xj = atompos[3*idx_neigh] + L1 * imgx_temp[j];
                            yj = atompos[3*idx_neigh+1] + L2 * imgy_temp[j];
                            zj = atompos[3*idx_neigh+2] + L3 * imgz_temp[j];
                            int neigh_atom_element_index = atom_type_to_indices[elem_typ];
                            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        }
						
                        x0 = xj - xi;
                        y0 = yj - yi;
                        z0 = zj - zi;

                        r0_sqr = x0*x0 + y0*y0 + z0*z0;

                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, m_desc, deriv);
                            sum_desc += m_desc[0];
							sum_dmiu_dxj[j+1] += deriv[0];
							sum_dmiu_dyj[j+1] += deriv[1];
							sum_dmiu_dzj[j+1] += deriv[2];
                        }
                    }
                    sum_square += group_coefficient * sum_desc * sum_desc;
					double dmdx, dmdy, dmdz;
                    for (int j = -1; j < nneigh[i]; ++j) {
                        dmdx = (sum_desc * sum_dmiu_dxj[j+1]) * group_coefficient * 2.0;
                        dmdy = (sum_desc * sum_dmiu_dyj[j+1]) * group_coefficient * 2.0;
                        dmdz = (sum_desc * sum_dmiu_dzj[j+1]) * group_coefficient * 2.0;
						if (j < 0) {
							gmp_str->dX_dX[ii][i][m] += dmdx;
							gmp_str->dX_dY[ii][i][m] += dmdy;
							gmp_str->dX_dZ[ii][i][m] += dmdz;
						}
						else {
							gmp_str->dX_dX[ii][(nlist->neighborList[i]).array[j]][m] += dmdx;
							gmp_str->dX_dY[ii][(nlist->neighborList[i]).array[j]][m] += dmdy;
							gmp_str->dX_dZ[ii][(nlist->neighborList[i]).array[j]][m] += dmdz;
						}
                        gmp_str->dX_dX[ii][i][m] -= dmdx;
                        gmp_str->dX_dY[ii][i][m] -= dmdy;
                        gmp_str->dX_dZ[ii][i][m] -= dmdz;						
                    }

                    free(sum_dmiu_dxj);
                    free(sum_dmiu_dyj);
                    free(sum_dmiu_dzj);
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu2_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu3_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu1_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu2_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu3_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu1_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu2_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu3_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));

                    for (int j=0; j<nneigh[i]+1; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                    }

                    double miu[3], deriv[9];

                    for (int j = -1; j < nneigh[i]; ++j) {
                        int neigh_atom_element_order;
                        if (j < 0) {
                            xj = xi;
                            yj = yi;
                            zj = zi;
                            int neigh_atom_element_index = atom_indices[i];
                            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        }
                        else{
                            idx_neigh = ntemp[j];
                            elem_typ = elemtyp_temp[j];
                            xj = atompos[3*idx_neigh] + L1 * imgx_temp[j];
                            yj = atompos[3*idx_neigh+1] + L2 * imgy_temp[j];
                            zj = atompos[3*idx_neigh+2] + L3 * imgz_temp[j];
                            int neigh_atom_element_index = atom_type_to_indices[elem_typ];
                            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        }

                        x0 = xj - xi;
                        y0 = yj - yi;
                        z0 = zj - zi;

                        r0_sqr = x0*x0 + y0*y0 + z0*z0;

                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu, deriv);
                            sum_miu1 += miu[0];
                            sum_miu2 += miu[1];
                            sum_miu3 += miu[2];
							sum_dmiu1_dxj[j+1] += deriv[0];
							sum_dmiu1_dyj[j+1] += deriv[1];
                            sum_dmiu1_dzj[j+1] += deriv[2];
                            sum_dmiu2_dxj[j+1] += deriv[3];
                            sum_dmiu2_dyj[j+1] += deriv[4];
                            sum_dmiu2_dzj[j+1] += deriv[5];
                            sum_dmiu3_dxj[j+1] += deriv[6];
                            sum_dmiu3_dyj[j+1] += deriv[7];
                            sum_dmiu3_dzj[j+1] += deriv[8];
                        }
					}
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
					
					double dmdx, dmdy, dmdz;
					for (int j = -1; j < nneigh[i]; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j+1] + sum_miu2 * sum_dmiu2_dxj[j+1] + sum_miu3 * sum_dmiu3_dxj[j+1]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j+1] + sum_miu2 * sum_dmiu2_dyj[j+1] + sum_miu3 * sum_dmiu3_dyj[j+1]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j+1] + sum_miu2 * sum_dmiu2_dzj[j+1] + sum_miu3 * sum_dmiu3_dzj[j+1]) * group_coefficient * 2.0;

                        if (j < 0) {
							gmp_str->dX_dX[ii][i][m] += dmdx;
							gmp_str->dX_dY[ii][i][m] += dmdy;
							gmp_str->dX_dZ[ii][i][m] += dmdz;
						}
						else {
							gmp_str->dX_dX[ii][(nlist->neighborList[i]).array[j]][m] += dmdx;
							gmp_str->dX_dY[ii][(nlist->neighborList[i]).array[j]][m] += dmdy;
							gmp_str->dX_dZ[ii][(nlist->neighborList[i]).array[j]][m] += dmdz;
						}

                        gmp_str->dX_dX[ii][i][m] -= dmdx;
                        gmp_str->dX_dY[ii][i][m] -= dmdy;
                        gmp_str->dX_dZ[ii][i][m] -= dmdz;
                    }

                    free(sum_dmiu1_dxj);
                    free(sum_dmiu1_dyj);
                    free(sum_dmiu1_dzj);
					free(sum_dmiu2_dxj);
                    free(sum_dmiu2_dyj);
                    free(sum_dmiu2_dzj);
					free(sum_dmiu3_dxj);
                    free(sum_dmiu3_dyj);
                    free(sum_dmiu3_dzj);
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

					double* sum_dmiu1_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu2_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu3_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu4_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu5_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu6_dxj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu1_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu2_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu3_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu4_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu5_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu6_dyj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu1_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu2_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu3_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu4_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu5_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));
                    double* sum_dmiu6_dzj = (double *) malloc(sizeof(double)*(nneigh[i]+1));

					for (int j=0; j<nneigh[i]+1; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;
					}
					
                    double miu[6], deriv[18];
                    int neigh_atom_element_order;
                    for (int j = -1; j < nneigh[i]; ++j) {
                        if (j < 0) {
                            xj = xi;
                            yj = yi;
                            zj = zi;
                            int neigh_atom_element_index = atom_indices[i];
                            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        }
                        else{
                            idx_neigh = ntemp[j];
                            elem_typ = elemtyp_temp[j];
                            xj = atompos[3*idx_neigh] + L1 * imgx_temp[j];
                            yj = atompos[3*idx_neigh+1] + L2 * imgy_temp[j];
                            zj = atompos[3*idx_neigh+2] + L3 * imgz_temp[j];
                            int neigh_atom_element_index = atom_type_to_indices[elem_typ];
                            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        }

                        x0 = xj - xi;
                        y0 = yj - yi;
                        z0 = zj - zi;

                        r0_sqr = x0*x0 + y0*y0 + z0*z0;

                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu, deriv);
                            sum_miu1 += miu[0];
                            sum_miu2 += miu[1];
                            sum_miu3 += miu[2];
                            sum_miu4 += miu[3];
                            sum_miu5 += miu[4];
                            sum_miu6 += miu[5];
							sum_dmiu1_dxj[j+1] += deriv[0];
                            sum_dmiu1_dyj[j+1] += deriv[1];
                            sum_dmiu1_dzj[j+1] += deriv[2];
                            sum_dmiu2_dxj[j+1] += deriv[3];
                            sum_dmiu2_dyj[j+1] += deriv[4];
                            sum_dmiu2_dzj[j+1] += deriv[5];
                            sum_dmiu3_dxj[j+1] += deriv[6];
                            sum_dmiu3_dyj[j+1] += deriv[7];
                            sum_dmiu3_dzj[j+1] += deriv[8];
                            sum_dmiu4_dxj[j+1] += deriv[9];
                            sum_dmiu4_dyj[j+1] += deriv[10];
                            sum_dmiu4_dzj[j+1] += deriv[11];
                            sum_dmiu5_dxj[j+1] += deriv[12];
                            sum_dmiu5_dyj[j+1] += deriv[13];
                            sum_dmiu5_dzj[j+1] += deriv[14];
                            sum_dmiu6_dxj[j+1] += deriv[15];
                            sum_dmiu6_dyj[j+1] += deriv[16];
                            sum_dmiu6_dzj[j+1] += deriv[17];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                    sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
					
					double dmdx, dmdy, dmdz;

                    for (int j = -1; j < nneigh[i]; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j+1] + sum_miu2 * sum_dmiu2_dxj[j+1] +
                                sum_miu3 * sum_dmiu3_dxj[j+1] + sum_miu4 * sum_dmiu4_dxj[j+1] +
                                sum_miu5 * sum_dmiu5_dxj[j+1] + sum_miu6 * sum_dmiu6_dxj[j+1]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j+1] + sum_miu2 * sum_dmiu2_dyj[j+1] +
                                sum_miu3 * sum_dmiu3_dyj[j+1] + sum_miu4 * sum_dmiu4_dyj[j+1] +
                                sum_miu5 * sum_dmiu5_dyj[j+1] + sum_miu6 * sum_dmiu6_dyj[j+1]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j+1] + sum_miu2 * sum_dmiu2_dzj[j+1] +
                                sum_miu3 * sum_dmiu3_dzj[j+1] + sum_miu4 * sum_dmiu4_dzj[j+1] +
                                sum_miu5 * sum_dmiu5_dzj[j+1] + sum_miu6 * sum_dmiu6_dzj[j+1]) * group_coefficient * 2.0;

						if (j < 0) {
							gmp_str->dX_dX[i][i][m] += dmdx;
							gmp_str->dX_dY[i][i][m] += dmdy;
							gmp_str->dX_dZ[i][i][m] += dmdz;
						}
						else {
							gmp_str->dX_dX[i][(nlist->neighborList[i]).array[j]][m] += dmdx;
							gmp_str->dX_dY[i][(nlist->neighborList[i]).array[j]][m] += dmdy;
							gmp_str->dX_dZ[i][(nlist->neighborList[i]).array[j]][m] += dmdz;
						}

                        gmp_str->dX_dX[i][i][m] -= dmdx;
                        gmp_str->dX_dY[i][i][m] -= dmdy;
                        gmp_str->dX_dZ[i][i][m] -= dmdz;
					}
                }
            }

            if (square != 0){
                gmp_str->X[ii][m] = sum_square;
            }
            else {
                gmp_str->X[ii][m] = sqrt(sum_square);
            }

        }
    }

	if (train == 0){
		scale_features(gmp_str, ftr_scale, nmcsh);
	}
	else {
		ftr_scale->mean_colWise = (double *) malloc(sizeof(double)*nmcsh);
		ftr_scale->stdev_colWise = (double *) malloc(sizeof(double)*nmcsh);
		double means[nmcsh], sums_sq[nmcsh];
		for (int i = 0; i < nmcsh; i++){
			means[i] = 0.;
			sums_sq[i] = 0.;
		}
		for (int i = 0; i < cal_num; i++){
			for (int j = 0; j < nmcsh; j++){
				means[j] += (gmp_str->X[i][j])/cal_num;
			}
		}
		for (int i = 0; i < cal_num; i++){
			for (int j = 0; j < nmcsh; j++){
				sums_sq[j] += pow((gmp_str->X[i][j])-means[j],2);
			}
		}
		for (int i = 0; i < nmcsh; i++){
			ftr_scale->mean_colWise[i] = means[i];
			double tmp_stdev = sqrt(sums_sq[i]/cal_num);
			if (tmp_stdev < .005) {
				ftr_scale->stdev_colWise[i] = 1.;
			} else {
				ftr_scale->stdev_colWise[i] = tmp_stdev;
			}
		}
		scale_features(gmp_str, ftr_scale, nmcsh);
	}

}

void scale_features(GMPObj *gmp_str, FeatureScaler *ftr_scale, int nmcsh){

	int cal_num = gmp_str->cal_num;
	int atom_num = gmp_str->natom;

	for (int i = 0; i < cal_num; i++){
		for (int j = 0; j < nmcsh; j++){
			double desc_entry = gmp_str->X[i][j], mean = ftr_scale->mean_colWise[j], stdev = ftr_scale->stdev_colWise[j];
			gmp_str->X[i][j] = (desc_entry-mean)/stdev;
			for (int k = 0; k < atom_num; k++){
				double dmdx = gmp_str->dX_dX[i][k][j], dmdy = gmp_str->dX_dY[i][k][j], dmdz = gmp_str->dX_dZ[i][k][j];
				gmp_str->dX_dX[i][k][j] = dmdx/stdev, gmp_str->dX_dY[i][k][j] = dmdy/stdev, gmp_str->dX_dZ[i][k][j] = dmdz/stdev;
			}
		}
	}
}