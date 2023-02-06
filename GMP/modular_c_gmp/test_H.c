#include <stdio.h>
#include <stdlib.h>
#include "mlff_types.h"
#include "gmp_descriptor.h"

int main(){

    double atom_pos[1*3] = {51.0,51.0,51.0};

    double *atom_pos_p = atom_pos;

    int atom_indices[1] = {1};
    int *atom_indices_p = atom_indices;
    //This is poorly named. It takes SPARC atom type {0,1,2, etc.} to atom type
    int atom_type_to_indices[1] = {1};
    int *atom_type_to_indices_p = atom_type_to_indices;
    int atom_num = 1;
    int cal_atoms[1] = {0};
    int *cal_atoms_p = cal_atoms;
    int cal_num = 1;

    int params_i[15][3] = {{0,1,0},
    {0,1,0},
    {0,1,0},
    {1,1,0},
    {1,1,0},
    {1,1,0},
    {2,1,0},
    {2,1,0},
    {2,1,0},
    {3,1,0},
    {3,1,0},
    {3,1,0},
    {4,1,0},
    {4,1,0},
    {4,1,0}};
    int *params_ip[15];
    int **params_ip_all;
    for (int i =0; i < 15; i ++){
            params_ip[i] = params_i[i];
        }
    params_ip_all = params_ip;

    double params_d[15][6] = {{1.00000000e-01,1.00000000e+00,6.34936359e+01,5.00000000e+01,
    5.0000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    5.0000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    5.0000000e+01,1.00000000e+00},
    {1.00000000e-01,1.00000000e+00,6.34936359e+01,5.00000000e+01,
    5.0000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    5.0000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    5.0000000e+01,1.00000000e+00},
    {1.00000000e-01,1.00000000e+00,6.34936359e+01,5.00000000e+01,
    5.0000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    5.0000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    5.0000000e+01,1.00000000e+00},
    {1.00000000e-01,1.00000000e+00,6.34936359e+01,5.00000000e+01,
    5.0000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    5.0000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    5.0000000e+01,1.00000000e+00},
    {1.00000000e-01,1.00000000e+00,6.34936359e+01,5.00000000e+01,
    5.0000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    5.0000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    5.0000000e+01,1.00000000e+00}};

    double *params_dp[15];
    double **params_dp_all;
    for (int i =0; i < 15; i ++){
            params_dp[i] = params_d[i];
        }
    params_dp_all = params_dp;

    int nmcsh = 15;

    double atom_gaussian[7][12] = {{1.18508086e+00,1.09138014e+00,2.33780444e+03,1.19153115e+01
    ,6.32016327e+01,5.80717093e+00,-2.39524586e+03,1.18145155e+01
    ,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00},
    {-2.11565175e+03,4.76796442e+00,-6.03186368e+02,1.49540036e+01
    ,2.14585000e+03,4.76800516e+00,1.02093893e+03,1.59805744e+01
    ,-4.42445584e+02,1.69594315e+01,1.61060434e+00,1.15246327e+00},
    {-1.89497893e+00,1.48991693e+01,3.41152646e-01,1.02616758e+00
    ,-3.80495085e-01,3.24116879e+01,2.19001711e+00,2.98127984e+00
    ,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00},
    {3.58183509e-01,1.74903622e+00,9.54354963e-01,1.15139501e+01
    ,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00
    ,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00},
    {1.24169482e+00,1.08471703e+00,2.42445030e+03,1.13353153e+01
    ,6.74490058e+01,6.01501714e+00,-2.48681769e+03,1.12428153e+01
    ,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00},
    {-2.18105985e+03,4.48450788e+00,-5.71952362e+02,1.53792106e+01
    ,2.20526103e+03,4.48450331e+00,1.24072188e+03,1.62906631e+01
    ,-6.88388280e+02,1.68542970e+01,1.38643778e+00,1.07260777e+00},
    {8.36619396e-01,9.13560275e-01,2.62426625e+03,8.31445424e+00
    ,7.03676197e+01,5.35612632e+00,-2.68985158e+03,8.27276549e+00
    ,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00}};

    double *atom_gaussian_p[7];
    double **atom_gaussian_all;
    for (int i =0; i < 7; i ++){
            atom_gaussian_p[i] = atom_gaussian[i];
        }
    atom_gaussian_all = atom_gaussian_p;

    int ngaussians[7] = {4,6,4,2,4,6,4};
    int *ngaussians_p = ngaussians;
    int element_index_to_order[120] = {0,3,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int *element_index_to_order_p = element_index_to_order;

    NeighList *nlist;
	GMPObj *gmp_str;
	int *Z, *BC, *atomtyp;
    //Types of element, total atom number, number of eachtype of element
    int Ntypes = 1, n_atom = 1, nAtomv[1] = {1};
	Z = (int *)malloc(Ntypes*sizeof(int));
	for (int i=0; i <Ntypes; i++){
		Z[i] = 1;//pSPARC->Znucl[i];
	}

	nlist = (NeighList *) malloc(sizeof(NeighList)*1);
	gmp_str = (GMPObj *) malloc(sizeof(GMPObj)*1);

	double *cell, F_estimate_max;
    double rcut = 50;

	BC = (int*)malloc(3*sizeof(int));
	cell =(double *) malloc(3*sizeof(double));
	atomtyp = (int *) malloc(n_atom*sizeof(int));
	BC[0] = 1; BC[1] = 1; BC[2] = 1;
	cell[0] = 102;
	cell[1] = 102;
	cell[2] = 102;
    
	int count = 0;
	for (int i=0; i < Ntypes; i++){
		for (int j=0; j < nAtomv[i]; j++){
			atomtyp[count] = i;
			count++;
		}
	}

	printf("sparc_mlff_interface_firstMD\n");
	// printf("--------------------\n");
	printf("Initializing nlist ... \n");
	initialize_nlist(nlist, n_atom, rcut, Ntypes);
	printf("Initialization of nlist done.\n");
	printf("Building nlist ... \n");
	build_nlist(rcut, Ntypes, n_atom, atom_pos_p, atomtyp, BC, cell, nlist);
	printf("Building nlist done. \n");

    printf("Initializing GMPObj. \n");
    initialize_gmpObj(gmp_str, nlist, cal_atoms, cal_num,params_ip_all, params_dp_all, nmcsh, atom_gaussian_all, ngaussians, element_index_to_order_p);
    printf("GMPObj initialized");
    
    printf("Building soapxx ... \n");
    build_gmpObj(gmp_str, nlist, nmcsh, atom_pos_p, params_ip_all, params_dp_all, atom_gaussian_all, ngaussians, element_index_to_order_p,atom_type_to_indices_p,atom_indices_p);
    printf("Building soap donexx. \n");
    for (int i=0; i < 1; i++){
        for (int j=0; j < 15; j++){
            printf("%d,%d: %f\n",i,j,gmp_str->X[i][j]);
        }
    }
    return 0;
}