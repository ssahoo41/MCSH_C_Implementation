#include <stdio.h>
#include <stdlib.h>
#include "mlff_types.h"
#include "gmp_descriptor.h"

int main(){

    double atom_pos[2*3] = {6.0,6.0,6.0,8.0,8.0,8.0};
    double *atom_pos_p = atom_pos;

    int atom_indices[2];
    for (int i = 0; i < 2; i ++) atom_indices[i] = 13;
    int *atom_indices_p = atom_indices;
    int atom_type_to_indices[1] = {13};
    int *atom_type_to_indices_p = atom_type_to_indices;
    int atom_num = 2;
    int cal_atoms[2];
    for (int i = 0; i < 2; i++){
        cal_atoms[i] = i;
    }
    int *cal_atoms_p = cal_atoms;
    int cal_num = 2;

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
    1.80000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    1.80000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    1.80000000e+01,1.00000000e+00},
    {1.00000000e-01,1.00000000e+00,6.34936359e+01,5.00000000e+01,
    1.80000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    1.80000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    1.80000000e+01,1.00000000e+00},
    {1.00000000e-01,1.00000000e+00,6.34936359e+01,5.00000000e+01,
    1.80000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    1.80000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    1.80000000e+01,1.00000000e+00},
    {1.00000000e-01,1.00000000e+00,6.34936359e+01,5.00000000e+01,
    1.80000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    1.80000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    1.80000000e+01,1.00000000e+00},
    {1.00000000e-01,1.00000000e+00,6.34936359e+01,5.00000000e+01,
    1.80000000e+01,1.00000000e+00},
    {6.30000000e-01,1.00000000e+00,2.53926805e-01,1.25976316e+00,
    1.80000000e+01,1.00000000e+00},
    {3.98000000e+00,1.00000000e+00,1.00711945e-03,3.15648595e-02,
    1.80000000e+01,1.00000000e+00}};

    double *params_dp[15];
    double **params_dp_all;
    for (int i =0; i < 15; i ++){
            params_dp[i] = params_d[i];
        }
    params_dp_all = params_dp;

    int nmcsh = 15;

    double atom_gaussian[1][8] = {{0.3730494069109985,0.7344318799434225,516.9378247134416,13.577027904536275,
    889298.9729349068,12.0315655354449,-889796.765560644,12.032681865001077}};

    double *atom_gaussian_p[1];
    double **atom_gaussian_all;
    for (int i =0; i < 1; i ++){
            atom_gaussian_p[i] = atom_gaussian[i];
        }
    atom_gaussian_all = atom_gaussian_p;

    int ngaussians[1] = {4};
    int *ngaussians_p = ngaussians;
    int element_index_to_order[120] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int *element_index_to_order_p = element_index_to_order;

    //printf("Made it to calc");

    NeighList *nlist;
	GMPObj *gmp_str;
    FeatureScaler *ftr_scale;
	int *Z, *BC, *atomtyp;
    int Ntypes = 1, n_atom = 2, nAtomv[1] = {2};
	Z = (int *)malloc(Ntypes*sizeof(int));
	for (int i=0; i <Ntypes; i++){
		Z[i] = 1;//pSPARC->Znucl[i];
	}

	nlist = (NeighList *) malloc(sizeof(NeighList)*1);
	gmp_str = (GMPObj *) malloc(sizeof(GMPObj)*1);
    ftr_scale = (FeatureScaler *) malloc(sizeof(FeatureScaler)*1);

	double *cell, F_estimate_max;
    double rcut = 18;

	BC = (int*)malloc(3*sizeof(int));
	cell =(double *) malloc(3*sizeof(double));
	atomtyp = (int *) malloc(n_atom*sizeof(int));
	BC[0] = 1; BC[1] = 1; BC[2] = 1;
	cell[0] = 12.15;
	cell[1] = 12.15;
	cell[2] = 12.15;
    
	int count = 0;
	for (int i=0; i < Ntypes; i++){
		for (int j=0; j < nAtomv[i]; j++){
			atomtyp[count] = i;
			count++;
		}
	}

	//printf("sparc_mlff_interface_firstMD\n");
	//printf("--------------------\n");
	//printf("Initializing nlist ... \n");
	initialize_nlist(nlist, n_atom, rcut, Ntypes);
	//printf("Initialization of nlist done.\n");
	//printf("Building nlist ... \n");
	build_nlist(rcut, Ntypes, n_atom, atom_pos_p, atomtyp, BC, cell, nlist);
	//printf("Building nlist done. \n");
    ftr_scale->train = 1;
    //printf("Initializing GMPObj. \n");
    initialize_gmpObj(gmp_str, nlist, cal_atoms, cal_num,params_ip_all, params_dp_all, nmcsh, atom_gaussian_all, ngaussians, element_index_to_order_p);
    //printf("GMPObj initialized");
    
    //printf("Building soapxx ... \n");
    build_gmpObj(gmp_str, nlist, ftr_scale, nmcsh, atom_pos_p, params_ip_all, params_dp_all, atom_gaussian_all, ngaussians, element_index_to_order_p,atom_type_to_indices_p,atom_indices_p);
    //printf("Building soap donexx. \n");
    for (int i=0; i < 2; i++){
        for (int j=0; j < 15; j++){
            //printf("%d,%d,%f",i,j,gmp_str->X[i][j]);
            for (int k=0; k<2; k++){
                printf(",%f,%f,%f",gmp_str->dX_dX[i][k][j],gmp_str->dX_dY[i][k][j],gmp_str->dX_dZ[i][k][j]);
            }
            printf("\n");
        }
    }
    return 0;
}