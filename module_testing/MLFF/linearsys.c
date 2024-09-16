#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "mkl.h"
#include "tools_mlff.h"
#include "spherical_harmonics.h"
#include "soap_descriptor.h"
#include "mlff_types.h"
#include "sparsification.h"
#include "regression.h"
#include "isddft.h"
#include "ddbp_tools.h"
#include "linearsys.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define au2GPa 29421.02648438959

/*
soap_kernel function computes the SOAP kernel between two descriptors {X2_str, X3_str} and {X2_tr, X3_tr}

[Input]
1. X2_str: pointer to the first descriptor (2-body)
2. X3_str: pointer to the first descriptor (3-body)
3. X2_tr: pointer to the second descriptor (2-body)
4. X3_tr: pointer to the second descriptor (3-body)
5. beta_2: weight to the 2-body term in the kernel
6. beta_3: weight to the 3-body term in the kernel
7. xi_3: exponent in the kernel
8. size_X2: length of the 2-body kernel
9. size_X3: length of the 3-body kernel
[Output]
1. kernel_val: Value of the kernel
*/


double soap_kernel(int kernel_typ, double *X2_str, double *X3_str, double *X2_tr, double *X3_tr,
			 double beta_2, double beta_3, double xi_3, int size_X2, int size_X3) {
	if (kernel_typ ==0)
		return soap_kernel_polynomial(X2_str, X3_str, X2_tr, X3_tr, beta_2, beta_3, xi_3, size_X2, size_X3);
	else if (kernel_typ==1)
		return soap_kernel_Gaussian(X2_str, X3_str, X2_tr, X3_tr, beta_2, beta_3, xi_3, size_X2, size_X3);
	else
		return soap_kernel_Laplacian(X2_str, X3_str, X2_tr, X3_tr, beta_2, beta_3, xi_3, size_X2, size_X3);
}

double der_soap_kernel(int kernel_typ, double *dX2_str, double *dX3_str, double *X2_str, double *X3_str, double *X2_tr, double *X3_tr,
			 double beta_2, double beta_3, double xi_3, int size_X2, int size_X3) {
	// printf("kernel_typ: %d, der_soap_kernel_polynomial: %f\n",kernel_typ,der_soap_kernel_polynomial(dX2_str, dX3_str, X2_str, X3_str, X2_tr, X3_tr, beta_2, beta_3, xi_3, size_X2, size_X3));
	// exit(1);
	if (kernel_typ ==0)
		return der_soap_kernel_polynomial(dX2_str, dX3_str, X2_str, X3_str, X2_tr, X3_tr, beta_2, beta_3, xi_3, size_X2, size_X3);
	else if (kernel_typ==1)
		return der_soap_kernel_Gaussian(dX2_str, dX3_str, X2_str, X3_str, X2_tr, X3_tr, beta_2, beta_3, xi_3, size_X2, size_X3);
	else
		return der_soap_kernel_Laplacian(dX2_str, dX3_str, X2_str, X3_str, X2_tr, X3_tr, beta_2, beta_3, xi_3, size_X2, size_X3);
}

double soap_kernel_polynomial(double *X2_str, double *X3_str, double *X2_tr, double *X3_tr,
			 double beta_2, double beta_3, double xi_3, int size_X2, int size_X3) {

	double norm_X3_str, norm_X3_tr, X3_str_temp[size_X3], X3_tr_temp[size_X3], kernel_val;
	int i;
	norm_X3_str = sqrt(dotProduct(X3_str, X3_str, size_X3));
	norm_X3_tr = sqrt(dotProduct(X3_tr, X3_tr, size_X3));
	for (i = 0; i<size_X3; i++){
		X3_str_temp[i] = X3_str[i]/norm_X3_str;
		X3_tr_temp[i] = X3_tr[i]/norm_X3_tr;
	}

	kernel_val = beta_2 * dotProduct(X2_str, X2_tr, size_X2) + beta_3 * pow(dotProduct(X3_tr_temp, X3_str_temp, size_X3), xi_3);
	return kernel_val;
}

double soap_kernel_Gaussian(double *X2_str, double *X3_str, double *X2_tr, double *X3_tr,
			 double beta_2, double beta_3, double xi_3, int size_X2, int size_X3) {

	double norm_X3_str, norm_X3_tr, X3_str_temp[size_X3], X3_tr_temp[size_X3], kernel_val;
	int i;

	double er_X2[size_X2], er_X3[size_X3];
	for (i=0; i <size_X2; i++)
		er_X2[i] = X2_str[i]-X2_tr[i];
	for (i=0; i <size_X3; i++)
		er_X3[i] = X3_str[i]-X3_tr[i];
	kernel_val = beta_2 * exp(-0.5*dotProduct(er_X2, er_X2, size_X2)) + beta_3 * exp(-0.5*dotProduct(er_X3, er_X3, size_X3));
	// kernel_val = beta_2 * dotProduct(X2_str, X2_tr, size_X2) + beta_3 * pow(dotProduct(X3_tr_temp, X3_str_temp, size_X3), xi_3);
	return kernel_val;
}

double soap_kernel_Laplacian(double *X2_str, double *X3_str, double *X2_tr, double *X3_tr,
			 double beta_2, double beta_3, double xi_3, int size_X2, int size_X3) {

	double norm_X3_str, norm_X3_tr, X3_str_temp[size_X3], X3_tr_temp[size_X3], kernel_val;
	int i;

	double er_X2[size_X2], er_X3[size_X3];
	for (i=0; i <size_X2; i++)
		er_X2[i] = X2_str[i]-X2_tr[i];
	for (i=0; i <size_X3; i++)
		er_X3[i] = X3_str[i]-X3_tr[i];
	kernel_val = beta_2 * exp(-0.5*sqrt(dotProduct(er_X2, er_X2, size_X2))) + beta_3 * exp(-0.5*sqrt(dotProduct(er_X3, er_X3, size_X3)));
	// kernel_val = beta_2 * dotProduct(X2_str, X2_tr, size_X2) + beta_3 * pow(dotProduct(X3_tr_temp, X3_str_temp, size_X3), xi_3);
	return kernel_val;
}

/*
der_soap_kernel function computes the derivative of the kernel w.r.t to some variable

[Input]
1. dX2_str: pointer to the derivative of first descriptor w.r.t to the given variable (2-body)
2. dX3_str: pointer to the derivative of first descriptor w.r.t to the given variable (3-body)
3. X2_str: pointer to the first descriptor (2-body)
4. X3_str: pointer to the first descriptor (3-body)
5. X2_tr: pointer to the second descriptor (2-body)
6. X3_tr: pointer to the second descriptor (3-body)
7. beta_2: weight to the 2-body term in the kernel
8. beta_3: weight to the 3-body term in the kernel
9. xi_3: exponent in the kernel
10. size_X2: length of the 2-body kernel
11. size_X3: length of the 3-body kernel
[Output]
1. der_val: derivative of the kernel
*/

double der_soap_kernel_Gaussian(double *dX2_str, double *dX3_str, double *X2_str, double *X3_str, double *X2_tr, double *X3_tr,
			 double beta_2, double beta_3, double xi_3, int size_X2, int size_X3) {

	double norm_X3_str, norm_X3_tr, const1, der_val, const2, X3_str_temp[size_X3], X3_tr_temp[size_X3], temp, temp1, temp0;



	double temp2[size_X2], temp3[size_X3];

	double er_X2[size_X2], er_X3[size_X3];
	for (int i=0; i <size_X2; i++)
		er_X2[i] = X2_str[i]-X2_tr[i];
	for (int i=0; i <size_X3; i++)
		er_X3[i] = X3_str[i]-X3_tr[i];

	double exp_temp3 = exp(-0.5*dotProduct(er_X3, er_X3, size_X3));
	double exp_temp2 = exp(-0.5*dotProduct(er_X2, er_X2, size_X2));

	for (int i=0; i <size_X2; i++)
		temp2[i] = -(X2_str[i]-X2_tr[i])*exp_temp2;
	for (int i=0; i <size_X3; i++)
		temp3[i] = -(X3_str[i]-X3_tr[i])*exp_temp3;

	der_val = beta_2 * dotProduct(temp2, dX2_str, size_X2) + beta_3 * dotProduct(temp3, dX3_str, size_X3);


	// temp = pow(norm_X3_str, -1-xi_3);

	// temp0 = dotProduct(X3_str, X3_tr_temp, size_X3);

	// temp1 = pow(temp0, xi_3-1);

	// const1 = -1 * beta_3* xi_3 * temp * temp1*temp0;
	// const2 = beta_3 * temp * norm_X3_str * xi_3 * temp1;

	// der_val = beta_2 * dotProduct(X2_tr, dX2_str, size_X2) + 
	// 		const1*dotProduct(X3_str_temp, dX3_str, size_X3) + const2*dotProduct(X3_tr_temp, dX3_str, size_X3);
	
	return der_val;
}

double der_soap_kernel_polynomial(double *dX2_str, double *dX3_str, double *X2_str, double *X3_str, double *X2_tr, double *X3_tr,
			 double beta_2, double beta_3, double xi_3, int size_X2, int size_X3) {

	double norm_X3_str, norm_X3_tr, const1, der_val, const2, X3_str_temp[size_X3], X3_tr_temp[size_X3], temp, temp1, temp0;

	norm_X3_str = sqrt(dotProduct(X3_str, X3_str, size_X3));
	norm_X3_tr = sqrt(dotProduct(X3_tr, X3_tr, size_X3));

	for (int i = 0; i<size_X3; i++){
		X3_str_temp[i] = X3_str[i]/norm_X3_str;
		X3_tr_temp[i] = X3_tr[i]/norm_X3_tr;
	}

	temp = pow(norm_X3_str, -1-xi_3);

	temp0 = dotProduct(X3_str, X3_tr_temp, size_X3);

	temp1 = pow(temp0, xi_3-1);

	const1 = -1 * beta_3* xi_3 * temp * temp1*temp0;
	const2 = beta_3 * temp * norm_X3_str * xi_3 * temp1;

	der_val = beta_2 * dotProduct(X2_tr, dX2_str, size_X2) + 
			const1*dotProduct(X3_str_temp, dX3_str, size_X3) + const2*dotProduct(X3_tr_temp, dX3_str, size_X3);
	
	return der_val;
}

double der_soap_kernel_Laplacian(double *dX2_str, double *dX3_str, double *X2_str, double *X3_str, double *X2_tr, double *X3_tr,
			 double beta_2, double beta_3, double xi_3, int size_X2, int size_X3) {

	double norm_X3_str, norm_X3_tr, const1, der_val, const2, X3_str_temp[size_X3], X3_tr_temp[size_X3], temp, temp1, temp0;


	double temp2[size_X2], temp3[size_X3];
	double er_X2[size_X2], er_X3[size_X3];
	for (int i=0; i <size_X2; i++)
		er_X2[i] = X2_str[i]-X2_tr[i];
	for (int i=0; i <size_X3; i++)
		er_X3[i] = X3_str[i]-X3_tr[i];

	double exp_temp2 = exp(-0.5*sqrt(dotProduct(er_X2, er_X2, size_X2)));
	double exp_temp3 = exp(-0.5*sqrt(dotProduct(er_X3, er_X3, size_X3)));

	for (int i=0; i <size_X2; i++){
		if ((X2_str[i]-X2_tr[i])>0)
			temp = -0.5;
		else if ((X2_str[i]-X2_tr[i])<0)
			temp = 0.5;
		else
			temp = 0.0;
		temp2[i] = temp*exp_temp2;
	}
	for (int i=0; i <size_X3; i++){
		if ((X3_str[i]-X3_tr[i])>0)
			temp = -0.5;
		else if ((X3_str[i]-X3_tr[i])<0)
			temp = 0.5;
		else
			temp = 0.0;
		temp3[i] = temp*exp_temp3;
	}

	der_val = beta_2 * dotProduct(temp2, dX2_str, size_X2) + beta_3 * dotProduct(temp3, dX3_str, size_X3);



	// temp = pow(norm_X3_str, -1-xi_3);

	// temp0 = dotProduct(X3_str, X3_tr_temp, size_X3);

	// temp1 = pow(temp0, xi_3-1);

	// const1 = -1 * beta_3* xi_3 * temp * temp1*temp0;
	// const2 = beta_3 * temp * norm_X3_str * xi_3 * temp1;

	// der_val = beta_2 * dotProduct(X2_tr, dX2_str, size_X2) + 
	// 		const1*dotProduct(X3_str_temp, dX3_str, size_X3) + const2*dotProduct(X3_tr_temp, dX3_str, size_X3);
	
	return der_val;
}


/*
init_MLFF function initializes and dynamically allocates memory to various objects in MLFF_Obj

[Input]
1. soap_str: SOAP structure to the first MD structure
2. mlff_str: MLFF_Obj structure to be initialized
3. n_str_max: max reference structure in training dataset
4. n_train_max: max local descriptor per element type in training dataset
[Output]
1. mlff_str: MLFF_Obj structure initialized
*/

void init_MLFF(MLFF_Obj *mlff_str, SPARC_OBJ *pSPARC) {

	int nelem = pSPARC->Ntypes, natom = pSPARC->n_atom, K_size_row, K_size_column, b_size, w_size,i, j;
	int size_X2, size_X3;
	int N_r = 72;

	// size_X2 = nelem * pSPARC->N_max_SOAP;
	// size_X3 = ((nelem * pSPARC->N_max_SOAP+1)*(nelem * pSPARC->N_max_SOAP))/2 * (pSPARC->L_max_SOAP+1);
	mlff_str->descriptor_typ = pSPARC->descriptor_typ_MLFF;
	if (pSPARC->descriptor_typ_MLFF==0){
		size_X2 = nelem * pSPARC->N_max_SOAP;
		size_X3 = ((nelem * pSPARC->N_max_SOAP+1)*(nelem * pSPARC->N_max_SOAP))/2 * (pSPARC->L_max_SOAP+1);
	}
	if (pSPARC->descriptor_typ_MLFF==1){
		size_X2 = pSPARC->N_max_SOAP;
		size_X3 = ((pSPARC->N_max_SOAP+1)*(pSPARC->N_max_SOAP))/2 * (pSPARC->L_max_SOAP+1);
	}
	


	mlff_str->n_str_max = pSPARC->n_str_max_mlff;
	mlff_str->n_train_max = pSPARC->n_train_max_mlff;
	mlff_str->n_str = 0;
	mlff_str->n_rows = 0;
	mlff_str->n_cols = 0;
	mlff_str->natm_train_total = 0;
	mlff_str->nelem = nelem;
	mlff_str->natom = natom;
	mlff_str->size_X2 = size_X2;
	mlff_str->size_X3 = size_X3;
	mlff_str->beta_2 = pSPARC->beta_2_SOAP;
	mlff_str->beta_3 = pSPARC->beta_3_SOAP;
	mlff_str->xi_3 = pSPARC->xi_3_SOAP;
	mlff_str->Nmax = pSPARC->N_max_SOAP;
	mlff_str->Lmax = pSPARC->L_max_SOAP;
	mlff_str->rcut = pSPARC->rcut_SOAP;
	mlff_str->F_tol = pSPARC->F_tol_SOAP;
	mlff_str->relative_scale_F = pSPARC->F_rel_scale;
	mlff_str->relative_scale_stress[0] = pSPARC->stress_rel_scale[0];
	mlff_str->relative_scale_stress[1] = pSPARC->stress_rel_scale[1];
	mlff_str->relative_scale_stress[2] = pSPARC->stress_rel_scale[2];
	mlff_str->relative_scale_stress[3] = pSPARC->stress_rel_scale[3];
	mlff_str->relative_scale_stress[4] = pSPARC->stress_rel_scale[4];
	mlff_str->relative_scale_stress[5] = pSPARC->stress_rel_scale[5];
	mlff_str->sigma_w = 100;
	mlff_str->sigma_v = 0.1;
	mlff_str->kernel_typ = pSPARC->kernel_typ_MLFF;
	
	mlff_str->rgrid = (double *) malloc(sizeof(double)* N_r);
 	mlff_str->h_nl = (double *) malloc(sizeof(double)* N_r* pSPARC->N_max_SOAP*(pSPARC->L_max_SOAP+1));
	mlff_str->dh_nl = (double *) malloc(sizeof(double)* N_r * pSPARC->N_max_SOAP*(pSPARC->L_max_SOAP+1));
	read_h_nl(pSPARC->N_max_SOAP, pSPARC->L_max_SOAP, mlff_str->rgrid, mlff_str->h_nl, mlff_str->dh_nl);


	mlff_str->cov_train = (double *) malloc(1*sizeof(double));
	mlff_str->natm_train_elemwise = (int *) malloc(nelem * sizeof(int));
	for (i=0; i<nelem; i++){
		mlff_str->natm_train_elemwise[i] = 0;
	}
	mlff_str->natm_typ_train = (int *)malloc(sizeof(int)*nelem * pSPARC->n_train_max_mlff);

	K_size_row = pSPARC->n_str_max_mlff * (3*natom + 1 + 6);
	K_size_column = nelem * pSPARC->n_train_max_mlff;
	b_size = pSPARC->n_str_max_mlff*(3*natom + 1 + 6);
	w_size = nelem * pSPARC->n_train_max_mlff;

	mlff_str->K_train = (double **) malloc(sizeof(double*)*K_size_row);
	for (i =0; i < K_size_row; i++){
		mlff_str->K_train[i] = (double *) malloc(sizeof(double)*K_size_column);
	}

	mlff_str->b_no_norm = (double *) malloc(sizeof(double)*b_size);
	mlff_str->weights = (double *) malloc(sizeof(double)*w_size);

	for (i = 0; i < K_size_row; i++)
		for (j=0; j < K_size_column; j++)
			mlff_str->K_train[i][j] = 0;

	for (i = 0; i < b_size; i++)
		mlff_str->b_no_norm[i] = 0;

	for (i = 0; i < w_size; i++)
		mlff_str->weights[i] = 0;

	mlff_str->soap_descriptor_strdataset = (SoapObj *) malloc(sizeof(SoapObj)*pSPARC->n_str_max_mlff);

	mlff_str->X2_traindataset = (double **) malloc(sizeof(double*)*pSPARC->n_train_max_mlff*nelem);
	mlff_str->X3_traindataset = (double **) malloc(sizeof(double*)*pSPARC->n_train_max_mlff*nelem);
	for (int i = 0; i < pSPARC->n_train_max_mlff*nelem; i++){
		mlff_str->X2_traindataset[i] = (double *) malloc(sizeof(double)*size_X2);
		mlff_str->X3_traindataset[i] = (double *) malloc(sizeof(double)*size_X3);
	}
}

/*
free_MLFF function frees the memory allocated to various objects in MLFF_Obj

[Input]
1. mlff_str: MLFF_Obj structure 
[Output]
1. mlff_str: MLFF_Obj structure 
*/

void free_MLFF(MLFF_Obj *mlff_str){

	int natom = mlff_str->natom, K_size_row = mlff_str->n_str_max*(3*natom + 1 + 6), nelem = mlff_str->nelem, i;
	for (i = 0; i < K_size_row; i++){
		free(mlff_str->K_train[i]);
	}
	free(mlff_str->K_train);
	free(mlff_str->b_no_norm);
	free(mlff_str->weights);
	free(mlff_str->natm_typ_train);
	free(mlff_str->natm_train_elemwise);
	free(mlff_str->cov_train);
	for (i=0; i < mlff_str->n_train_max*nelem; i++){
		free(mlff_str->X2_traindataset[i]);
		free(mlff_str->X3_traindataset[i]);
	}
	free(mlff_str->X2_traindataset);
	free(mlff_str->X3_traindataset);

	 for (i = 0; i < mlff_str->n_str; i++){
	 	delete_soapObj(mlff_str->soap_descriptor_strdataset+i);
	 }
	 free(mlff_str->soap_descriptor_strdataset);
}


/*
copy_descriptors function copies the content of one SoapObj to another

[Input]
1. soap_str: SoapObj structure to be copied
[Output]
1. soap_str_MLFF: SoapObj structure where it needs to be copied
*/

void copy_descriptors(SoapObj *soap_str_MLFF, SoapObj *soap_str){
	int i, j, k, l;
	soap_str_MLFF->natom = soap_str->natom;
	soap_str_MLFF->rcut = soap_str->rcut;
	soap_str_MLFF->cell[0] = soap_str->cell[0];
	soap_str_MLFF->cell[1] = soap_str->cell[1];
	soap_str_MLFF->cell[2] = soap_str->cell[2];
	soap_str_MLFF->Lmax = soap_str->Lmax;
	soap_str_MLFF->Nmax = soap_str->Nmax;
	soap_str_MLFF->size_X2 = soap_str->size_X2;
	soap_str_MLFF->size_X3 = soap_str->size_X3;
	soap_str_MLFF->beta_2 = soap_str->beta_2;
	soap_str_MLFF->beta_3 = soap_str->beta_3;
	soap_str_MLFF->xi_3 = soap_str->xi_3;
	soap_str_MLFF->nelem = soap_str->nelem;
	soap_str_MLFF->N_rgrid = soap_str->N_rgrid;

	for (i = 0; i < soap_str->natom; i++)
		soap_str_MLFF->Nneighbors[i] = soap_str->Nneighbors[i];

	for (i = 0; i < soap_str->natom; i++)
		soap_str_MLFF->unique_Nneighbors[i] = soap_str->unique_Nneighbors[i];

	for (i = 0; i < soap_str->nelem; i++)
		soap_str_MLFF->natom_elem[i] = soap_str->natom_elem[i];

	for (i = 0; i < soap_str->natom; i++)
		for (j = 0; j < soap_str->nelem; j++)
			soap_str_MLFF->unique_Nneighbors_elemWise[i][j] = soap_str->unique_Nneighbors_elemWise[i][j];

	for (i = 0; i < soap_str->natom; i++){
		soap_str_MLFF->neighborList[i].len = soap_str->neighborList[i].len;
		soap_str_MLFF->neighborList[i].array = (int*)malloc(sizeof(int)*soap_str->neighborList[i].len);
		soap_str_MLFF->unique_neighborList[i].len = soap_str->unique_neighborList[i].len;
		soap_str_MLFF->unique_neighborList[i].array = (int*)malloc(sizeof(int)*soap_str->unique_neighborList[i].len);
		for (l = 0; l < soap_str->neighborList[i].len; l++){
			soap_str_MLFF->neighborList[i].array[l] = soap_str->neighborList[i].array[l];
		}
		for (l = 0; l < soap_str->unique_neighborList[i].len; l++){
			soap_str_MLFF->unique_neighborList[i].array[l] = soap_str->unique_neighborList[i].array[l];
		}
		for (j = 0; j < soap_str->nelem; j++){
			soap_str_MLFF->unique_neighborList_elemWise[i][j].len = soap_str->unique_neighborList_elemWise[i][j].len;
			soap_str_MLFF->unique_neighborList_elemWise[i][j].array = (int*)malloc(sizeof(int)*soap_str->unique_neighborList_elemWise[i][j].len);
			for (k = 0; k < soap_str->unique_neighborList_elemWise[i][j].len; k++){
				soap_str_MLFF->unique_neighborList_elemWise[i][j].array[k] = soap_str->unique_neighborList_elemWise[i][j].array[k];
			}
		}
	}


	for (i = 0; i < soap_str->natom; i++){
		for (j = 0; j < soap_str->size_X2; j++){
			soap_str_MLFF->X2[i][j] = soap_str->X2[i][j];
		}
		for (j = 0; j < soap_str->size_X3; j++){
			soap_str_MLFF->X3[i][j] = soap_str->X3[i][j];
		}
		int uniq_natms = uniqueEle((soap_str->neighborList[i]).array, soap_str->Nneighbors[i]);
		for (j = 0; j < 1+uniq_natms; j++){
			for (k = 0; k < soap_str->size_X2; k++){
				soap_str_MLFF->dX2_dX[i][j][k] = soap_str->dX2_dX[i][j][k];
				soap_str_MLFF->dX2_dY[i][j][k] = soap_str->dX2_dY[i][j][k];
				soap_str_MLFF->dX2_dZ[i][j][k] = soap_str->dX2_dZ[i][j][k];
			}
			for (k = 0; k < soap_str->size_X3; k++){
				soap_str_MLFF->dX3_dX[i][j][k] = soap_str->dX3_dX[i][j][k];
				soap_str_MLFF->dX3_dY[i][j][k] = soap_str->dX3_dY[i][j][k];
				soap_str_MLFF->dX3_dZ[i][j][k] = soap_str->dX3_dZ[i][j][k];
			}
		}
		for (j = 0; j < 6; j++){
			for (k = 0; k < soap_str->size_X2; k++){
				soap_str_MLFF->dX2_dF[i][j][k] = soap_str->dX2_dF[i][j][k];
			}
			for (k = 0; k < soap_str->size_X3; k++){
				soap_str_MLFF->dX3_dF[i][j][k] = soap_str->dX3_dF[i][j][k];
			}
		}

	}
}

/*
add_firstMD function updates the MLFF_Obj by updating design matrix, b vector etc. for the first MD

[Input]
1. soap_str: SoapObj structure of the first MD
2. nlist: NeighList strcuture of the first MD
3. mlff_str: MLFF_Obj structure
4. E: energy per atom of first MD structure (Ha/atom)
5. F: atomic foces of first MD structure (Ha/bohr) [ColMajor]
6. stress: stress of first MD structure (GPa)
[Output]
1. mlff_str: MLFF_Obj structure
*/


void add_firstMD(SoapObj *soap_str, NeighList *nlist, MLFF_Obj *mlff_str, double E, double* F, double* stress) {
	int row_idx, col_idx, count, atom_idx, size_X2 = soap_str->size_X2, size_X3 = soap_str->size_X3;
	int natom = soap_str->natom, nelem = soap_str->nelem, i, j, iel, iatm_el_tr, iatm_str, neighs, istress;
	double beta_2 = soap_str->beta_2, beta_3 = soap_str->beta_3, xi_3 = soap_str->xi_3;
	// for (int i = 0; i<3*natom; i++)
	// 	F[i] = -F[i];
	// double *stress;
	// stress = (double *)malloc(6*sizeof(double));
	// stress[0] = stress1[0];
	// stress[1] = stress1[3];
	// stress[2] = stress1[5];
	// stress[3] = stress1[1];
	// stress[4] = stress1[4];
	// stress[5] = stress1[2];
// calculate mean and std deviation to do the normalization

	int kernel_typ = mlff_str->kernel_typ;
	mlff_str->mu_E  = E;
	mlff_str->std_E  = 1;
	mlff_str->std_F  = sqrt(get_variance(F, 3*natom));
	for (int i = 0; i < 6; i++){
		mlff_str->mu_stress[i]  = au2GPa*stress[i];
		mlff_str->std_stress[i]  = 1;
	}
// populate bvec
	mlff_str->b_no_norm[0] = E;

	for (i = 0; i < natom; i++){
		mlff_str->b_no_norm[3*i+1] = F[3*i];
		mlff_str->b_no_norm[3*i+2] = F[3*i+1];
		mlff_str->b_no_norm[3*i+3] = F[3*i+2];
	}
	for (i = 0; i < 6; i++){
		mlff_str->b_no_norm[1+3*natom+i] = au2GPa*stress[i];
	}

	double temp = (1.0/ (double) soap_str->natom);
	int *cum_natm_elem, *cum_natm_elem1;
	cum_natm_elem = (int *)malloc(sizeof(int)*nelem);
	cum_natm_elem1 = (int *)malloc(sizeof(int)*nelem);
	cum_natm_elem[0] = 0;
	for (int i = 1; i < nelem; i++){
		cum_natm_elem[i] = soap_str->natom_elem[i-1];
	}



	dyArray *highrank_ID_descriptors;
	highrank_ID_descriptors = (dyArray *) malloc(sizeof(dyArray)*nelem);
	for (int i=0; i <nelem; i++){
		init_dyarray(&highrank_ID_descriptors[i]);
		SOAP_CUR_sparsify(kernel_typ, &soap_str->X2[cum_natm_elem[i]], &soap_str->X3[cum_natm_elem[i]],
				 soap_str->natom_elem[i], size_X2, size_X3, beta_2, beta_3, xi_3, &highrank_ID_descriptors[i]);
	}
	
	for (int i = 0; i < nelem; i++){
		// mlff_str->natm_train_elemwise[i] = soap_str->natom_elem[i];
		mlff_str->natm_train_elemwise[i] = (highrank_ID_descriptors[i]).len;
	}

	count=0;
	for (int i = 0; i < nelem; i++){
		// for (int j = 0; j < soap_str->natom_elem[i]; j++){
		for (int j = 0; j < (highrank_ID_descriptors[i]).len; j++){
			mlff_str->natm_typ_train[count] = i;
			for(int jj = 0; jj < size_X2; jj++){
				// mlff_str->X2_traindataset[count][jj] = soap_str->X2[count][jj];
				mlff_str->X2_traindataset[count][jj] = soap_str->X2[cum_natm_elem[i]+(highrank_ID_descriptors[i]).array[j]][jj];
			}
			for(int jj = 0; jj < size_X3; jj++){
				// mlff_str->X3_traindataset[count][jj] = soap_str->X3[count][jj];
				mlff_str->X3_traindataset[count][jj] = soap_str->X3[cum_natm_elem[i]+(highrank_ID_descriptors[i]).array[j]][jj];
			}
			count++;
		}
	}


// // copy the X2 and X3 into train history data 
// 	for (i = 0; i < natom; i++){
// 		for(j = 0; j < size_X2; j++){
// 			mlff_str->X2_traindataset[i][j] = soap_str->X2[i][j];
// 		}
// 		for(j = 0; j < size_X3; j++){
// 			mlff_str->X3_traindataset[i][j] = soap_str->X3[i][j];
// 		}
// 	}

// copy X2, X3, dX2, dX3 for the first structure into the mlff obj 
	initialize_soapObj(mlff_str->soap_descriptor_strdataset, nlist, soap_str->Lmax, soap_str->Nmax, soap_str->beta_3, soap_str->xi_3);
	copy_descriptors(mlff_str->soap_descriptor_strdataset, soap_str);

// updating other MLFF parameters such as number of structures, number of training environment, element typ of training env
	mlff_str->n_str = 1;
	mlff_str->natm_train_total = count;
	// mlff_str->natm_train_total = soap_str->natom;

	mlff_str->n_rows = 3*soap_str->natom + 1 + 6;
	// mlff_str->n_cols = soap_str->natom;
	mlff_str->n_cols = count;


// Energy term to populate K_train matrix
	cum_natm_elem1[0] = 0;
	for (int i = 1; i < nelem; i++){
		cum_natm_elem1[i] = soap_str->natom_elem[i-1];
	}

	for (iel = 0; iel < nelem; iel++){
		// for (iatm_el_tr = 0; iatm_el_tr < soap_str->natom_elem[iel]; iatm_el_tr++){
		for (iatm_el_tr = 0; iatm_el_tr < (highrank_ID_descriptors[iel]).len; iatm_el_tr++){
			col_idx = cum_natm_elem1[iel] + iatm_el_tr;
			// if (iel > 0)
			// 	col_idx = soap_str->natom_elem[iel-1] + iatm_el_tr;
			for (iatm_str = 0; iatm_str < soap_str->natom_elem[iel]; iatm_str++){
				atom_idx = cum_natm_elem[iel] + iatm_str;
				// if (iel>0)
				// 	atom_idx = soap_str->natom_elem[iel-1] + iatm_str;
				// mlff_str->K_train[0][col_idx] +=  temp * soap_kernel(soap_str->X2[atom_idx], soap_str->X3[atom_idx],
				// 										 soap_str->X2[col_idx], soap_str->X3[col_idx],
				// 								 		 beta_2, beta_3, xi_3, size_X2, size_X3);
				mlff_str->K_train[0][col_idx] +=  temp * soap_kernel(kernel_typ, soap_str->X2[atom_idx], soap_str->X3[atom_idx],
														 mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
												 		 beta_2, beta_3, xi_3, size_X2, size_X3);			
			}
		}
	}


// Force terms to populate K_train matrix
	for (iel = 0; iel < nelem; iel++){
		// for (iatm_el_tr = 0; iatm_el_tr < soap_str->natom_elem[iel]; iatm_el_tr++){
		for (iatm_el_tr = 0; iatm_el_tr < (highrank_ID_descriptors[iel]).len; iatm_el_tr++){
			col_idx = cum_natm_elem1[iel] + iatm_el_tr;
			// if (iel > 0)
			// 	col_idx = soap_str->natom_elem[iel-1] + iatm_el_tr;
			for (iatm_str = 0; iatm_str < soap_str->natom_elem[iel]; iatm_str++){
				atom_idx = cum_natm_elem[iel] + iatm_str;
				// if (iel > 0)
				// 	atom_idx = soap_str->natom_elem[iel-1]+iatm_str;

				row_idx = 3*atom_idx+1;
				// x-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
				// mlff_str->K_train[row_idx][col_idx] +=  
				// 	der_soap_kernel(soap_str->dX2_dX[atom_idx][0], soap_str->dX3_dX[atom_idx][0],
				// 	 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
				// 	soap_str->X2[col_idx], soap_str->X3[col_idx],
				// 	beta_2, beta_3, xi_3, size_X2, size_X3);
				mlff_str->K_train[row_idx][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dX[atom_idx][0], soap_str->dX3_dX[atom_idx][0],
					 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);
				// y-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
				// mlff_str->K_train[row_idx+1][col_idx] +=  
				// 	der_soap_kernel(soap_str->dX2_dY[atom_idx][0], soap_str->dX3_dY[atom_idx][0],
				// 	 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
				// 	soap_str->X2[col_idx], soap_str->X3[col_idx],
				// 	beta_2, beta_3, xi_3, size_X2, size_X3);
					mlff_str->K_train[row_idx+1][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dY[atom_idx][0], soap_str->dX3_dY[atom_idx][0],
					 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);
				// z-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
				// mlff_str->K_train[row_idx+2][col_idx] +=  
				// 	der_soap_kernel(soap_str->dX2_dZ[atom_idx][0], soap_str->dX3_dZ[atom_idx][0],
				// 	 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
				// 	soap_str->X2[col_idx], soap_str->X3[col_idx],
				// 	beta_2, beta_3, xi_3, size_X2, size_X3);
					mlff_str->K_train[row_idx+2][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dZ[atom_idx][0], soap_str->dX3_dZ[atom_idx][0],
					 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);

				for (neighs =0; neighs < soap_str->unique_Nneighbors[atom_idx]; neighs++){
					row_idx = 3*soap_str->unique_neighborList[atom_idx].array[neighs]+1; // Possible source of error
					// x-component (w.r.t neighs neighbour)
					// mlff_str->K_train[row_idx][col_idx] +=  
					// der_soap_kernel(soap_str->dX2_dX[atom_idx][1+neighs], soap_str->dX3_dX[atom_idx][1+neighs],
					// 	 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					// 	soap_str->X2[col_idx], soap_str->X3[col_idx],
					// 	beta_2, beta_3, xi_3, size_X2, size_X3);
					mlff_str->K_train[row_idx][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dX[atom_idx][1+neighs], soap_str->dX3_dX[atom_idx][1+neighs],
						 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
						beta_2, beta_3, xi_3, size_X2, size_X3);
					// y-component (w.r.t neighs neighbour)
					// mlff_str->K_train[row_idx+1][col_idx] +=  
					// der_soap_kernel(soap_str->dX2_dY[atom_idx][1+neighs], soap_str->dX3_dY[atom_idx][1+neighs], 
					// 	soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					// 	soap_str->X2[col_idx], soap_str->X3[col_idx],
					// 	beta_2, beta_3, xi_3, size_X2, size_X3);
					mlff_str->K_train[row_idx+1][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dY[atom_idx][1+neighs], soap_str->dX3_dY[atom_idx][1+neighs], 
						soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
						beta_2, beta_3, xi_3, size_X2, size_X3);
					// z-component (w.r.t neighs neighbour)
					// mlff_str->K_train[row_idx+2][col_idx] +=  
					// der_soap_kernel(soap_str->dX2_dZ[atom_idx][1+neighs], soap_str->dX3_dZ[atom_idx][1+neighs], 
					// 	soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					// 	soap_str->X2[col_idx], soap_str->X3[col_idx],
					// 	beta_2, beta_3, xi_3, size_X2, size_X3);
					mlff_str->K_train[row_idx+2][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dZ[atom_idx][1+neighs], soap_str->dX3_dZ[atom_idx][1+neighs], 
						soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
						beta_2, beta_3, xi_3, size_X2, size_X3);
				}

			}
		}
	}

	double volume = soap_str->cell[0] *soap_str->cell[1]*soap_str->cell[2];

// Stress terms to populate K_train matrix

	for (iel = 0; iel < nelem; iel++){
		// for (int iatm_el_tr = 0; iatm_el_tr < soap_str->natom_elem[iel]; iatm_el_tr++){
		for (iatm_el_tr = 0; iatm_el_tr < (highrank_ID_descriptors[iel]).len; iatm_el_tr++){
			col_idx = cum_natm_elem1[iel] + iatm_el_tr;
			// if (iel > 0)
			// 	col_idx = soap_str->natom_elem[iel-1]+iatm_el_tr;
			for (iatm_str = 0; iatm_str < soap_str->natom_elem[iel]; iatm_str++){
				atom_idx = cum_natm_elem[iel] + iatm_str;
				// if (iel > 0)
				// 	atom_idx = soap_str->natom_elem[iel-1] + iatm_str;
				for (istress =0; istress < 6; istress++){
					row_idx = 3*soap_str->natom+1+istress;

					// mlff_str->K_train[row_idx][col_idx] +=  au2GPa*(1.0/volume)*
					// der_soap_kernel(soap_str->dX2_dF[atom_idx][istress], soap_str->dX3_dF[atom_idx][istress],
					// 	 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					// 	soap_str->X2[col_idx], soap_str->X3[col_idx],
					// 	beta_2, beta_3, xi_3, size_X2, size_X3);
					mlff_str->K_train[row_idx][col_idx] +=  au2GPa*(1.0/volume)*
					der_soap_kernel(kernel_typ, soap_str->dX2_dF[atom_idx][istress], soap_str->dX3_dF[atom_idx][istress],
						 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
						beta_2, beta_3, xi_3, size_X2, size_X3);

				}

			}
		}
	}

	free(cum_natm_elem);
	free(cum_natm_elem1);
	// free(stress);

}

/*
add_newstr_rows function updates the MLFF_Obj by updating design matrix, b vector etc. for a new reference structure

[Input]
1. soap_str: SoapObj structure of the reference structure to be added
2. nlist: NeighList strcuture of the first MD
3. mlff_str: MLFF_Obj structure
4. E: energy per atom of the reference structure to be added (Ha/atom)
5. F: atomic foces of the reference structure to be added (Ha/bohr) [ColMajor]
6. stress: stress of the reference structure to be added (GPa)
[Output]
1. mlff_str: MLFF_Obj structure
*/

void add_newstr_rows(SoapObj *soap_str, NeighList *nlist, MLFF_Obj *mlff_str, double E, double *F, double *stress) {
	int row_idx, col_idx, atom_idx, natom = soap_str->natom, nelem = soap_str->nelem, num_Fterms_exist, num_Fterms_newstr, i, iel, iatm_el_tr, iatm_str, neighs, istress;
	int size_X2 = soap_str->size_X2, size_X3 = soap_str->size_X3;
	double beta_2 = soap_str->beta_2, beta_3 = soap_str->beta_3, xi_3 = soap_str->xi_3;
	double old_mu_stress, old_std_stress, old_mu_E, old_std_E, old_std_F, newstr_std_F;

	int kernel_typ = mlff_str->kernel_typ;
	// double *stress;
	// stress = (double *)malloc(6*sizeof(double));
	// stress[0] = stress1[0];
	// stress[1] = stress1[3];
	// stress[2] = stress1[5];
	// stress[3] = stress1[1];
	// stress[4] = stress1[4];
	// stress[5] = stress1[2];


	// for (int i = 0; i<3*natom; i++)
	// 	F[i] = -F[i];
	// Updating the mean and std deviation of E, F, stress to do the normalization

	old_mu_E = mlff_str->mu_E;
	old_std_E = mlff_str->std_E;
	old_std_F = mlff_str->std_F;
	mlff_str->mu_E  = old_mu_E *(mlff_str->n_str/(mlff_str->n_str+1.0)) + E/(mlff_str->n_str+1.0);

	// printf("mlff_str->n_str: %d",mlff_str->n_str);
	if (mlff_str->n_str==1){
		double E_mean = 0.5*(old_mu_E+E);
		mlff_str->std_E = sqrt(((E-E_mean)*(E-E_mean) +(old_mu_E-E_mean)*(old_mu_E-E_mean))/mlff_str->n_str);

	} else if (mlff_str->n_str==2) {
		double E0 = mlff_str->b_no_norm[0], E1 = mlff_str->b_no_norm[7+3*natom], E2 = E;
		double E_mean = (1.0/3.0)*(E0+E1+E2);
		mlff_str->std_E = sqrt((E0-E_mean)*(E0-E_mean)+(E1-E_mean)*(E1-E_mean)+(E2-E_mean)*(E2-E_mean)/2.0);
	}else {

		mlff_str->std_E = sqrt(((mlff_str->n_str-1.0)*old_std_E*old_std_E + (E-mlff_str->mu_E)*(E-old_mu_E))/mlff_str->n_str);
	}

	newstr_std_F = sqrt(get_variance(F, 3*natom));
	num_Fterms_exist = mlff_str->n_rows - mlff_str->n_str*7;
	num_Fterms_newstr = 3*soap_str->natom;

	mlff_str->std_F = sqrt((double) (num_Fterms_exist-1)/(double) (num_Fterms_exist+num_Fterms_newstr-1.0) *(old_std_F*old_std_F)
						+ (double) (num_Fterms_newstr-1)/(double) (num_Fterms_exist+num_Fterms_newstr-1.0) * (newstr_std_F*newstr_std_F));

	for (i=0; i<6; i++){
		old_mu_stress = mlff_str->mu_stress[i];
		old_std_stress = mlff_str->std_stress[i];
		mlff_str->mu_stress[i]  = old_mu_stress*(mlff_str->n_str/(mlff_str->n_str+1.0)) + au2GPa*stress[i]/(double) (mlff_str->n_str+1.0);
		if (mlff_str->n_str==1){
			double mean_stress = 0.5*(au2GPa*stress[i]+old_mu_stress);
			mlff_str->std_stress[i] = sqrt(((au2GPa*stress[i]-mean_stress)*(au2GPa*stress[i]-mean_stress) + (old_mu_stress-mean_stress)*(old_mu_stress-mean_stress))/(double)mlff_str->n_str);
		} else if (mlff_str->n_str==2){
			double st0 = mlff_str->b_no_norm[1+3*natom+i], st1 = mlff_str->b_no_norm[7+3*natom+1+3*natom+i], st2 =au2GPa*stress[i];
			double mean_stress = (1.0/3.0)*(st0+st1+st2);
			mlff_str->std_stress[i] = sqrt((st0-mean_stress)*(st0-mean_stress)+(st1-mean_stress)*(st1-mean_stress)+(st2-mean_stress)*(st2-mean_stress)/2.0);
		}else {
			mlff_str->std_stress[i] = sqrt(((double) (mlff_str->n_str-1.0)*old_std_stress*old_std_stress 
				+ (au2GPa*stress[i]-mlff_str->mu_stress[i])*(au2GPa*stress[i]-old_mu_stress))/(double) mlff_str->n_str);
		}
	}

	printf("std_E: %f, std_F: %f\n",mlff_str->std_E,mlff_str->std_F);
	printf("std_stress: %f,%f,%f,%f,%f,%f\n",mlff_str->std_stress[0],mlff_str->std_stress[1],mlff_str->std_stress[2],
		mlff_str->std_stress[3],mlff_str->std_stress[4],mlff_str->std_stress[5]);
	// populate bvec (not normalized)
	mlff_str->b_no_norm[mlff_str->n_rows] = E;
	for (i = 0; i < natom; i++){
		mlff_str->b_no_norm[mlff_str->n_rows+3*i+1] = F[3*i];
		mlff_str->b_no_norm[mlff_str->n_rows+3*i+2] = F[3*i+1];
		mlff_str->b_no_norm[mlff_str->n_rows+3*i+3] = F[3*i+2];
	}
	for (i = 0; i < 6; i++){
		mlff_str->b_no_norm[1 + mlff_str->n_rows + 3*natom + i] = au2GPa*stress[i];
	}

	// copy the X2 and X3 into train history data 
		// Not required in "add_str" - only do in "add-train"

	// copy X2, X3, dX2, dX3 for the first structure into the mlff obj 
	initialize_soapObj(mlff_str->soap_descriptor_strdataset + mlff_str->n_str, nlist, soap_str->Lmax, soap_str->Nmax, soap_str->beta_3, soap_str->xi_3);
	copy_descriptors(mlff_str->soap_descriptor_strdataset+mlff_str->n_str, soap_str);

	int *cum_natm_elem;
	cum_natm_elem = (int *)malloc(sizeof(int)*nelem);
	cum_natm_elem[0] = 0;
	for (int i = 1; i < nelem; i++){
		cum_natm_elem[i] = soap_str->natom_elem[i-1];
	}

	int **cum_natm_ele_cols;
	cum_natm_ele_cols = (int **)malloc(sizeof(int*)*nelem);
	for (int i = 0; i <nelem; i++){
		cum_natm_ele_cols[i] = (int *)malloc(sizeof(int)*mlff_str->natm_train_elemwise[i]);
	}

	for (int i = 0; i < nelem; i++){
		int count=0;
		for (int j=0; j < mlff_str->natm_train_total; j++){
			if (mlff_str->natm_typ_train[j] == i){
				cum_natm_ele_cols[i][count] = j;
				count++;
			}
		}
	}

	// Energy term to populate K_train matrix
	double temp = (1.0/ (double) natom);
	for (iel = 0; iel < nelem; iel++){
		for (iatm_el_tr = 0; iatm_el_tr < mlff_str->natm_train_elemwise[iel]; iatm_el_tr++){
			// col_idx = cum_natm_elem[iel] + iatm_el_tr;
			col_idx = cum_natm_ele_cols[iel][iatm_el_tr];
			// if (iel > 0)
			// 	col_idx = mlff_str->natm_train_elemwise[iel-1] + iatm_el_tr;
			for (iatm_str = 0; iatm_str < soap_str->natom_elem[iel]; iatm_str++){
				atom_idx = cum_natm_elem[iel] + iatm_str;
				// if (iel > 0)
				// 	atom_idx = soap_str->natom_elem[iel-1] + iatm_str;
				mlff_str->K_train[mlff_str->n_rows][col_idx] +=  temp* soap_kernel(kernel_typ, soap_str->X2[atom_idx], soap_str->X3[atom_idx],
															mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
												 			beta_2, beta_3, xi_3, size_X2, size_X3);
			}
		}
	}

	// Force terms to populate K_train matrix

	for (iel = 0; iel < nelem; iel++){
		for (iatm_el_tr = 0; iatm_el_tr < mlff_str->natm_train_elemwise[iel]; iatm_el_tr++){
			// col_idx = cum_natm_elem[iel] + iatm_el_tr;
			col_idx = cum_natm_ele_cols[iel][iatm_el_tr];
			// if (iel > 0)
			// 	col_idx = mlff_str->natm_train_elemwise[iel-1] + iatm_el_tr;
			for (iatm_str = 0; iatm_str < soap_str->natom_elem[iel]; iatm_str++){
				atom_idx = cum_natm_elem[iel] + iatm_str;
				// if (iel>0)
				// 	atom_idx = soap_str->natom_elem[iel-1] + iatm_str;

				row_idx = mlff_str->n_rows + 3*atom_idx+1;
				// x-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
				mlff_str->K_train[row_idx][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dX[atom_idx][0], soap_str->dX3_dX[atom_idx][0],
					 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);
				// y-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
				mlff_str->K_train[row_idx+1][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dY[atom_idx][0], soap_str->dX3_dY[atom_idx][0],
					 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);
				// z-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
				mlff_str->K_train[row_idx+2][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dZ[atom_idx][0], soap_str->dX3_dZ[atom_idx][0],
					 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					mlff_str->X3_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);

				for (neighs =0; neighs < soap_str->unique_Nneighbors[atom_idx]; neighs++){
					row_idx = mlff_str->n_rows + 3*soap_str->unique_neighborList[atom_idx].array[neighs] + 1; // Possible source of error
					// x-component (w.r.t neighs neighbour)
					mlff_str->K_train[row_idx][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dX[atom_idx][1+neighs], soap_str->dX3_dX[atom_idx][1+neighs],
						 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
						beta_2, beta_3, xi_3, size_X2, size_X3);
					// y-component (w.r.t neighs neighbour)
					mlff_str->K_train[row_idx+1][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dY[atom_idx][1+neighs], soap_str->dX3_dY[atom_idx][1+neighs], 
						soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
						beta_2, beta_3, xi_3, size_X2, size_X3);
					// z-component (w.r.t neighs neighbour)
					mlff_str->K_train[row_idx+2][col_idx] +=  
					der_soap_kernel(kernel_typ, soap_str->dX2_dZ[atom_idx][1+neighs], soap_str->dX3_dZ[atom_idx][1+neighs], 
						soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
						beta_2, beta_3, xi_3, size_X2, size_X3);
				}

			}
		}
	}


	// Stress terms to populate K_train matrix
	double volume = soap_str->cell[0] *soap_str->cell[1]*soap_str->cell[2];
	for (iel = 0; iel < nelem; iel++){
		for (iatm_el_tr = 0; iatm_el_tr < mlff_str->natm_train_elemwise[iel]; iatm_el_tr++){
			// col_idx = cum_natm_elem[iel] + iatm_el_tr;
			col_idx = cum_natm_ele_cols[iel][iatm_el_tr];
			// if (iel > 0)
			// 	col_idx = mlff_str->natm_train_elemwise[iel-1] + iatm_el_tr;
			for (iatm_str = 0; iatm_str < soap_str->natom_elem[iel]; iatm_str++){
				atom_idx = cum_natm_elem[iel] + iatm_str;
				// if (iel>0)
				// 	atom_idx = soap_str->natom_elem[iel-1] + iatm_str;
				for (istress =0; istress < 6; istress++){
					row_idx = mlff_str->n_rows+1+3*soap_str->natom+istress;

					mlff_str->K_train[row_idx][col_idx] +=  au2GPa*(1.0/volume)*
					der_soap_kernel(kernel_typ, soap_str->dX2_dF[atom_idx][istress], soap_str->dX3_dF[atom_idx][istress],
						 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
						beta_2, beta_3, xi_3, size_X2, size_X3);

				}

			}
		}
	}
	// updating other MLFF parameters such as number of structures, number of training environment, element typ of training env
	mlff_str->n_str += 1;
	mlff_str->n_rows += 3*soap_str->natom + 1 + 6;
	free(cum_natm_elem);
	for (int i = 0; i <nelem; i++){
		free(cum_natm_ele_cols[i]);
	}
	free(cum_natm_ele_cols);
	// free(stress);

}


/*
calculate_Kpredict function calculate the design matrix for prediction for a new structure

[Input]
1. soap_str: SoapObj structure of the new structure
2. nlist: NeighList strcuture of the new structure
3. mlff_str: MLFF_Obj structure
[Output]
1. K_predict: design prediction matrix
*/

void calculate_Kpredict(SoapObj *soap_str, NeighList *nlist, MLFF_Obj *mlff_str, double **K_predict){

 	int row_idx, col_idx, atom_idx, iel, iatm_el_tr, iatm_str, istress, natom = soap_str->natom, nelem = soap_str->nelem;
	int size_X2 = soap_str->size_X2, size_X3 = soap_str->size_X3;
	double beta_2 = soap_str->beta_2, beta_3 = soap_str->beta_3, xi_3 = soap_str->xi_3, E_scale, F_scale, stress_scale[6] ;

 	//int nrows = 3*soap_str->natom + 1 + 6;
 	//int ncols = mlff_str->n_cols;
 	int kernel_typ = mlff_str->kernel_typ;

 	E_scale = mlff_str->E_scale;
 	F_scale = mlff_str->F_scale * mlff_str->relative_scale_F;
 	stress_scale[0] = mlff_str->stress_scale[0] * mlff_str->relative_scale_stress[0];
 	stress_scale[1] = mlff_str->stress_scale[1] * mlff_str->relative_scale_stress[1];
 	stress_scale[2] = mlff_str->stress_scale[2] * mlff_str->relative_scale_stress[2];
 	stress_scale[3] = mlff_str->stress_scale[3] * mlff_str->relative_scale_stress[3];
 	stress_scale[4] = mlff_str->stress_scale[4] * mlff_str->relative_scale_stress[4];
 	stress_scale[5] = mlff_str->stress_scale[5] * mlff_str->relative_scale_stress[5];

 	int *cum_natm_elem;
	cum_natm_elem = (int *)malloc(sizeof(int)*nelem);
	cum_natm_elem[0] = 0;
	for (int i = 1; i < nelem; i++){
		cum_natm_elem[i] = soap_str->natom_elem[i-1];
	}

	int **cum_natm_ele_cols;
	cum_natm_ele_cols = (int **)malloc(sizeof(int*)*nelem);
	for (int i = 0; i <nelem; i++){
		cum_natm_ele_cols[i] = (int *)malloc(sizeof(int)*mlff_str->natm_train_elemwise[i]);
	}

	for (int i = 0; i < nelem; i++){
		int count=0;
		for (int j=0; j < mlff_str->natm_train_total; j++){
			if (mlff_str->natm_typ_train[j] == i){
				cum_natm_ele_cols[i][count] = j;
				count++;
			}
		}
	}


 	// Energy term to populate K_predict matrix
	for (iel = 0; iel < nelem; iel++){
		for (iatm_el_tr = 0; iatm_el_tr < mlff_str->natm_train_elemwise[iel]; iatm_el_tr++){
			// col_idx = cum_natm_elem[iel] + iatm_el_tr;
			col_idx = cum_natm_ele_cols[iel][iatm_el_tr];
			// if (iel > 0)
			// 	col_idx = mlff_str->natm_train_elemwise[iel-1] + iatm_el_tr;
			for (iatm_str = 0; iatm_str < soap_str->natom_elem[iel]; iatm_str++){
				atom_idx = cum_natm_elem[iel] + iatm_str;
				// if (iel > 0)
				// 	atom_idx = soap_str->natom_elem[iel-1] + iatm_str;
				K_predict[0][col_idx] +=  E_scale*(1.0/soap_str->natom)*
												soap_kernel(kernel_typ, soap_str->X2[atom_idx], soap_str->X3[atom_idx],
												mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
												 beta_2, beta_3, xi_3, size_X2, size_X3);
			}
		}
	}


	// Force terms to populate K_predict matrix

	for (iel = 0; iel < nelem; iel++){
		for (iatm_el_tr = 0; iatm_el_tr < mlff_str->natm_train_elemwise[iel]; iatm_el_tr++){
			// col_idx = cum_natm_elem[iel] + iatm_el_tr;
			col_idx = cum_natm_ele_cols[iel][iatm_el_tr];
			// if (iel > 0)
			// 	col_idx = mlff_str->natm_train_elemwise[iel-1] + iatm_el_tr;
			for (iatm_str = 0; iatm_str < soap_str->natom_elem[iel]; iatm_str++){
				atom_idx = cum_natm_elem[iel] + iatm_str;
				// if (iel > 0)
				// 	atom_idx = soap_str->natom_elem[iel-1]+iatm_str;

				row_idx = 3*atom_idx+1;
				// x-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
				K_predict[row_idx][col_idx] +=  
					F_scale*der_soap_kernel(kernel_typ, soap_str->dX2_dX[atom_idx][0], soap_str->dX3_dX[atom_idx][0],
					 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);
				// y-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
				K_predict[row_idx+1][col_idx] +=  
					F_scale*der_soap_kernel(kernel_typ, soap_str->dX2_dY[atom_idx][0], soap_str->dX3_dY[atom_idx][0],
					 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);
				// z-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
				K_predict[row_idx+2][col_idx] +=  
					F_scale*der_soap_kernel(kernel_typ, soap_str->dX2_dZ[atom_idx][0], soap_str->dX3_dZ[atom_idx][0],
					 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
					mlff_str->X3_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);

				for (int neighs =0; neighs < soap_str->unique_Nneighbors[atom_idx]; neighs++){
					row_idx = 3*soap_str->unique_neighborList[atom_idx].array[neighs]+1; // Possible source of error
					// x-component (w.r.t neighs neighbour)
					K_predict[row_idx][col_idx] +=  
					F_scale*der_soap_kernel(kernel_typ, soap_str->dX2_dX[atom_idx][1+neighs], soap_str->dX3_dX[atom_idx][1+neighs],
						 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);
					// y-component (w.r.t neighs neighbour)
					K_predict[row_idx+1][col_idx] +=  
					F_scale*der_soap_kernel(kernel_typ, soap_str->dX2_dY[atom_idx][1+neighs], soap_str->dX3_dY[atom_idx][1+neighs], 
						soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);
					// z-component (w.r.t neighs neighbour)
					K_predict[row_idx+2][col_idx] +=  
					F_scale*der_soap_kernel(kernel_typ, soap_str->dX2_dZ[atom_idx][1+neighs], soap_str->dX3_dZ[atom_idx][1+neighs], 
						soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
					beta_2, beta_3, xi_3, size_X2, size_X3);
				}

			}
		}
	}
	double volume = soap_str->cell[0] *soap_str->cell[1]*soap_str->cell[2];
// Stress terms to populate K_predict matrix
	for (iel = 0; iel < nelem; iel++){
		for (iatm_el_tr = 0; iatm_el_tr < mlff_str->natm_train_elemwise[iel]; iatm_el_tr++){
			// col_idx = cum_natm_elem[iel] + iatm_el_tr;
			col_idx = cum_natm_ele_cols[iel][iatm_el_tr];
			// if (iel > 0)
			// 	col_idx = mlff_str->natm_train_elemwise[iel-1] + iatm_el_tr;
			for (iatm_str = 0; iatm_str < soap_str->natom_elem[iel]; iatm_str++){
				atom_idx = cum_natm_elem[iel] + iatm_str;
				// if (iel > 0)
				// 	atom_idx = soap_str->natom_elem[iel-1] + iatm_str;
				for (istress = 0; istress < 6; istress++){
					row_idx = 3*natom + 1 + istress;

					K_predict[row_idx][col_idx] +=  au2GPa*(1.0/volume)*
					stress_scale[istress]*der_soap_kernel(kernel_typ, soap_str->dX2_dF[atom_idx][istress], soap_str->dX3_dF[atom_idx][istress],
						 soap_str->X2[atom_idx], soap_str->X3[atom_idx],
						mlff_str->X2_traindataset[col_idx], mlff_str->X3_traindataset[col_idx],
						beta_2, beta_3, xi_3, size_X2, size_X3);

				}

			}
		}
	}
	free(cum_natm_elem);
	for (int i = 0; i <nelem; i++){
		free(cum_natm_ele_cols[i]);
	}
	free(cum_natm_ele_cols);

 }

/*
add_newtrain_cols function updates the MLFF_Obj by updating design matrix columns etc. for a new local confiugration

[Input]
1. mlff_str: MLFF_Obj structure
2. X2: 2-body ddescriptor of the new local confiugration
3. X3: 3-body ddescriptor of the new local confiugration
4. elem_typ: Element type of the new local confiugration
[Output]
1. mlff_str: MLFF_Obj structure
*/

 void add_newtrain_cols(double *X2, double *X3, int elem_typ, MLFF_Obj *mlff_str){
	int row_idx, col_idx, atom_idx, istr, iatm_str, neighs, j, istress, nelem = mlff_str->nelem;
	int size_X2 = mlff_str->size_X2, size_X3 = mlff_str->size_X3;
	double beta_2 = mlff_str->beta_2, beta_3 = mlff_str->beta_3, xi_3 = mlff_str->xi_3;
	int kernel_typ = mlff_str->kernel_typ;
	// copy the X2 and X3 into train history data 

	for(j = 0; j < size_X2; j++){
		mlff_str->X2_traindataset[mlff_str->n_cols][j] = X2[j];
	}
	for(j = 0; j < size_X3; j++){
		mlff_str->X3_traindataset[mlff_str->n_cols][j] = X3[j];
	}

	int *cum_natm_elem;
	cum_natm_elem = (int *)malloc(sizeof(int)*nelem);
	cum_natm_elem[0] = 0;
	for (int i = 1; i < nelem; i++){
		cum_natm_elem[i] = (mlff_str->soap_descriptor_strdataset+i)->natom_elem[i-1];
	}

	// Energy term to populate A matrix
	col_idx = mlff_str->n_cols;
	row_idx = 0;
	for (istr = 0; istr < mlff_str->n_str; istr++){
		for (iatm_str=0; iatm_str < (mlff_str->soap_descriptor_strdataset+istr)->natom_elem[elem_typ]; iatm_str++){
			atom_idx = cum_natm_elem[elem_typ] + iatm_str;
			// if (elem_typ > 0)
			// 	atom_idx = (mlff_str->soap_descriptor_strdataset+istr)->natom_elem[elem_typ-1] + iatm_str;
			mlff_str->K_train[row_idx][col_idx] +=  (1.0/(mlff_str->soap_descriptor_strdataset+istr)->natom)*
											soap_kernel(kernel_typ, (mlff_str->soap_descriptor_strdataset+istr)->X2[atom_idx], (mlff_str->soap_descriptor_strdataset+istr)->X3[atom_idx],
											X2, X3, beta_2, beta_3, xi_3, size_X2, size_X3);
		}
		row_idx += 1+3*(mlff_str->soap_descriptor_strdataset+istr)->natom+6;
	}

	// Force term to populate K_train matrix
	row_idx = 0;
	int row_index_F;
	for (istr = 0; istr < mlff_str->n_str; istr++){
		for (iatm_str = 0; iatm_str < (mlff_str->soap_descriptor_strdataset+istr)->natom_elem[elem_typ]; iatm_str++){
			atom_idx = cum_natm_elem[elem_typ] + iatm_str;
			// if (elem_typ > 0)
			// 	atom_idx = (mlff_str->soap_descriptor_strdataset+istr)->natom_elem[elem_typ-1] + iatm_str;
			row_index_F = row_idx + 3*atom_idx+1;
			// x-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
			mlff_str->K_train[row_index_F][col_idx] +=  
				der_soap_kernel(kernel_typ, (mlff_str->soap_descriptor_strdataset+istr)->dX2_dX[atom_idx][0], (mlff_str->soap_descriptor_strdataset+istr)->dX3_dX[atom_idx][0],
				 (mlff_str->soap_descriptor_strdataset+istr)->X2[atom_idx], (mlff_str->soap_descriptor_strdataset+istr)->X3[atom_idx],
				X2, X3, beta_2, beta_3, xi_3, size_X2, size_X3);
			// y-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
			mlff_str->K_train[row_index_F+1][col_idx] +=  
				der_soap_kernel(kernel_typ, (mlff_str->soap_descriptor_strdataset+istr)->dX2_dY[atom_idx][0], (mlff_str->soap_descriptor_strdataset+istr)->dX3_dY[atom_idx][0],
				 (mlff_str->soap_descriptor_strdataset+istr)->X2[atom_idx], (mlff_str->soap_descriptor_strdataset+istr)->X3[atom_idx],
				X2, X3, beta_2, beta_3, xi_3, size_X2, size_X3);
			// z-component (w.r.t itself) "because an atom is not considered it's neighbour hence dealt outside the neighs loop"
			mlff_str->K_train[row_index_F+2][col_idx] +=  
				der_soap_kernel(kernel_typ, (mlff_str->soap_descriptor_strdataset+istr)->dX2_dZ[atom_idx][0], (mlff_str->soap_descriptor_strdataset+istr)->dX3_dZ[atom_idx][0],
				 (mlff_str->soap_descriptor_strdataset+istr)->X2[atom_idx], (mlff_str->soap_descriptor_strdataset+istr)->X3[atom_idx],
				X2, X3, beta_2, beta_3, xi_3, size_X2, size_X3);

			for (neighs = 0; neighs < (mlff_str->soap_descriptor_strdataset+istr)->unique_Nneighbors[atom_idx]; neighs++){
				row_index_F = row_idx + 3*(mlff_str->soap_descriptor_strdataset+istr)->unique_neighborList[atom_idx].array[neighs]+1;
				mlff_str->K_train[row_index_F][col_idx] +=  
				der_soap_kernel(kernel_typ, (mlff_str->soap_descriptor_strdataset+istr)->dX2_dX[atom_idx][1+neighs], (mlff_str->soap_descriptor_strdataset+istr)->dX3_dX[atom_idx][1+neighs],
					 (mlff_str->soap_descriptor_strdataset+istr)->X2[atom_idx], (mlff_str->soap_descriptor_strdataset+istr)->X3[atom_idx],
					X2, X3, beta_2, beta_3, xi_3, size_X2, size_X3);
				// y-component (w.r.t neighs neighbour)
				mlff_str->K_train[row_index_F+1][col_idx] +=  
				der_soap_kernel(kernel_typ, (mlff_str->soap_descriptor_strdataset+istr)->dX2_dY[atom_idx][1+neighs], (mlff_str->soap_descriptor_strdataset+istr)->dX3_dY[atom_idx][1+neighs],
					 (mlff_str->soap_descriptor_strdataset+istr)->X2[atom_idx], (mlff_str->soap_descriptor_strdataset+istr)->X3[atom_idx],
					X2, X3, beta_2, beta_3, xi_3, size_X2, size_X3);
				// z-component (w.r.t neighs neighbour)
				mlff_str->K_train[row_index_F+2][col_idx] +=  
				der_soap_kernel(kernel_typ, (mlff_str->soap_descriptor_strdataset+istr)->dX2_dZ[atom_idx][1+neighs], (mlff_str->soap_descriptor_strdataset+istr)->dX3_dZ[atom_idx][1+neighs],
					 (mlff_str->soap_descriptor_strdataset+istr)->X2[atom_idx], (mlff_str->soap_descriptor_strdataset+istr)->X3[atom_idx],
					X2, X3, beta_2, beta_3, xi_3, size_X2, size_X3);
			}		
		}
		row_idx += 1+3*(mlff_str->soap_descriptor_strdataset+istr)->natom+6;
	}



	// Stress term to populate K_train matrix
	row_idx = 0;
	int row_index_stress;
	for (istr = 0; istr < mlff_str->n_str; istr++){
		double volume = (mlff_str->soap_descriptor_strdataset+istr)->cell[0]*(mlff_str->soap_descriptor_strdataset+istr)->cell[1]*(mlff_str->soap_descriptor_strdataset+istr)->cell[2];
		for (iatm_str = 0; iatm_str < (mlff_str->soap_descriptor_strdataset+istr)->natom_elem[elem_typ]; iatm_str++){
			atom_idx = cum_natm_elem[elem_typ] + iatm_str;
			// if (elem_typ > 0)
			// 	atom_idx = (mlff_str->soap_descriptor_strdataset+istr)->natom_elem[elem_typ-1] + iatm_str;
			for (istress = 0; istress < 6; istress++){
				row_index_stress = row_idx + 3*(mlff_str->soap_descriptor_strdataset+istr)->natom + 1 + istress;
				mlff_str->K_train[row_index_stress][col_idx] +=  au2GPa*(1.0/volume)*
				der_soap_kernel(kernel_typ, (mlff_str->soap_descriptor_strdataset+istr)->dX2_dF[atom_idx][istress], (mlff_str->soap_descriptor_strdataset+istr)->dX3_dF[atom_idx][istress],
					 (mlff_str->soap_descriptor_strdataset+istr)->X2[atom_idx], (mlff_str->soap_descriptor_strdataset+istr)->X3[atom_idx],
					X2, X3, beta_2, beta_3, xi_3, size_X2, size_X3);
			}
		}
		row_idx += 1+3*(mlff_str->soap_descriptor_strdataset+istr)->natom+6;
	}		


	// updating other MLFF parameters such as number of structures, number of training environment, element typ of training env
	mlff_str->natm_typ_train[mlff_str->n_cols] = elem_typ;
	mlff_str->natm_train_total += 1;
	mlff_str->n_cols += 1;
	mlff_str->natm_train_elemwise[elem_typ] += 1;
	free(cum_natm_elem);

}

/*
remove_str_rows function removes a given reference structure from the training dataset

[Input]
1. mlff_str: MLFF_Obj structure
2. str_ID: ID of the reference structure in training dataset
[Output]
1. mlff_str: MLFF_Obj structure
*/

void remove_str_rows(MLFF_Obj *mlff_str, int str_ID){
	int start_idx = 0, end_idx = 0, i, j, istr, rows_to_delete;

	for (i = 0; i < str_ID; i++){
		start_idx += 3*(mlff_str->soap_descriptor_strdataset+i)->natom + 1 + 6;
	}
	end_idx = start_idx + 3*(mlff_str->soap_descriptor_strdataset+str_ID)->natom + 1 + 6;

	rows_to_delete = end_idx - start_idx;
	

	for (istr = str_ID + 1; istr < mlff_str->n_str; istr++){
		copy_descriptors(mlff_str->soap_descriptor_strdataset + istr - 1, mlff_str->soap_descriptor_strdataset + istr );
		// mlff_str->soap_descriptor_strdataset[istr] = mlff_str->soap_descriptor_strdataset[istr + 1];
	}
	delete_soapObj(mlff_str->soap_descriptor_strdataset + mlff_str->n_str);
	mlff_str->n_str = mlff_str->n_str - 1;

	for (i = start_idx; i < end_idx; i++) {
		mlff_str->b_no_norm[i] = mlff_str->b_no_norm[i + rows_to_delete];
		for (j = 0; j < mlff_str->n_cols; j++){
			mlff_str->K_train[i][j] = mlff_str->K_train[i + rows_to_delete][j];
		}
	}

	for (i = mlff_str->n_rows - rows_to_delete; i < mlff_str->n_rows; i++) {
		mlff_str->b_no_norm[i] = 0.0;
		for (j = 0; j < mlff_str->n_cols; j++){
			mlff_str->K_train[i][j] = 0.0;
		}
	}


	mlff_str->n_rows = mlff_str->n_rows - rows_to_delete;
}

/*
remove_train_cols function removes a given local confiugration from the training dataset

[Input]
1. mlff_str: MLFF_Obj structure
2. col_ID: ID of the local confiugration in training dataset
[Output]
1. mlff_str: MLFF_Obj structure
*/
void remove_train_cols(MLFF_Obj *mlff_str, int col_ID){

	int i, j;
	for (i = col_ID; i < mlff_str->n_cols-1; i++){
		for (j = 0; j < mlff_str->n_rows; j++){
			mlff_str->K_train[j][i] = mlff_str->K_train[j][i+1];
		}
	}


	for (j = 0; j < mlff_str->n_rows; j++){
		mlff_str->K_train[j][mlff_str->n_cols] = 0.0;
	}

	for (i =col_ID; i < mlff_str->n_cols-1; i++){
		for (j=0; j < mlff_str->size_X3; j++){
			mlff_str->X3_traindataset[i][j] = mlff_str->X3_traindataset[i+1][j];
		}
		for (j=0; j < mlff_str->size_X2; j++){
			mlff_str->X2_traindataset[i][j] = mlff_str->X2_traindataset[i+1][j];
		}
	}

	for (j=0; j < mlff_str->size_X3; j++){
		mlff_str->X3_traindataset[mlff_str->n_cols][j] = 0.0;
	}
	for (j=0; j < mlff_str->size_X2; j++){
		mlff_str->X2_traindataset[mlff_str->n_cols][j] = 0.0;
	}


	mlff_str->natm_train_elemwise[mlff_str->natm_typ_train[col_ID]] = mlff_str->natm_train_elemwise[mlff_str->natm_typ_train[col_ID]] -1;
	mlff_str->natm_train_total = mlff_str->natm_train_total -1;
	for (i =col_ID; i < mlff_str->n_cols-1; i++){
		mlff_str->natm_typ_train[i] = mlff_str->natm_typ_train[i+1];
	}
	
	mlff_str->n_cols = mlff_str->n_cols - 1;
}

// int main(){
// 	// int N=10, N1 = 20;
// 	// double *X2_str, *X2_tr, *X3_str, *X3_tr, *dX2_str, *dX3_str, kernel_val, der_val;
// 	// double beta_2 = 0.0, beta_3 = 1.0, xi_3 = 4.0;
// 	// time_t t;

// 	// X2_str = (double *) malloc(N*sizeof(double));
// 	// X2_tr = (double *) malloc(N*sizeof(double));
// 	// X3_str = (double *) malloc(N1*sizeof(double));
// 	// X3_tr = (double *) malloc(N1*sizeof(double));
// 	// dX2_str = (double *) malloc(N*sizeof(double));
// 	// dX3_str = (double *) malloc(N1*sizeof(double));
// 	// srand((unsigned) time(&t));

// 	// for (int i = 0; i < N; i++){
// 	// 	X2_str[i] = 0.1+i*0.01;
// 	// 	X2_tr[i] = 0.2+i*0.015;
// 	// 	dX2_str[i] = 0.156+i*0.029;
// 	// }

// 	// for (int i = 0; i < N1; i++){
// 	// 	X3_str[i] = 0.3+i*0.01;
// 	// 	X3_tr[i] = 0.4+i*0.01;
// 	// 	dX3_str[i] = 0.13+i*0.023;
// 	// }

// 	// kernel_val = soap_kernel(X2_str, X3_str, X2_tr, X3_tr, beta_2, beta_3, xi_3, N, N1);
// 	// der_val = der_soap_kernel(dX2_str, dX3_str, X2_str, X3_str, X2_tr, X3_tr, beta_2, beta_3, xi_3, N, N1);

// 	// printf("Kernel values: %10.9f, der_kernel_value: %10.9f\n", kernel_val, der_val);
// 	// free(X2_str);
// 	// free(X2_tr);
// 	// free(X3_str);
// 	// free(X3_tr);
// 	// free(dX3_str);
// 	// free(dX2_str);

// 	clock_t start, end;
// 	start = clock();
// 	int N=12, L=8, N_r = 72;
// 	int i;
// 	double *rgrid, *h_nl, *dh_nl;
// 	NeighList nlist;
// 	SoapObj soap_str;

// 	rgrid = (double *) malloc(sizeof(double)*N_r);
// 	h_nl = (double *) malloc(sizeof(double)*N_r*N*(L+1));
// 	dh_nl = (double *) malloc(sizeof(double)*N_r*N*(L+1));

// 	read_h_nl(N, L, rgrid, h_nl, dh_nl);

// 	double rcut=12;
// 	int nelem=1, natom = 64;
// 	double *atompos;
// 	int *atomtyp;
// 	double BC[3] = {1,1,1}, *cell;
// 	cell = (double *)malloc(sizeof(double)*3);
// 	cell[0] = 19.8;
// 	cell[1] = 19.8;
// 	cell[2] = 19.8;
// 	atompos = (double *)malloc(sizeof(double)*3*natom);
// 	atomtyp = (int *)malloc(sizeof(int)*natom);

// 	for (int i=0; i<natom; i++){
// 		atomtyp[i]=0;
// 	}
	
// 	FILE *fp;
// 	char a1[512];

// 	fp = fopen("atompos.txt","r");
// 	for (int i =0; i < natom; i++){
// 		fscanf(fp, "%s", a1);
// 		atompos[3*i] = atof(a1);
// 		fscanf(fp, "%s", a1);
// 		atompos[3*i+1] = atof(a1);
// 		fscanf(fp, "%s", a1);
// 		atompos[3*i+2] = atof(a1);
// 	}

// 	fclose(fp);

// 	double E=-20, *stress, *force;
// 	stress = (double *)malloc(6*sizeof(double));
// 	force = (double *)malloc(natom*3*sizeof(double));
// 	for (int i=0; i<6; i++){
// 		stress[i] =6+i;
// 	}
// 	for (int i=0; i<natom*3; i++){
// 		force[i] =0.001+0.002*i;
// 	}



// 	MLFF_Obj mlff_str;

// 	initialize_nlist(&nlist, natom,  rcut,  nelem);

// 	build_nlist(rcut, nelem, natom, atompos, atomtyp, BC, cell, &nlist);

// 	double beta_3 = 1.0;
// 	double xi_3 = 4.0;
// 	soap_str.N_rgrid = N_r;

// 	build_soapObj(&soap_str, &nlist,  rgrid, h_nl, dh_nl, atompos, N, L, beta_3, xi_3);

// 	init_MLFF(&soap_str, &mlff_str, 10, 600);
// 	add_firstMD(&soap_str, &nlist, &mlff_str, E, force, stress);
// 	add_newstr_rows(&soap_str, &nlist, &mlff_str, E, force, stress);
// 	double **K_p;
// 	K_p = (double **) malloc(sizeof(double*) * (3*natom+6+1));
// 	for (int i=0; i<(3*natom+6+1); i++){
// 		K_p[i] = (double *) malloc(sizeof(double)*mlff_str.n_cols);
// 	}
// 	calculate_Kpredict(&soap_str, &nlist, &mlff_str, K_p);
// 	remove_str_rows(&mlff_str, 0);
// 	remove_train_cols(&mlff_str, 5);

// 	end = clock();
//     double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

//     printf("time taken: %f\n",cpu_time_used);


// 	free(rgrid);
// 	free(h_nl);
// 	free(dh_nl);

// 	delete_soapObj(&soap_str);
// 	clear_nlist(&nlist);
// 	free_MLFF(&mlff_str);
// 	free(stress);
// 	free(force);
// 	free(K_p);

// 	return 0;

// }
