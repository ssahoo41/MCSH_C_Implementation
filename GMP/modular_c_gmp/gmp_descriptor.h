#ifndef CALCULATE_GMPORDERNORM_NODERIV_H
#define CALCULATE_GMPORDERNORM_NODERIV_H

#include "mlff_types.h"

void initialize_nlist(NeighList* nlist, const int natom, const double rcut, const int nelem );
void build_nlist(const double rcut, const int nelem, const int natom, const double * const atompos,
				int * atomtyp, int* BC, double* cell, NeighList* nlist);
void initialize_gmpObj(GMPObj *gmp_str, NeighList *nlist, int *cal_atoms, int cal_num, int **params_i, double **params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order);
int uniqueEle(int* a, int n);
//void delete_gmpObj(GMPObj *gmp_str);
void build_gmpObj(GMPObj *gmp_str, NeighList *nlist, FeatureScaler *ftr_scale, int nmcsh, double *atompos, int **params_i, double **params_d, double** atom_gaussian, int* ngaussians, int* element_index_to_order, int* atom_type_to_indices, int* atom_indices_p);
void scale_features(GMPObj *gmp_str, FeatureScaler *ftr_scale, int nmcsh);
/*
int lin_search(int *arr, int n, int x);
void clear_nlist(NeighList* nlist);

*/
#endif