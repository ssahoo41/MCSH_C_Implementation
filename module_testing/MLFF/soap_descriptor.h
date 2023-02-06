#ifndef SOAP_DESCRIPTOR_H
#define SOAP_DESCRIPTOR_H

#include "mlff_types.h"

void read_h_nl(const int N, const int L, double *rgrid, double *h_nl, double *dh_nl);
void initialize_nlist(NeighList* nlist, const int natom, const double rcut, const int nelem );
int lin_search(int *arr, int n, int x);
void build_nlist(const double rcut, const int nelem, const int natom, const double * const atompos,
				 int * atomtyp, int* BC, double* cell, NeighList* nlist);
void clear_nlist(NeighList* nlist);
int uniqueEle(int* a, int n);
void initialize_soapObj(SoapObj *soap_str, NeighList *nlist, int Lmax, int Nmax, double beta_3, double xi_3);
void delete_soapObj(SoapObj *soap_str);
void build_soapObj(SoapObj *soap_str, NeighList *nlist, double* rgrid, double* h_nl, double* dh_nl, double *atompos, int Nmax, int Lmax, double beta_3, double xi_3);


void print_descriptors(double *X2, double *X3, double **dX2_dX, double **dX2_dY, double **dX2_dZ, double **dX2_dF, double **dX3_dF, double **dX3_dX, double **dX3_dY, double **dX3_dZ, int size_X2, int size_X3, int neighs);
void print_cnlm(double complex **cnlm, double complex ***dcnlm_dX, double complex ***dcnlm_dY, double complex ***dcnlm_dZ, double complex ***dcnlm_dF, int size_cnlm, int* neighs, int nelem);

void build_soapObj_wZ(SoapObj *soap_str, NeighList *nlist, double* rgrid, double* h_nl, double* dh_nl, double *atompos, int Nmax, int Lmax, double beta_3, double xi_3, int* Z);
void delete_soapObj_wZ(SoapObj *soap_str);
void initialize_soapObj_wZ(SoapObj *soap_str, NeighList *nlist, int Lmax, int Nmax, double beta_3, double xi_3);
#endif
