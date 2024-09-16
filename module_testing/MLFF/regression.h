#ifndef REGRESSION_H
#define REGRESSION_H

#include "mlff_types.h"

void mlff_train_SVD(MLFF_Obj *mlff_str);
void mlff_predict(double *K_predict, MLFF_Obj *mlff_str, double *E,  double* F, double* stress, double* error_bayesian, int natoms );
void mlff_train_Bayesian(MLFF_Obj *mlff_str);
void hyperparameter_Bayesian(double *A, double *AtA, double *Atb, double *b, MLFF_Obj *mlff_str);
double get_regularization_min(double *A, int size);
#endif
