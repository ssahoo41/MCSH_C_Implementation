#ifndef MLFF_READ_WRITE_H
#define MLFF_READ_WRITE_H

#include "mlff_types.h"
#include "isddft.h"

void intialize_print_MLFF(MLFF_Obj *mlff_str, SPARC_OBJ *pSPARC);
void print_new_ref_atom_MLFF(MLFF_Obj *mlff_str, int elem_typ, int nimg, double *X2, double *X3);
void print_new_ref_structure_MLFF(MLFF_Obj *mlff_str, int nstr, SoapObj *soap_str, double *atompos, double Etot, double *force, double *stress);
void print_restart_MLFF(MLFF_Obj *mlff_str, SPARC_OBJ *pSPARC, int nref_str, int *nref_atoms);

#endif
