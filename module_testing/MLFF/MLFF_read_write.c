#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include "mkl.h"
#include "tools_mlff.h"
#include "spherical_harmonics.h"
#include "soap_descriptor.h"
#include "mlff_types.h"
#include "linearsys.h"
#include "sparsification.h"
#include "ddbp_tools.h"
#include "tools.h"
#include "regression.h"
#define au2GPa 29421.02648438959

void intialize_print_MLFF(MLFF_Obj *mlff_str, SPARC_OBJ *pSPARC){

    printf("intialize_print_MLFF called ...\n");
	FILE *fptr;
    char str1[512] = "MLFF_data_reference_atoms.txt"; 
    strcpy(mlff_str->ref_atom_name, str1);
	// mlff_str->ref_atom_name = 'MLFF_data_reference_atoms.txt';
	fptr = fopen(mlff_str->ref_atom_name,"w");

   if(fptr == NULL)
   {
      printf("Error opening the file!");   
      exit(1);             
   }
   fprintf(fptr, "SPARC MLFF on-the-fly trained\n");
   // fprintf(fptr, "nelem:\n")
   // fprintf("%d\n",mlff_str->nelem);
   // fprintf(fptr, "Znucl:\n"):
   // for (int i=0; i < mlff_str->nelem; i++)
   //    fprintf(fptr, "%d\n",pSPARC->Znucl[i]);
   // fprintf(fptr, "beta_2:\n");
   // fprintf(fptr, "%f\n", mlff_str->beta_2);
   // fprintf(fptr, "beta_3:\n");
   // fprintf(fptr, "%f\n", mlff_str->beta_3);
   // fprintf(fptr, "xi_3:\n");
   // fprintf(fptr, "%f\n", mlff_str->xi_3);
   // fprintf(fptr, "N_max:\n");
   // fprintf(fptr, "%f\n", mlff_str->NMAX);
   // fprintf(fptr, "L_max:\n");
   // fprintf(fptr, "%f\n", mlff_str->LMAX);
   // fprintf(fptr, "Rcut:\n");
   // fprintf(fptr, "%f\n", mlff_str->rcut);
   // fprintf(fptr, "size_X2:\n");
   // fprintf(fptr, "%d\n",mlff_str->size_X2);
   // fprintf(fptr, "size_X3:\n");
   // fprintf(fptr, "%d\n",mlff_str->size_X3);

   fprintf(fptr, "SOAP descriptors:\n");
   fclose(fptr);
   char str2[512] = "MLFF_data_reference_structures.txt"; 
    strcpy(mlff_str->ref_str_name, str2);
   // mlff_str->ref_str_name = "MLFF_data_reference_structures.txt";
   fptr = fopen(mlff_str->ref_str_name,"w");
   if(fptr == NULL)
   {
      printf("Error opening the file!");   
      exit(1);             
   }
   fprintf(fptr, "SPARC MLFF on-the-fly trained\n");
   // fprintf(fptr, "nelem:\n")
   // fprintf("%d\n",mlff_str->nelem);
   // fprintf(fptr, "Znucl:\n"):
   // for (int i=0; i < mlff_str->nelem; i++)
   //    fprintf(fptr, "%d\n",pSPARC->Znucl[i]);
   // fprintf(fptr, "beta_2:\n");
   // fprintf(fptr, "%f\n", mlff_str->beta_2);
   // fprintf(fptr, "beta_3:\n");
   // fprintf(fptr, "%f\n", mlff_str->beta_3);
   // fprintf(fptr, "xi_3:\n");
   // fprintf(fptr, "%f\n", mlff_str->xi_3);
   // fprintf(fptr, "N_max:\n");
   // fprintf(fptr, "%f\n", mlff_str->NMAX);
   // fprintf(fptr, "L_max:\n");
   // fprintf(fptr, "%f\n", mlff_str->LMAX);
   // fprintf(fptr, "Rcut:\n");
   // fprintf(fptr, "%f\n", mlff_str->rcut);
   // fprintf(fptr, "size_X2:\n");
   // fprintf(fptr, "%d\n",mlff_str->size_X2);
   // fprintf(fptr, "size_X3:\n");
   // fprintf(fptr, "%d\n",mlff_str->size_X3);
   fprintf(fptr, "Structures:\n");
   fclose(fptr);

   printf("intialize_print_MLFF returned ...\n");
}

void print_new_ref_atom_MLFF(MLFF_Obj *mlff_str, int elem_typ, int nimg, double *X2, double *X3){
	FILE *fptr;
    // char str1[512] = "MLFF_data_reference_atoms.txt"; 
    // strcpy(mlff_str->ref_atom_name, str1);
   // mlff_str->ref_atom_name = "MLFF_data_reference_atoms.txt";
	fptr = fopen(mlff_str->ref_atom_name,"a");

   if(fptr == NULL)
   {
      printf("Error opening the file!");   
      exit(1);             
   }
   fprintf(fptr, "elem_typ: %d, reference_atom_no: %d\n",elem_typ, nimg);
   for (int i=0; i < mlff_str->size_X2; i++)
   	fprintf(fptr, "%10.9f ",X2[i]);
   fprintf(fptr,"\n");
   for (int i=0; i < mlff_str->size_X3; i++)
   	fprintf(fptr, "%10.9f ",X3[i]);
   fprintf(fptr,"\n");
   fclose(fptr);

}

void print_new_ref_structure_MLFF(MLFF_Obj *mlff_str, int nstr, SoapObj *soap_str, double *atompos, double Etot, double *force, double *stress){
    FILE *fptr;
	fptr = fopen(mlff_str->ref_str_name,"a");

   if(fptr == NULL)
   {
      printf("Error opening the file!");   
      exit(1);             
   }
   fprintf(fptr, "structure_no: %d\n", nstr);
   fprintf(fptr, "CELL: \n");
   fprintf(fptr, "%f %f %f\n", soap_str->cell[0], soap_str->cell[1], soap_str->cell[2]);
   fprintf(fptr,"natom: \n");
   fprintf(fptr, "%d\n", soap_str->natom);
   fprintf(fptr,"natom_elem:\n");
   for (int i=0; i <soap_str->nelem; i++)
   	fprintf(fptr,"%d\n",soap_str->natom_elem[i]);
   fprintf(fptr,"Atom_positions: \n");
   for (int i=0; i < soap_str->natom; i++)
   	fprintf(fptr, "%10.9f %10.9f %10.9f\n",atompos[3*i],atompos[3*i+1],atompos[3*i+2]);

   fprintf(fptr, "Etot(Ha):\n");
   fprintf(fptr, "%10.9f\n",Etot);
   fprintf(fptr, "F(Ha/bohr):\n");
   for (int i=0; i <soap_str->natom; i++)
   	fprintf(fptr, "%10.9f %10.9f %10.9f\n", force[3*i],force[3*i+1],force[3*i+2]);
   fprintf(fptr, "Stress(GPa)\n");
   fprintf(fptr,"%10.9f %10.9f %10.9f\n",au2GPa*stress[0],au2GPa*stress[3], au2GPa*stress[5]);
   fprintf(fptr,"%10.9f %10.9f %10.9f\n",au2GPa*stress[3],au2GPa*stress[1], au2GPa*stress[4]);
   fprintf(fptr,"%10.9f %10.9f %10.9f\n",au2GPa*stress[5],au2GPa*stress[4], au2GPa*stress[2]);
   fclose(fptr);

}

void print_restart_MLFF(MLFF_Obj *mlff_str, SPARC_OBJ *pSPARC, int nref_str, int *nref_atoms){

	FILE *fptr;
    char str1[512] = "MLFF_RESTART.txt"; 
    strcpy(mlff_str->restart_name, str1);
   // mlff_str->restart_name = "MLFF_RESTART.txt";
	fptr = fopen(mlff_str->restart_name,"w");

   if(fptr == NULL)
   {
      printf("Error opening the file!");   
      exit(1);             
   }
   fprintf(fptr, "SPARC MLFF on-the-fly trained\n");
   fprintf(fptr, "nelem:\n");
   fprintf(fptr,"%d\n",mlff_str->nelem);
   fprintf(fptr, "Znucl:\n");
   for (int i=0; i < mlff_str->nelem; i++)
   	fprintf(fptr, "%d\n",pSPARC->Znucl[i]);
   fprintf(fptr, "beta_2:\n");
   fprintf(fptr, "%f\n", mlff_str->beta_2);
   fprintf(fptr, "beta_3:\n");
   fprintf(fptr, "%f\n", mlff_str->beta_3);
   fprintf(fptr, "xi_3:\n");
   fprintf(fptr, "%f\n", mlff_str->xi_3);
   fprintf(fptr, "N_max:\n");
   fprintf(fptr, "%d\n", mlff_str->Nmax);
   fprintf(fptr, "L_max:\n");
   fprintf(fptr, "%d\n", mlff_str->Lmax);
   fprintf(fptr, "Rcut:\n");
   fprintf(fptr, "%f\n", mlff_str->rcut);
   fprintf(fptr, "size_X2:\n");
   fprintf(fptr, "%d\n",mlff_str->size_X2);
   fprintf(fptr, "size_X3:\n");
   fprintf(fptr, "%d\n",mlff_str->size_X3);
   fprintf(fptr, "mu_E:\n");
   fprintf(fptr,"%10.9f\n",mlff_str->mu_E);
   fprintf(fptr, "std_E:\n");
   fprintf(fptr,"%10.9f\n",mlff_str->std_E);
   fprintf(fptr, "std_F:\n");
   fprintf(fptr, "%10.9f\n",mlff_str->std_F);
   fprintf(fptr, "mu_Stress:\n");
   fprintf(fptr, "%10.9f\n %10.9f\n %10.9f\n %10.9f\n %10.9f\n %10.9f\n",mlff_str->mu_stress[0],mlff_str->mu_stress[1],mlff_str->mu_stress[2],
            mlff_str->mu_stress[3],mlff_str->mu_stress[4],mlff_str->mu_stress[5]);
   fprintf(fptr, "std_Stress:\n");
   fprintf(fptr, "%10.9f\n %10.9f\n %10.9f\n %10.9f\n %10.9f\n %10.9f\n",mlff_str->std_stress[0],mlff_str->std_stress[1],mlff_str->std_stress[2],
            mlff_str->std_stress[3],mlff_str->std_stress[4],mlff_str->std_stress[5]);

   fprintf(fptr, "E_scale:\n");
   fprintf(fptr, "%f\n",mlff_str->E_scale);
   fprintf(fptr, "F_scale:\n");
   fprintf(fptr, "%f\n",mlff_str->F_scale);
   fprintf(fptr, "stress_scale:\n");
   fprintf(fptr, "%f %f %f %f %f %f\n",mlff_str->stress_scale[0],mlff_str->stress_scale[1],mlff_str->stress_scale[2],
                  mlff_str->stress_scale[3],mlff_str->stress_scale[4],mlff_str->stress_scale[5]);

   fprintf(fptr, "N_ref_str: \n");
   fprintf(fptr, "%d\n", mlff_str->n_str);

   fprintf(fptr, "N_ref_atoms_elemwise: \n");
   for (int i=0; i < mlff_str->nelem; i++)
   	fprintf(fptr, "%d\n", mlff_str->natm_train_elemwise[i]);

   fprintf(fptr, "n_rows: \n");
   fprintf(fptr, "%d\n", mlff_str->n_rows);

   fprintf(fptr, "n_cols: \n");
   fprintf(fptr, "%d\n", mlff_str->n_cols);

   fprintf(fptr, "K_train:\n");
   for (int i = 0; i < mlff_str->n_rows; i++){
      for (int j=0; j< mlff_str->n_cols; j++){
         fprintf(fptr, "%10.9f ",mlff_str->K_train[i][j]);
      }
      fprintf(fptr, "\n");
   }

   fprintf(fptr, "b:\n");
   for (int i = 0; i < mlff_str->n_rows; i++){
      fprintf(fptr, "%10.9f\n",mlff_str->b_no_norm[i]);
   }

   fprintf(fptr, "weights_regression:\n");
   for (int i = 0; i < mlff_str->n_cols; i++){
      fprintf(fptr, "%10.9f\n",mlff_str->weights[i]);
   }

   fclose(fptr);
}

void read_MLFF_files(MLFF_Obj *mlff_str, SPARC_OBJ *pSPARC){
   FILE *fptr;
   int info;
   char a1[512], str[512];

   // mlff_str->restart_name = "MLFF_RESTART.txt";
   fptr = fopen(mlff_str->restart_name,"r");

   if(fptr == NULL)
   {
      printf("Error opening the file!");   
      exit(1);             
   }

   fgets(a1, sizeof (a1), fptr);
   // fgets(a1, sizeof (a1), fp);//Number of element types:

   // info = fscanf(fp, "%s", a1);
   // int nelem = atoi(a1);
   // if (nelem ~= pSPARC->Ntypes){
   //    printf("Number of element types in MLFF file and the input file for DFT are not the same!");
   //    exit(1); 
   // }
   // fgets(a1, sizeof (a1), fp); // Elements:
   int nelem;
   int Znucl[pSPARC->Ntypes];
   double K_elem;
   while (!feof(fptr)){
      fscanf(fptr,"%s",str);
      fscanf(fptr, "%*[^\n]\n");
      // enable commenting with '#'
     if (str[0] == '#' || str[0] == '\n'|| strcmpi(str,"undefined") == 0) {
         fscanf(fptr, "%*[^\n]\n"); // skip current line
         continue;
     }
     if (strcmpi(str,"nelem:") == 0) {
         fscanf(fptr,"%d", &nelem);
         if (nelem != pSPARC->Ntypes){
            printf("Number of element types in MLFF file and the input file for DFT are not the same!");
            exit(1); 
         }
         fscanf(fptr, "%*[^\n]\n");

     } else if (strcmpi(str,"Znucl:") == 0){
         for (int i=0; i < pSPARC->Ntypes; i++){
            fscanf(fptr,"%d", &Znucl[i]);
            if (Znucl[i] != pSPARC->Znucl[i]){
               printf("Atomic number of the element types in MLFF file and the input file for DFT are not the same!");
               exit(1);
            }
            fscanf(fptr, "%*[^\n]\n");
         }
     } else if (strcmpi(str,"beta_2:") == 0){
         fscanf(fptr,"%lf", &mlff_str->beta_2);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"beta_3:") == 0){
         fscanf(fptr,"%lf", &mlff_str->beta_3);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"xi_3:") == 0){
         fscanf(fptr,"%lf", &mlff_str->xi_3);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"N_max:") == 0){
         fscanf(fptr,"%d", &mlff_str->Nmax);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"L_max:") == 0){
         fscanf(fptr,"%d", &mlff_str->Lmax);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"Rcut:") == 0){
         fscanf(fptr,"%lf", &mlff_str->rcut);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"size_X2:") == 0){
         fscanf(fptr,"%d", &mlff_str->size_X2);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"size_X3:") == 0){
         fscanf(fptr,"%d", &mlff_str->size_X3);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"mu_E:") == 0){
         fscanf(fptr,"%lf", &mlff_str->mu_E);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"std_E:") == 0){
         fscanf(fptr,"%lf", &mlff_str->std_E);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"std_F:") == 0){
         fscanf(fptr,"%lf", &mlff_str->std_F);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"mu_Stress:") == 0){
         for (int i=0; i < 6; i++){
            fscanf(fptr,"%lf", &mlff_str->mu_stress[i]);
            fscanf(fptr, "%*[^\n]\n");
         }  
     } else if (strcmpi(str,"std_Stress:") == 0){
         for (int i=0; i < 6; i++){
            fscanf(fptr,"%lf", &mlff_str->std_stress[i]);
            fscanf(fptr, "%*[^\n]\n");
         }   
     } else if (strcmpi(str,"E_scale:") == 0){
         fscanf(fptr,"%lf", &mlff_str->E_scale);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"F_scale:") == 0){
         fscanf(fptr,"%lf", &mlff_str->F_scale);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"stress_scale:") == 0){
         for (int i=0; i <6; i++){
            fscanf(fptr,"%lf", &mlff_str->stress_scale[i]);
            fscanf(fptr, "%*[^\n]\n");
         }
     } else if (strcmpi(str,"N_ref_str:") == 0){
         fscanf(fptr,"%d", &mlff_str->n_str);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"N_ref_atoms_elemwise:") == 0){
      int natm_train_total=0;
         for (int i = 0; i < pSPARC->Ntypes; i++){
            fscanf(fptr,"%d", &mlff_str->natm_train_elemwise[i]);
            natm_train_total += mlff_str->natm_train_elemwise[i];
            fscanf(fptr, "%*[^\n]\n");
         }
         mlff_str->natm_train_total = natm_train_total;    
     } else if (strcmpi(str,"n_rows:") == 0){
         fscanf(fptr,"%d", &mlff_str->n_rows);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"n_cols:") == 0){
         fscanf(fptr,"%d", &mlff_str->n_cols);
         fscanf(fptr, "%*[^\n]\n");
     } else if (strcmpi(str,"K_train:") == 0){
         for (int i = 0; i < mlff_str->n_rows; i++){
            for (int j = 0; j < mlff_str->n_cols; j++){
               fscanf(fptr,"%lf", &K_elem);
               mlff_str->K_train[i][j] = K_elem;
            }
            fscanf(fptr, "%*[^\n]\n");
         }       
     } else if (strcmpi(str,"b:") == 0){
         for (int i = 0; i < mlff_str->n_rows; i++){
            fscanf(fptr,"%lf", &K_elem);
            mlff_str->b_no_norm[i] = K_elem;
            fscanf(fptr, "%*[^\n]\n");
         }
     } else if (strcmpi(str,"weights_regression:") == 0){
         for (int i = 0; i < mlff_str->n_cols; i++){
            fscanf(fptr,"%lf", &K_elem);
            mlff_str->weights[i] = K_elem;
            fscanf(fptr, "%*[^\n]\n");
         }
     } else {
            printf("\nCannot recognize input variable identifier: \"%s\"\n",str);
            exit(EXIT_FAILURE);
     }
   }

   fclose(fptr);

   fptr = fopen(mlff_str->ref_atom_name,"r");
   if(fptr == NULL)
   {
      printf("Error opening the file!");   
      exit(1);             
   }
   fgets(a1, sizeof (a1), fptr);

   int atm_typ, nimg;
   int count=0;
   while (!feof(fptr)){
      fscanf(fptr,"%s",str);
      fscanf(fptr, "%*[^\n]\n");
      // enable commenting with '#'
     if (str[0] == '#' || str[0] == '\n'|| strcmpi(str,"undefined") == 0) {
         fscanf(fptr, "%*[^\n]\n"); // skip current line
         continue;
     }

     if (strcmpi(str,"elem_typ:") == 0) {
         fscanf(fptr,"%d", &atm_typ);
         fscanf(fptr,"%s", a1);
         fscanf(fptr,"%d", &nimg);
         fscanf(fptr, "%*[^\n]\n");
         mlff_str->natm_typ_train[count] = atm_typ;
         for (int i=0; i < mlff_str->size_X2; i++){
            fscanf(fptr,"%lf", &mlff_str->X2_traindataset[count][i]);
         }
         fscanf(fptr, "%*[^\n]\n");
         for (int i=0; i < mlff_str->size_X3; i++){
            fscanf(fptr,"%lf", &mlff_str->X3_traindataset[count][i]);
         }
         count++;
     } 
   }
   fclose(fptr);


   fptr = fopen(mlff_str->ref_str_name,"r");
   if(fptr == NULL)
   {
      printf("Error opening the file!");   
      exit(1);             
   }
   fgets(a1, sizeof (a1), fptr);

   int str_no, natom, natom_elem[pSPARC->Ntypes], BC[3] ={1,1,1};
   double cell[3], atompos[3*pSPARC->n_atom];
   NeighList *nlist;
   // SoapObj *soap_str;
   nlist = (NeighList *) malloc(sizeof(NeighList)*1);
   // soap_str = (SoapObj *) malloc(sizeof(SoapObj)*1);
   int *atomtyp;
   int count0=0;
   while (!feof(fptr)){
      fscanf(fptr,"%s",str);
      fscanf(fptr, "%*[^\n]\n");
      // enable commenting with '#'
     if (str[0] == '#' || str[0] == '\n'|| strcmpi(str,"undefined") == 0) {
         fscanf(fptr, "%*[^\n]\n"); // skip current line
         continue;
     }

     if (strcmpi(str,"structure_no:") == 0) {
         fscanf(fptr,"%d", &str_no);
         fscanf(fptr, "%*[^\n]\n");
         fscanf(fptr,"%s",str);
         fscanf(fptr,"%lf", &cell[0]);
         fscanf(fptr,"%lf", &cell[1]);
         fscanf(fptr,"%lf", &cell[2]);
         fscanf(fptr, "%*[^\n]\n");
         fscanf(fptr,"%s",str);
         fscanf(fptr,"%d", &natom);
         if (natom != pSPARC->n_atom){
            printf("natoms not matching between MLFF file and DFT ion file!\n");
            exit(1);
         }
         fscanf(fptr,"%s",str);
         for (int i=0; i < pSPARC->Ntypes; i++){
            fscanf(fptr,"%d", &natom_elem[i]);
         }
         fscanf(fptr, "%*[^\n]\n");
         fscanf(fptr,"%s",str);
         for (int i=0; i < pSPARC->n_atom; i++){
            fscanf(fptr,"%lf", &atompos[3*i]);
            fscanf(fptr,"%lf", &atompos[3*i+1]);
            fscanf(fptr,"%lf", &atompos[3*i+2]);
            fscanf(fptr, "%*[^\n]\n");
         }

         atomtyp = (int *) malloc(sizeof(int)*natom);
         int count1=0;
         for (int i = 0; i < pSPARC->Ntypes; i++){
            for (int j = 0; j < natom_elem[i]; j++){
               atomtyp[count] = i;
               count++;
            }
         }
         initialize_nlist(nlist, natom, mlff_str->rcut, pSPARC->Ntypes);
         build_nlist(mlff_str->rcut, pSPARC->Ntypes, pSPARC->n_atom, atompos, atomtyp, BC, cell, nlist);
         initialize_soapObj(&mlff_str->soap_descriptor_strdataset[count0], nlist, mlff_str->Lmax, mlff_str->Nmax, mlff_str->beta_3, mlff_str->xi_3);
         build_soapObj(&mlff_str->soap_descriptor_strdataset[count0], nlist, mlff_str->rgrid, mlff_str->h_nl, mlff_str->dh_nl, pSPARC->atom_pos, mlff_str->Nmax,
                mlff_str->Lmax, mlff_str->beta_3, mlff_str->xi_3);
         clear_nlist(nlist);
         free(atomtyp);
         count0++;
     } 

   }

   free(nlist);

}

