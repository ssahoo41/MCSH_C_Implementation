#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <ctype.h>
#include "mkl.h"
#include "tools_mlff.h"
#include "mlff_types.h"
#include "local_type_tools.h"
#include "mlff_types.h"
#include "gmp_deriv_descriptor.h"

#define L_STRING 512
#define L_PSD 4096
#define L_FILENAME 4096
#define MAX_PATH_LEN 10240
#define L_ATMTYPE 8
#define RED   "\x1B[31m"
#define RESET "\x1B[0m"
//Consider how to handle Bohr vs Angstrom! !IMPORTANT!
int check_num_input(FILE *fp, void *array, char TYPE);
void check_below_entries(FILE *ion_fp, char *tag);
void check_below_entries(FILE *ion_fp, char *tag);
void read_ion(READ_Obj *read_obj);
void read_input(READ_Obj *read_obj);
void read_GMP_inpt(MLFF_Obj *mlff_str, READ_Obj *read_obj);
int strcmpi(const char *p1, const char *p2);
// Write main function here
int main(int argc, char **argv){
    MLFF_Obj *mlff_str = (MLFF_Obj*) malloc(sizeof(MLFF_Obj));
    READ_Obj *read_obj = (READ_Obj*) malloc(sizeof(READ_Obj));
    GMPObj *gmp_str = (GMPObj*) malloc(sizeof(GMPObj));
    NeighList *nlist = (NeighList*) malloc(sizeof(NeighList));
    FeatureScaler *ftr_scale = (FeatureScaler*) malloc(sizeof(FeatureScaler));
    read_obj->filename = argv[1];
    printf("Memory allocated for structures and filename read\n");
    printf("Filename: %s\n",read_obj->filename);

    read_input(read_obj);
    printf("Read input\n");
    read_ion(read_obj);
    printf("Read ion\n");
    read_GMP_inpt(mlff_str,read_obj);
    printf("Read GMP input\n");
    int nelem = read_obj->Ntypes;
	int natom = read_obj->n_atom;

    initialize_nlist(nlist, read_obj->n_atom, mlff_str->rcut, read_obj->Ntypes);
	printf("Initialization of nlist done.\n");

    int *BC = (int*)malloc(3*sizeof(int));
	double *cell =(double *) malloc(3*sizeof(double));
	int *atomtyp = (int *) malloc(read_obj->n_atom*sizeof(int));

	BC[0] = read_obj->BCx; BC[1] = read_obj->BCy; BC[2] = read_obj->BCz;
	cell[0] = read_obj->range_x;
	cell[1] = read_obj->range_y;
	cell[2] = read_obj->range_z;

    int count = 0;
	for (int i=0; i < read_obj->Ntypes; i++){
		for (int j=0; j < read_obj->nAtomv[i]; j++){
			atomtyp[count] = i;
			count++;
		}
	}

	build_nlist(mlff_str->rcut, read_obj->Ntypes, read_obj->n_atom, read_obj->atom_pos, atomtyp, BC, cell, nlist);
	printf("Building nlist done. \n");
    initialize_gmpObj(gmp_str, nlist, mlff_str->cal_atoms, mlff_str->cal_num, mlff_str->params_i, 
		mlff_str->params_d, mlff_str->nmcsh, mlff_str->atom_gaussians_p, 
		mlff_str->ngaussians, mlff_str->element_index_to_order_p);
    printf("GMPObj initialized\n");
    
    build_gmpObj(gmp_str, nlist, ftr_scale, mlff_str->nmcsh, read_obj->atom_pos, 
		mlff_str->params_i, mlff_str->params_d, mlff_str->atom_gaussians_p, 
		mlff_str->ngaussians, mlff_str->element_index_to_order_p, mlff_str->atom_type_to_indices_p, 
		mlff_str->atom_indices_p);
    printf("Building gmp donexx. \n");

    FILE *fptr;
    fptr = fopen("num_out.txt","w");
    for (int i=0; i < 108; i++){
        for (int j=0; j < 15; j++){
            fprintf(fptr,"%d,%d,%f",i,j,gmp_str->X[i][j]);
            for (int k=0; k<108; k++){
                fprintf(fptr,",%f,%f,%f",gmp_str->dX_dX[i][k][j],gmp_str->dX_dY[i][k][j],gmp_str->dX_dZ[i][k][j]);
            }
            fprintf(fptr,"\n");
        }
    }
    fclose(fptr);
    return 0;
}
//End of main function

void read_GMP_inpt(MLFF_Obj *mlff_str, READ_Obj *read_obj){
    FILE *input_fp = fopen("GMP_input.in", "r");
	if (NULL == input_fp) printf("GMP input file cannot be opened\n");

	char *str            = malloc(L_STRING * sizeof(char));
    char *str_tmp = malloc(L_PSD * sizeof(char));
    char *pseudos;
    double *params_di = (double*) malloc(sizeof(double)*4);  // why hardcode 4? Ask Lucas
    int    *params_ii = (int*) malloc(sizeof(int)*2); ;  // why hardcode 2? Ask Lucas
    int *elem, *ngaussians;
    int nsig, nang, Ntypes, maxprim;
    Ntypes = read_obj->Ntypes;

	while (!feof(input_fp)){
        int i = 0;
        int count = fscanf(input_fp,"%s",str);
        if (count < 0) continue;  // for some specific cases
        // enable commenting with '#'
        if (str[0] == '#' || str[0] == '\n'|| strcmpi(str,"undefined") == 0) {
        	
            fscanf(input_fp, "%*[^\n]\n"); // skip current line
            continue;
        }
        // Read number of radial probes - integer
        if (strcmpi(str,"N_SIGMA:") == 0){
        	
            fscanf(input_fp,"%d",&nsig);
            fscanf(input_fp, "%*[^\n]\n");
        // Read number of angular probes - integer
        } else if (strcmpi(str,"N_ANGULAR:") == 0){
        	
            fscanf(input_fp,"%d",&nang);
            fscanf(input_fp, "%*[^\n]\n");
        // Read sigmas - array of doubles of length nsig
        } else if (strcmpi(str,"SIGMA:") == 0) {
        	
            for (int i=0; i<nsig;i++){
                fscanf(input_fp,"%lf",&params_di[i]);
            }
            fscanf(input_fp, "%*[^\n]\n");
        // Read cutoff distance - double
        } else if (strcmpi(str,"CUTOFF:") == 0) {
        	
            fscanf(input_fp,"%lf",&params_di[3]);
            fscanf(input_fp, "%*[^\n]\n");
        // Read squared (1) or square root (0) features - int
        } else if (strcmpi(str,"SQUARED:") == 0){
        	
            fscanf(input_fp,"%d",&params_ii[0]);
            fscanf(input_fp, "%*[^\n]\n");
        // Read solid (1) or surface (0) harmonics - int
        } else if (strcmpi(str,"SOLID:") == 0){
        	
            fscanf(input_fp,"%d",&params_ii[1]);
            fscanf(input_fp, "%*[^\n]\n");
        // Read array of elements - array of ints corresponding to Ntypes
        } else if (strcmpi(str,"ELEMENTS:") == 0){
        	
            elem = (int*) malloc(sizeof(int)*Ntypes);
            for (int i=0; i<Ntypes; i++){
                fscanf(input_fp,"%d",&elem[i]);
            }
            fscanf(input_fp, "%*[^\n]\n");
        } else if (strcmpi(str,"PSEUDODENSITY:") == 0){
        	pseudos = (char*) malloc(Ntypes*L_PSD*sizeof(char));
            #define STR_(X) #X
            #define STR(X) STR_(X)

            for (int i=0; i<Ntypes; i++){
                fscanf(input_fp,"%" STR(MAX_PATH_LEN) "s",str_tmp); // read at most MAX_PATH_LEN chars

                if (strlen(str_tmp) > L_PSD) {
                    printf("\n[FATAL] PSEUDO_DENS: path length (%ld) exceeds maximum length (%d)\n", strlen(str_tmp), L_PSD);
                    printf("Please provide a shorter path to your pseudopotential\n");
                    exit(EXIT_FAILURE);
                }
                snprintf(&pseudos[i*L_PSD], L_PSD, "%s", str_tmp);

            }
            // free(str_tmp);
            #undef STR
            #undef STR_
            #undef MAX_PATH_LEN
            fscanf(input_fp, "%*[^\n]\n");


        // Read list of number of primitives for each atom - array of ints
        // Current implementation requries relatively strict ordering of inputs in .inpt file
        } else if (strcmpi(str,"N_PRIMITIVES:") == 0) {
        	
            ngaussians = (int*) malloc(sizeof(int)*Ntypes);
            for (int i=0; i<Ntypes; i++){
                fscanf(input_fp,"%d",&ngaussians[i]);
            }
            fscanf(input_fp, "%*[^\n]\n");
        // Read max number of prmitives (could use built in function, but just an input for now) - int
        } else if (strcmpi(str,"MAX_PRIMITIVES:") == 0) {
        	
            fscanf(input_fp,"%d",&maxprim);
            fscanf(input_fp, "%*[^\n]\n");
        }
    }
    
    double *sigmas = (double*) malloc(nsig*sizeof(double));
    for (int i=0; i<nsig; i++){
        sigmas[i] = params_di[i];
    }
    // Store cutoff
    double rcut = params_di[nsig];
    // Get descriptor vector size from nsig and nang probes
    int nmcsh = nsig*nang;
    // Allocate memory for GMP double params and populate 2D array
    fclose(input_fp);
    
    double **params_d = (double**) malloc(sizeof(double*)*nmcsh);
    for (int i = 0; i < nmcsh; i++) params_d[i] = (double*) malloc(sizeof(double)*6); // 6 is hardcoded here
    for (int i = 0; i < nang; i++){
        for (int j = 0; j < nsig; j++){
            params_d[i*3+j][0] = sigmas[j];
            params_d[i*3+j][1] = 1.0;
            params_d[i*3+j][2] = pow(1/(sigmas[j]*sqrt(2.0*M_PI)),3);
            params_d[i*3+j][3] = 1/(2*pow(sigmas[j],2));
            params_d[i*3+j][4] = rcut;
            params_d[i*3+j][5] = 1.0;
        }
    }
 	
    // Allocate memory for GMP integer params and populate 2D array
	int **params_i = (int**) malloc(sizeof(int*)*nmcsh);
    for (int i=0; i<nmcsh; i++) params_i[i] = (int*) malloc(sizeof(int)*3);  // 3 is hardcoded here
	for (int i = 0; i < nang; i++){
        for (int j = 0; j < nsig; j++){
            params_i[i*3+j][0] = i;
            params_i[i*3+j][1] = params_ii[0];
            params_i[i*3+j][2] = params_ii[1];
            
        }
    }	

    // Allocate memory for gaussian primitive parameters and populate array. 2D array determined by Ntypes and max number of primitives
    // Primitives read from file which at current must be in working directory
    double **atom_gaussians_p = (double**) malloc(sizeof(double*)*Ntypes);
    for (int i=0; i<Ntypes; i++) atom_gaussians_p[i] = (double*) malloc(sizeof(double)*maxprim*2);
    // FILE *input_fp;
    double temp_read, temp_read1;
    for (int i=0; i<Ntypes; i++){
        // FILE *input_fp = fopen("Al_pseudodensity_4.g", "r");
        input_fp = fopen(pseudos + i*L_PSD, "r");
        printf("%s\n",pseudos + i*L_PSD);
        if (NULL == input_fp) printf("File cannot Al_ps be opened\n");
        for (int j=0; j<ngaussians[i]; j++){
            // fscanf(input_fp,"%f %f",&atom_gaussians_p[i][2*j],&atom_gaussians_p[i][2*j+1]);
            fscanf(input_fp, "%lf\t%lf", &temp_read, &temp_read1);
			atom_gaussians_p[i][2*j] = temp_read;
			atom_gaussians_p[i][2*j+1] = temp_read1;
            printf("%f %f\n",atom_gaussians_p[i][2*j],atom_gaussians_p[i][2*j+1]);
            fscanf(input_fp, "%*[^\n]\n");
        }
        
        fclose(input_fp);
    }

    

    // No change needed - 120 just corresponds to the atomic numbers for known elements with some space at beginning and end
    int element_index_to_order[120] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int *element_index_to_order_p = element_index_to_order;
    for (int i=0; i<Ntypes; i++) element_index_to_order_p[elem[i]] = i; 

    mlff_str->element_index_to_order_p = (int *)malloc(120*sizeof(int));
    for (int i=0; i<120; i++) mlff_str->element_index_to_order_p[i] = element_index_to_order_p[i]; 



    // Array of atomic numbers
    // Hard coded 108 should be replaced with number of atoms from read_obj and the corresponding hard coded elemental number should be replaced 
    int atom_indices[read_obj->n_atom];
	int count = 0;
	for (int i=0; i < read_obj->Ntypes; i++){
		for (int j=0; j < read_obj->nAtomv[i]; j++){
			atom_indices[count] = elem[i];
			count++;
		}
	}


    // for (int i = 0; i < 108; i ++) atom_indices[i] = 13;
    int *atom_indices_p = atom_indices;
    // SPARC neighbor list seems to store atom type sequentially from 0. This converts SPARC atom type to atomic number.
    // Hard coding should be altered to match SPARC's readin of .ion file. 
    int atom_type_to_indices[read_obj->Ntypes];
    for (int i=0; i <read_obj->Ntypes; i++)
    	atom_type_to_indices[i] = elem[i];
    // int atom_type_to_indices[1] = {13};
    int *atom_type_to_indices_p = atom_type_to_indices;
    // Indices of atoms requiring descriptor calculation
    // This is needed for current factoring of GMP code. Hard coding should be replaced with information from read_obj object
    int cal_atoms[read_obj->n_atom];
    for (int i = 0; i < read_obj->n_atom; i++){
        cal_atoms[i] = i;
    }
    int *cal_atoms_p = cal_atoms;
    // Replace cal_num with read_obj information
    int cal_num = read_obj->n_atom;

    //FeatureScaler *ftr_scale;
    //ftr_scale = (FeatureScaler *) malloc(sizeof(FeatureScaler)*1);

	mlff_str->rcut = rcut;
	mlff_str->N_SIGMA = nsig;
	mlff_str->N_ANGULAR = nang;
	mlff_str->cal_atoms = (int *)malloc(read_obj->n_atom * sizeof(int));
	for (int i = 0; i < read_obj->n_atom; i++){
        mlff_str->cal_atoms[i] = cal_atoms[i];
    }
    printf("read_obj->n_atom %d\n",read_obj->n_atom);
    mlff_str->cal_num = read_obj->n_atom;

    mlff_str->params_i = (int**) malloc(sizeof(int*)*nmcsh);
    for (int i=0; i < nmcsh; i++)
    	mlff_str->params_i[i] = (int*) malloc(sizeof(int)*3);
   	for (int i=0; i < nang; i++){
   		for (int j=0; j<nsig; j++){
   			mlff_str->params_i[i*3+j][0] = params_i[i*3+j][0];
            mlff_str->params_i[i*3+j][1] = params_i[i*3+j][1];
            mlff_str->params_i[i*3+j][2] = params_i[i*3+j][2];
   		}
   	}

   	mlff_str->params_d = (double**) malloc(sizeof(double*)*nmcsh);
    for (int i = 0; i < nmcsh; i++)
    	mlff_str->params_d[i] = (double*) malloc(sizeof(double)*6);
	for (int i=0; i<nang; i++){
		for (int j=0; j < nsig; j++){
			mlff_str->params_d[i*3+j][0] = params_d[i*3+j][0];
            mlff_str->params_d[i*3+j][1] = params_d[i*3+j][1];
            mlff_str->params_d[i*3+j][2] = params_d[i*3+j][2];
            mlff_str->params_d[i*3+j][3] = params_d[i*3+j][3];
            mlff_str->params_d[i*3+j][4] = params_d[i*3+j][4];
            mlff_str->params_d[i*3+j][5] = params_d[i*3+j][5];
		}
	}

	mlff_str->nmcsh = nmcsh;

	mlff_str->atom_gaussians_p = (double**) malloc(sizeof(double*)*Ntypes);
    for (int i=0; i < Ntypes; i++)
    	mlff_str->atom_gaussians_p[i] = (double*) malloc(sizeof(double)*maxprim*2);
    for (int i=0; i < Ntypes; i++){
    	for (int j=0; j < ngaussians[i]; j++){
    		mlff_str->atom_gaussians_p[i][2*j] = atom_gaussians_p[i][2*j];
    		mlff_str->atom_gaussians_p[i][2*j+1] = atom_gaussians_p[i][2*j+1];
    	}
    }
    mlff_str->ngaussians = (int*) malloc(sizeof(int)*Ntypes);
    for (int i = 0; i < Ntypes; i++){
    	mlff_str->ngaussians[i] = ngaussians[i]; 
    }

	mlff_str->ftr_scale = (FeatureScaler *) malloc(sizeof(FeatureScaler)*1);
	(mlff_str->ftr_scale)->train = 1;

	mlff_str->atom_type_to_indices_p = (int *) malloc(read_obj->Ntypes*sizeof(int));
	for (int i=0; i <read_obj->Ntypes; i++)
		mlff_str->atom_type_to_indices_p[i] = atom_type_to_indices_p[i];

	mlff_str->atom_indices_p = (int *) malloc(read_obj->n_atom*sizeof(int));
	for (int i = 0; i < read_obj->n_atom; i++)
		mlff_str->atom_indices_p[i] = atom_indices_p[i];

	free(str_tmp);
	free(params_di);
	free(params_ii);
	free(pseudos);
	for (int i = 0; i < nmcsh; i++) free(params_d[i]);
	free(params_d);
	for (int i=0; i<nmcsh; i++) free(params_i[i]);
	free(params_i);
	for (int i=0; i<Ntypes; i++) free(atom_gaussians_p[i]);
	free(atom_gaussians_p);
}

void read_input(READ_Obj *read_obj) {
    char *input_filename = malloc(L_STRING * sizeof(char));
    char *str            = malloc(L_STRING * sizeof(char));
    char *temp           = malloc(L_STRING * sizeof(char));
    int i;
    int Flag_cell = 0;
    int Flag_latvec_scale = 0;
    snprintf(input_filename, L_STRING, "%s.inpt", read_obj->filename);
    
    FILE *input_fp = fopen(input_filename,"r");
    
    if (input_fp == NULL) {
        printf("\nCannot open file \"%s\"\n",input_filename);
        exit(EXIT_FAILURE);
    }
    
    while (!feof(input_fp)) {
        fscanf(input_fp,"%s",str);
        
        // enable commenting with '#'
        if (str[0] == '#' || str[0] == '\n'|| strcmpi(str,"undefined") == 0) {
            fscanf(input_fp, "%*[^\n]\n"); // skip current line
            continue;
        }

        // check variable name and assign value
        if (strcmpi(str,"CELL:") == 0) {
            Flag_cell = 1;
            fscanf(input_fp,"%lf", &read_obj->range_x);
            fscanf(input_fp,"%lf", &read_obj->range_y);
            fscanf(input_fp,"%lf", &read_obj->range_z);
            //printf("ranges: %lf,%lf,%lf\n",read_obj->range_x,read_obj->range_y,read_obj->range_z);
            fscanf(input_fp, "%*[^\n]\n");
        } else if (strcmpi(str,"BC:") == 0) { 
            // read BC in 1st lattice direction
            fscanf(input_fp,"%s",temp);
            if (strcmpi(temp,"p") == 0) {
                read_obj->BCx = 1;
            } else if (strcmpi(temp,"d") == 0) {
                read_obj->BCx = 0;
            } else {
                printf("Cannot recognize boundary condition: %s\n", temp);
                exit(EXIT_FAILURE);
            }
            // read BC in 2nd lattice direction
            fscanf(input_fp,"%s",temp);
            if (strcmpi(temp,"p") == 0) {
                read_obj->BCy = 1;
            } else if (strcmpi(temp,"d") == 0) {
                read_obj->BCy = 0;
            } else {
                printf("Cannot recognize boundary condition: %s\n", temp);
                exit(EXIT_FAILURE);
            }
            // read BC in 3rd lattice direction
            fscanf(input_fp,"%s",temp);
            if (strcmpi(temp,"p") == 0) {
                read_obj->BCz = 1;
            } else if (strcmpi(temp,"d") == 0) {
                read_obj->BCz = 0;
            } else {
                printf("Cannot recognize boundary condition: %s\n", temp);
                exit(EXIT_FAILURE);
            }
            snprintf(str, L_STRING, "undefined");    // initialize str
            //printf("BCs: %d,%d,%d",read_obj->BCx,read_obj->BCy,read_obj->BCz);
            fscanf(input_fp, "%*[^\n]\n");
        } else {
            printf("\nCannot recognize input variable identifier: \"%s\"\n",str);
            fscanf(input_fp, "%*[^\n]\n");
            //exit(EXIT_FAILURE);
        }
    }

    free(input_filename);
    free(str);
    free(temp);
    fclose(input_fp);
} 

/**
 * @brief   Read ion file.
 */
 //TODO: add n cell replica
void read_ion(READ_Obj *read_obj) {
    char *ion_filename = malloc(L_STRING * sizeof(char));
    char *str          = malloc(L_STRING * sizeof(char));
    int i, ityp, typcnt, atmcnt_coord, atmcnt_relax, atmcnt_spin, *atmcnt_cum, n_atom;
    snprintf(ion_filename, L_STRING, "%s.ion", read_obj->filename);
    FILE *ion_fp = fopen(ion_filename,"r");
    
    if (ion_fp == NULL) {
        printf("\nCannot open file \"%s\"\n",ion_filename);
        exit(EXIT_FAILURE);
    }
    
    /* first identify total number of atom types */
    typcnt = 0;
    while (!feof(ion_fp)) {
        fscanf(ion_fp,"%s",str);
        if (strcmpi(str, "ATOM_TYPE:") == 0) {
            typcnt++;
            fscanf(ion_fp, "%*[^\n]\n"); 
        } else {
            // skip current line
            fscanf(ion_fp, "%*[^\n]\n"); 
        }
    }
    
    if (typcnt < 1) {
        printf("\nPlease provide at least one type of atoms!\n");
        exit(EXIT_FAILURE);
    }
    read_obj->Ntypes = typcnt;
    //printf("Ntypes: %d",read_obj->Ntypes);
    
    // reset file pointer to the start of the file
    fseek(ion_fp, 0L, SEEK_SET);  // returns 0 if succeeded, can check status
    
    // allocate memory for nAtomv and atmcnt_cum
    atmcnt_cum = (int *)malloc( (read_obj->Ntypes+1) * sizeof(int));
    if (atmcnt_cum == NULL) {
        printf("\nCannot allocate memory for \"atmcnt_cum\"!\n");
        exit(EXIT_FAILURE);
    }
    
    // allocate memory 
    read_obj->atomType = (char *)calloc( read_obj->Ntypes * L_ATMTYPE, sizeof(char) ); 
    read_obj->nAtomv = (int *)malloc( read_obj->Ntypes * sizeof(int) );
    if (read_obj->atomType == NULL || read_obj->nAtomv == NULL) {
        printf("\nmemory cannot be allocated\n");
        exit(EXIT_FAILURE);
    }

    /* find total number of atoms */
    n_atom = 0;    // totoal num of atoms
    typcnt = -1;    // atom type count    
    atmcnt_cum[0] = 0;
    while (!feof(ion_fp)) {
        fscanf(ion_fp,"%s",str);
        if (strcmpi(str, "ATOM_TYPE:") == 0) {
            typcnt++;
            fscanf(ion_fp, "%s", &read_obj->atomType[L_ATMTYPE*typcnt]);
        } else if (strcmpi(str, "N_TYPE_ATOM:") == 0) {
            fscanf(ion_fp, "%d", &read_obj->nAtomv[typcnt]);
            //printf("natom_v: %d\n",read_obj->nAtomv[typcnt]);
            fscanf(ion_fp, "%*[^\n]\n"); 
            n_atom += read_obj->nAtomv[typcnt];
            atmcnt_cum[typcnt+1] = n_atom;
        } else {
            // skip current line
            fscanf(ion_fp, "%*[^\n]\n"); 
            continue;
        }
    }

    if (n_atom < 1) {
        printf("\nPlease provide at least one atom!\n");
        exit(EXIT_FAILURE);
    }    
    printf("n_atom from ion %d",n_atom);
    read_obj->n_atom = n_atom;
    //printf("n_atoms: %d\n",read_obj->n_atom);
    
    // allocate memory for atom positions, atom relax constraints and atom spin
    // currently no way to account for fractional coordinates - have to assume its one or the other. Print warning
    read_obj->atom_pos = (double *)malloc(3*n_atom*sizeof(double));
    if (read_obj->atom_pos == NULL) {
        printf("\nCannot allocate memory for atom positions, atom relax constraints and atom spin!\n");
        exit(EXIT_FAILURE);
    }
    
    // reset file pointer to the start of the file
    fseek(ion_fp, 0L, SEEK_SET);  // returns 0 if succeeded, can check status

    // reset temp var
    typcnt = -1;
    atmcnt_coord = 0;
    atmcnt_relax = 0;
    atmcnt_spin = 0;

    // allocate the size of the Isfrac vector which stores the coordinate type of each atom type
    read_obj->IsFrac = (int *)calloc( read_obj->Ntypes, sizeof(int) );
    
    // variables for checking number of inputs in a row
    int nums_read, array_read_int[10];
    double array_read_double[10];

    int read_flag;

    while (fscanf(ion_fp,"%s",str) != EOF) {
        // enable commenting with '#'
        if (str[0] == '#' || str[0] == '\n') {
            fscanf(ion_fp, "%*[^\n]\n"); // skip current line
            continue;
        }

        if (strcmpi(str, "ATOM_TYPE:") == 0) {
            typcnt++; 
            fscanf(ion_fp, "%*[^\n]\n"); // skip current line
        } else if (strcmpi(str, "N_TYPE_ATOM:") == 0) {
            fscanf(ion_fp, "%*[^\n]\n"); // skip current line
        } else if (strcmpi(str, "COORD:") == 0) {
            // fscanf(ion_fp, "%*[^\n]\n"); // skip current line
            check_below_entries(ion_fp, "COORD");
            //printf("coords:\n");
            for (i = 0; i < read_obj->nAtomv[typcnt]; i++) {
                nums_read = check_num_input(ion_fp, (void *) array_read_double, 'D');
                if (nums_read == -1) { i --; continue; }   // This is comment
                if (nums_read == 0) {
                    printf(RED "ERROR: Number of atom coordinates is less than number of atoms for atom type %d.\n" RESET, typcnt+1);
                    exit(EXIT_FAILURE);
                } else if (nums_read != 3)  { 
                    printf(RED "ERROR: please provide 3 coordinates on x y z for each atom of atom type %d in a row.\n" RESET, typcnt+1);
                    exit(EXIT_FAILURE);
                }
                read_obj->atom_pos[3*atmcnt_coord] = array_read_double[0];
                read_obj->atom_pos[3*atmcnt_coord+1] = array_read_double[1];
                read_obj->atom_pos[3*atmcnt_coord+2] = array_read_double[2];
                //printf("%lf,%lf,%lf\n",read_obj->atom_pos[3*atmcnt_coord],read_obj->atom_pos[3*atmcnt_coord+1],read_obj->atom_pos[3*atmcnt_coord+2]);
                atmcnt_coord++;  
            }
        } else if (strcmpi(str, "COORD_FRAC:") == 0) {
            // fscanf(ion_fp, "%*[^\n]\n"); // skip current line
            check_below_entries(ion_fp, "COORD_FRAC");
            //printf("coords:\n");
            for (i = 0; i < read_obj->nAtomv[typcnt]; i++) {
                nums_read = check_num_input(ion_fp, (void *) array_read_double, 'D');
                if (nums_read == -1) { i --; continue; }   // This is comment
                if (nums_read == 0) {
                    printf(RED "ERROR: Number of atom coordinates is less than number of atoms for atom type %d.\n" RESET, typcnt+1);
                    exit(EXIT_FAILURE);
                } else if (nums_read != 3)  { 
                    printf(RED "ERROR: please provide 3 coordinates on x y z for each atom of atom type %d in a row.\n" RESET, typcnt+1);
                    exit(EXIT_FAILURE);
                }
                read_obj->atom_pos[3*atmcnt_coord] = array_read_double[0];
                read_obj->atom_pos[3*atmcnt_coord+1] = array_read_double[1];
                read_obj->atom_pos[3*atmcnt_coord+2] = array_read_double[2];
                read_obj->atom_pos[3*atmcnt_coord] *= read_obj->range_x;
                read_obj->atom_pos[3*atmcnt_coord+1] *= read_obj->range_y;
                read_obj->atom_pos[3*atmcnt_coord+2] *= read_obj->range_z;
                //printf("%lf,%lf,%lf\n",read_obj->atom_pos[3*atmcnt_coord],read_obj->atom_pos[3*atmcnt_coord+1],read_obj->atom_pos[3*atmcnt_coord+2]);
                atmcnt_coord++;  
            }
            read_obj->IsFrac[typcnt] = 1;
        } else {
            printf("\nCannot recognize input variable identifier: \"%s\"\n",str);
            fscanf(ion_fp, "%*[^\n]\n"); 
            //exit(EXIT_FAILURE);
        }
    }
    
    if (atmcnt_coord != n_atom) {
        printf("the number of coordinates provided is inconsistent "
               "with the given number of atoms!\n");
        exit(EXIT_FAILURE);
    }

    free(atmcnt_cum);
    free(ion_filename);
    free(str);
    fclose(ion_fp);
}

void check_below_entries(FILE *ion_fp, char *tag) 
{
    int i;
    char *str = malloc(L_STRING * sizeof(char));

    str[0] = '\0';
    fscanf(ion_fp, "%[^\n]%*c", str);
    /*
    for (i = 0; i < strlen(str); i++) {
        if (isdigit(str[i])) {
            printf(RED"ERROR: Please remove the data in the same line as the %s tag in the ION file. All entries should be strictly below the %s tag.\n"RESET, tag, tag);
            exit(EXIT_FAILURE);
        }
    }
    */
    free(str);
}

int check_num_input(FILE *fp, void *array, char TYPE)
{
    int nums_now, bytes_now;
    int bytes_consumed = 0, nums_read = 0;
    char *str = malloc(L_STRING * sizeof(char));
    str[0] = '\0';

    if (TYPE != 'I' && TYPE != 'D') {
        printf("ERROR: Unknown type\n");
        exit(EXIT_FAILURE);
    }

    fscanf(fp,"%s",str);
    fseek ( fp , -strlen(str) , SEEK_CUR );
    
    if (str[0] == '#') {
        fscanf(fp, "%*[^\n]\n"); // skip current line
        free(str);
        return -1;
    }
    
    fscanf(fp, "%[^\n]%*c", str);

    if (TYPE == 'I') {
        while ( ( nums_now = 
                sscanf( str + bytes_consumed, "%d%n", (int *)array + nums_read, & bytes_now )
                ) > 0 && nums_read < 10) {
            bytes_consumed += bytes_now;
            nums_read += nums_now;
        }
    } else if (TYPE == 'D') {
        while ( ( nums_now = 
                sscanf( str + bytes_consumed, "%lf%n", (double *)array + nums_read, & bytes_now )
                ) > 0 && nums_read < 10) {
            bytes_consumed += bytes_now;
            nums_read += nums_now;
        }
    }
    
    free(str);
    return nums_read;
}

int strcmpi(const char *p1, const char *p2)
{
    const unsigned char *s1 = (const unsigned char *) p1;
    const unsigned char *s2 = (const unsigned char *) p2;
    unsigned char c1, c2;

    // in case p1 and p2 are pointed to the same string
    if (s1 == s2) return 0;

    do {
        // convert both strs to lower case first
        c1 = tolower(*s1++);
        c2 = tolower(*s2++);
        if (c1 == '\0') return c1 - c2;
    } while (c1 == c2);

    return c1 - c2;
}