#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>
#include "mlff_types.h"
#include "gmp_deriv_descriptor.h"

#define TEMP_TOL 1e-12
#define min(x,y) ((x)<(y)?(x):(y))
#define L_STRING 512
#define L_PSD 4096
#define MAX_PATH_LEN 10240
// M_PI needs to be defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    // Initialize variable in global scope and allocate memory for vars of known size
    char *str            = malloc(L_STRING * sizeof(char));
    char *str_tmp = malloc(MAX_PATH_LEN * sizeof(char));
    char *pseudos;
    double *params_di = (double*) malloc(sizeof(double)*4);
    int    *params_ii = (int*) malloc(sizeof(int)*2);
    int *elem, *ngaussians;
    int nsig,nang,Ntypes,maxprim;

    // Read GMP input file - this input is designed to be compatible with .inpt file for SPARC. Code is written analogously with some overlap
    FILE *input_fp = fopen("GMP.txt", "r");
    if (NULL == input_fp) printf("File cannont be opened\n");
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
                fscanf(input_fp, "%*[^\n]\n");
            }
        // Read number of unique elements - int
        } else if (strcmpi(str,"N_TYPES:") == 0){
            fscanf(input_fp,"%d",&Ntypes);
            pseudos = (char*) malloc(Ntypes*L_PSD*sizeof(char));
            fscanf(input_fp, "%*[^\n]\n");
        // Read array of pseudodensity files and store in character array
        } else if (strcmpi(str,"PSEUDODENSITY:") == 0){
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
            free(str_tmp);
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
                fscanf(input_fp, "%*[^\n]\n");
            }
        // Read max number of prmitives (could use built in function, but just an input for now) - int
        } else if (strcmpi(str,"MAX_PRIMITIVES:") == 0) {
            fscanf(input_fp,"%d",&maxprim);
            fscanf(input_fp, "%*[^\n]\n");
        }
    }
    fclose(input_fp);

    //Begin putting inputs from file into arrays
    // Store radial probes in sgma array
    double *sigmas = (double*) malloc(nsig*sizeof(double));
    for (int i=0; i<nsig; i++){
        sigmas[i] = params_di[i];
    }
    // Store cutoff
    double rcut = params_di[nsig];
    // Get descriptor vector size from nsig and nang probes
    int nmcsh = nsig*nang;
    // Allocate memory for GMP double params and populate 2D array
    double **params_d = (double**) malloc(sizeof(double*)*nmcsh);
    for (int i = 0; i < nmcsh; i++) params_d[i] = (double*) malloc(sizeof(double)*6);
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
    for (int i=0; i<nmcsh; i++) params_i[i] = (int*) malloc(sizeof(int)*3);
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
    for (int i=0; i<Ntypes; i++){
        FILE *input_fp = fopen(pseudos + i*L_PSD, "r");
        if (NULL == input_fp) printf("File cannont be opened\n");
        for (int j=0; j<ngaussians[i]; j++){
            fscanf(input_fp,"%lf%lf",&atom_gaussians_p[i][2*j],&atom_gaussians_p[i][2*j+1]);
        }
        fscanf(input_fp, "%*[^\n]\n");
        fclose(input_fp);
    }

    // No change needed - 120 just corresponds to the atomic numbers for known elements with some space at beginning and end
    int element_index_to_order[120] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int *element_index_to_order_p = element_index_to_order;
    for (int i=0; i<Ntypes; i++) element_index_to_order_p[elem[i]] = i; 

    // Manually input atom positions from Al supercell. Call pSPARC to generate this array or pass directly to GMP
    double atom_pos[108*3] = {0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,2.024999999999999911e+00,2.024999999999999911e+00,2.024999999999999911e+00,0.000000000000000000e+00,2.024999999999999911e+00,2.024999999999999911e+00,2.024999999999999911e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,4.049999999999999822e+00,0.000000000000000000e+00,2.024999999999999911e+00,6.074999999999999289e+00,2.024999999999999911e+00,0.000000000000000000e+00,6.074999999999999289e+00,2.024999999999999911e+00,2.024999999999999911e+00,4.049999999999999822e+00,0.000000000000000000e+00,0.000000000000000000e+00,8.099999999999999645e+00,0.000000000000000000e+00,2.024999999999999911e+00,1.012500000000000000e+01,2.024999999999999911e+00,0.000000000000000000e+00,1.012500000000000000e+01,2.024999999999999911e+00,2.024999999999999911e+00,8.099999999999999645e+00,0.000000000000000000e+00,4.049999999999999822e+00,0.000000000000000000e+00,0.000000000000000000e+00,6.074999999999999289e+00,2.024999999999999911e+00,2.024999999999999911e+00,4.049999999999999822e+00,2.024999999999999911e+00,2.024999999999999911e+00,6.074999999999999289e+00,0.000000000000000000e+00,0.000000000000000000e+00,4.049999999999999822e+00,4.049999999999999822e+00,0.000000000000000000e+00,6.074999999999999289e+00,6.074999999999999289e+00,2.024999999999999911e+00,4.049999999999999822e+00,6.074999999999999289e+00,2.024999999999999911e+00,6.074999999999999289e+00,4.049999999999999822e+00,0.000000000000000000e+00,4.049999999999999822e+00,8.099999999999999645e+00,0.000000000000000000e+00,6.074999999999999289e+00,1.012500000000000000e+01,2.024999999999999911e+00,4.049999999999999822e+00,1.012500000000000000e+01,2.024999999999999911e+00,6.074999999999999289e+00,8.099999999999999645e+00,0.000000000000000000e+00,8.099999999999999645e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.012500000000000000e+01,2.024999999999999911e+00,2.024999999999999911e+00,8.099999999999999645e+00,2.024999999999999911e+00,2.024999999999999911e+00,1.012500000000000000e+01,0.000000000000000000e+00,0.000000000000000000e+00,8.099999999999999645e+00,4.049999999999999822e+00,0.000000000000000000e+00,1.012500000000000000e+01,6.074999999999999289e+00,2.024999999999999911e+00,8.099999999999999645e+00,6.074999999999999289e+00,2.024999999999999911e+00,1.012500000000000000e+01,4.049999999999999822e+00,0.000000000000000000e+00,8.099999999999999645e+00,8.099999999999999645e+00,0.000000000000000000e+00,1.012500000000000000e+01,1.012500000000000000e+01,2.024999999999999911e+00,8.099999999999999645e+00,1.012500000000000000e+01,2.024999999999999911e+00,1.012500000000000000e+01,8.099999999999999645e+00,4.049999999999999822e+00,0.000000000000000000e+00,0.000000000000000000e+00,4.049999999999999822e+00,2.024999999999999911e+00,2.024999999999999911e+00,6.074999999999999289e+00,0.000000000000000000e+00,2.024999999999999911e+00,6.074999999999999289e+00,2.024999999999999911e+00,0.000000000000000000e+00,4.049999999999999822e+00,0.000000000000000000e+00,4.049999999999999822e+00,4.049999999999999822e+00,2.024999999999999911e+00,6.074999999999999289e+00,6.074999999999999289e+00,0.000000000000000000e+00,6.074999999999999289e+00,6.074999999999999289e+00,2.024999999999999911e+00,4.049999999999999822e+00,4.049999999999999822e+00,0.000000000000000000e+00,8.099999999999999645e+00,4.049999999999999822e+00,2.024999999999999911e+00,1.012500000000000000e+01,6.074999999999999289e+00,0.000000000000000000e+00,1.012500000000000000e+01,6.074999999999999289e+00,2.024999999999999911e+00,8.099999999999999645e+00,4.049999999999999822e+00,4.049999999999999822e+00,0.000000000000000000e+00,4.049999999999999822e+00,6.074999999999999289e+00,2.024999999999999911e+00,6.074999999999999289e+00,4.049999999999999822e+00,2.024999999999999911e+00,6.074999999999999289e+00,6.074999999999999289e+00,0.000000000000000000e+00,4.049999999999999822e+00,4.049999999999999822e+00,4.049999999999999822e+00,4.049999999999999822e+00,6.074999999999999289e+00,6.074999999999999289e+00,6.074999999999999289e+00,4.049999999999999822e+00,6.074999999999999289e+00,6.074999999999999289e+00,6.074999999999999289e+00,4.049999999999999822e+00,4.049999999999999822e+00,4.049999999999999822e+00,8.099999999999999645e+00,4.049999999999999822e+00,6.074999999999999289e+00,1.012500000000000000e+01,6.074999999999999289e+00,4.049999999999999822e+00,1.012500000000000000e+01,6.074999999999999289e+00,6.074999999999999289e+00,8.099999999999999645e+00,4.049999999999999822e+00,8.099999999999999645e+00,0.000000000000000000e+00,4.049999999999999822e+00,1.012500000000000000e+01,2.024999999999999911e+00,6.074999999999999289e+00,8.099999999999999645e+00,2.024999999999999911e+00,6.074999999999999289e+00,1.012500000000000000e+01,0.000000000000000000e+00,4.049999999999999822e+00,8.099999999999999645e+00,4.049999999999999822e+00,4.049999999999999822e+00,1.012500000000000000e+01,6.074999999999999289e+00,6.074999999999999289e+00,8.099999999999999645e+00,6.074999999999999289e+00,6.074999999999999289e+00,1.012500000000000000e+01,4.049999999999999822e+00,4.049999999999999822e+00,8.099999999999999645e+00,8.099999999999999645e+00,4.049999999999999822e+00,1.012500000000000000e+01,1.012500000000000000e+01,6.074999999999999289e+00,8.099999999999999645e+00,1.012500000000000000e+01,6.074999999999999289e+00,1.012500000000000000e+01,8.099999999999999645e+00,8.099999999999999645e+00,0.000000000000000000e+00,0.000000000000000000e+00,8.099999999999999645e+00,2.024999999999999911e+00,2.024999999999999911e+00,1.012500000000000000e+01,0.000000000000000000e+00,2.024999999999999911e+00,1.012500000000000000e+01,2.024999999999999911e+00,0.000000000000000000e+00,8.099999999999999645e+00,0.000000000000000000e+00,4.049999999999999822e+00,8.099999999999999645e+00,2.024999999999999911e+00,6.074999999999999289e+00,1.012500000000000000e+01,0.000000000000000000e+00,6.074999999999999289e+00,1.012500000000000000e+01,2.024999999999999911e+00,4.049999999999999822e+00,8.099999999999999645e+00,0.000000000000000000e+00,8.099999999999999645e+00,8.099999999999999645e+00,2.024999999999999911e+00,1.012500000000000000e+01,1.012500000000000000e+01,0.000000000000000000e+00,1.012500000000000000e+01,1.012500000000000000e+01,2.024999999999999911e+00,8.099999999999999645e+00,8.099999999999999645e+00,4.049999999999999822e+00,0.000000000000000000e+00,8.099999999999999645e+00,6.074999999999999289e+00,2.024999999999999911e+00,1.012500000000000000e+01,4.049999999999999822e+00,2.024999999999999911e+00,1.012500000000000000e+01,6.074999999999999289e+00,0.000000000000000000e+00,8.099999999999999645e+00,4.049999999999999822e+00,4.049999999999999822e+00,8.099999999999999645e+00,6.074999999999999289e+00,6.074999999999999289e+00,1.012500000000000000e+01,4.049999999999999822e+00,6.074999999999999289e+00,1.012500000000000000e+01,6.074999999999999289e+00,4.049999999999999822e+00,8.099999999999999645e+00,4.049999999999999822e+00,8.099999999999999645e+00,8.099999999999999645e+00,6.074999999999999289e+00,1.012500000000000000e+01,1.012500000000000000e+01,4.049999999999999822e+00,1.012500000000000000e+01,1.012500000000000000e+01,6.074999999999999289e+00,8.099999999999999645e+00,8.099999999999999645e+00,8.099999999999999645e+00,0.000000000000000000e+00,8.099999999999999645e+00,1.012500000000000000e+01,2.024999999999999911e+00,1.012500000000000000e+01,8.099999999999999645e+00,2.024999999999999911e+00,1.012500000000000000e+01,1.012500000000000000e+01,0.000000000000000000e+00,8.099999999999999645e+00,8.099999999999999645e+00,4.049999999999999822e+00,8.099999999999999645e+00,1.012500000000000000e+01,6.074999999999999289e+00,1.012500000000000000e+01,8.099999999999999645e+00,6.074999999999999289e+00,1.012500000000000000e+01,1.012500000000000000e+01,4.049999999999999822e+00,8.099999999999999645e+00,8.099999999999999645e+00,8.099999999999999645e+00,8.099999999999999645e+00,1.012500000000000000e+01,1.012500000000000000e+01,1.012500000000000000e+01,8.099999999999999645e+00,1.012500000000000000e+01,1.012500000000000000e+01,1.012500000000000000e+01,8.099999999999999645e+00};

    double *atom_pos_p = atom_pos;
    // Array of atomic numbers
    // Hard coded 108 should be replaced with number of atoms from pSPARC and the corresponding hard coded elemental number should be replaced 
    int atom_indices[108];
    for (int i = 0; i < 108; i ++) atom_indices[i] = 13;
    int *atom_indices_p = atom_indices;
    // SPARC neighbor list seems to store atom type sequentially from 0. This converts SPARC atom type to atomic number.
    // Hard coding should be altered to match SPARC's readin of .ion file. 
    int atom_type_to_indices[1] = {13};
    int *atom_type_to_indices_p = atom_type_to_indices;
    // Indices of atoms requiring descriptor calculation
    // This is needed for current factoring of GMP code. Hard coding should be replaced with information from pSPARC object
    int cal_atoms[108];
    for (int i = 0; i < 108; i++){
        cal_atoms[i] = i;
    }
    int *cal_atoms_p = cal_atoms;
    // Replace cal_num with pSPARC information
    int cal_num = 108;

    NeighList *nlist;
	GMPObj *gmp_str;
    FeatureScaler *ftr_scale;
	int *Z, *BC, *atomtyp;
    // atoms in system and ordered list of number of atoms corresponding to each atom type. Replace with SPARC input
    int n_atom = 108, nAtomv[1] = {108};
	Z = (int *)malloc(Ntypes*sizeof(int));
	for (int i=0; i <Ntypes; i++){
		Z[i] = 1;//pSPARC->Znucl[i];
	}

	nlist = (NeighList *) malloc(sizeof(NeighList)*1);
	gmp_str = (GMPObj *) malloc(sizeof(GMPObj)*1);
    ftr_scale = (FeatureScaler *) malloc(sizeof(FeatureScaler)*1);

	double *cell, F_estimate_max;
    // BC - boundary conditions 1 = periodic
	BC = (int*)malloc(3*sizeof(int));
    //cell - unit cell
	cell =(double *) malloc(3*sizeof(double));
    // atomtyp - array storing atoms types 
	atomtyp = (int *) malloc(n_atom*sizeof(int));
    // Boundary conditions and cell params - should be readable from SPARC
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

	initialize_nlist(nlist, n_atom, rcut, Ntypes);
    //printf("%d,%lf,%d",n_atom,rcut,Ntypes);
	build_nlist(rcut, Ntypes, n_atom, atom_pos_p, atomtyp, BC, cell, nlist);

    // specifying if descriptor calculation is for training or prediction - determines if scaling factors are computed for standardization
    // This requires smarter implementation
    ftr_scale->train = 1;

    initialize_gmpObj(gmp_str, nlist, cal_atoms, cal_num,params_i, params_d, nmcsh, atom_gaussians_p, ngaussians, element_index_to_order_p);
    build_gmpObj(gmp_str, nlist, ftr_scale, nmcsh, atom_pos_p, params_i, params_d, atom_gaussians_p, ngaussians, element_index_to_order_p,atom_type_to_indices_p,atom_indices_p);
    
    for (int i=0; i<n_atom; i++){
        for (int j=0; j < 15; j++){
            printf("%.10f,",gmp_str->X[i][j]);
        
            for (int k=0; k<n_atom; k++){
                printf("%.10f,",gmp_str->dX_dX[i][k][j]);
            }
            printf("\n");
        }        
    }
    
    return 0;
}

int strcmpi(const char *p1, const char *p2){
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