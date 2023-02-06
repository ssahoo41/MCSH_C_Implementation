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
#include "linearsys.h"
#include "sparsification.h"
#include "ddbp_tools.h"
#include "tools.h"
#include "regression.h"
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))


/*
mlff_train_SVD solves the linear system in training using SVD

[Input]
1. mlff_str: MLFF_Obj structure
[Output]
1. mlff_str: MLFF_Obj structure
*/

void mlff_train_SVD(MLFF_Obj *mlff_str){

	int info, m = mlff_str->n_rows, n = mlff_str->n_cols, lda=mlff_str->n_rows;
    double ldu=mlff_str->n_rows, ldvt=mlff_str->n_cols;
    int i, j,k;
	double superb[min(m,n)-1];
	double *s, *u, *vt;
	double *a_scaled, *a_scaled1, *b_scaled, *b_scaled1;

    ///////////////////////////////////////////////////////////////////////////////////
                                //  Scaling applied to design matrices
    ///////////////////////////////////////////////////////////////////////////////////
    mlff_str->E_scale = 1.0;
    mlff_str->F_scale = mlff_str->std_E/mlff_str->std_F;
    mlff_str->stress_scale[0] = mlff_str->std_E/mlff_str->std_stress[0];
    mlff_str->stress_scale[1] = mlff_str->std_E/mlff_str->std_stress[1];
    mlff_str->stress_scale[2] = mlff_str->std_E/mlff_str->std_stress[2];
    mlff_str->stress_scale[3] = mlff_str->std_E/mlff_str->std_stress[3];
    mlff_str->stress_scale[4] = mlff_str->std_E/mlff_str->std_stress[4];
    mlff_str->stress_scale[5] = mlff_str->std_E/mlff_str->std_stress[5];


	a_scaled = (double *) malloc(lda*n * sizeof(double));
    a_scaled1 = (double *) malloc(lda*n * sizeof(double));
    b_scaled = (double *) malloc(lda * sizeof(double));
    b_scaled1 = (double *) malloc(lda * sizeof(double));
	s = (double *) malloc(n * sizeof(double));
	// u = (double *) malloc(ldu*m * sizeof(double));
	// vt = (double *) malloc(ldvt*n * sizeof(double));

    //int isE=0;
    //int isF=0;
    //int isStress = {0,0,0,0,0,0};
    double scale =0;
	for (i = 0; i < m; i++ ){
        int quot = i%(7+3*mlff_str->natom);
        if (quot==0){
            scale = mlff_str->E_scale;
            b_scaled[i] = (1.0/mlff_str->std_E)*(mlff_str->b_no_norm[i] - mlff_str->mu_E);
        } else if (quot>0 && quot <= 3*mlff_str->natom){
            scale = mlff_str->F_scale;
            b_scaled[i] = (1.0/mlff_str->std_F)*(mlff_str->b_no_norm[i]);
        }else{
            scale = mlff_str->stress_scale[quot-3*mlff_str->natom-1];
            b_scaled[i] = (1.0/mlff_str->std_stress[quot-3*mlff_str->natom-1])*(mlff_str->b_no_norm[i]);
        }
        b_scaled1[i] = b_scaled[i];//mlff_str->b_no_norm[i];
        // scale=1.0;
		for (j = 0; j < n; j++){
			a_scaled[j*m+i] = scale * mlff_str->K_train[i][j];
            a_scaled1[j*m+i] = a_scaled[j*m+i];
		}
	}


    //////////////////////////////////
   //  double acheck[50]={0.337122644398882,   0.748151592823709,   0.442678269775446,   0.800068480224308,   0.144954798223727,
   // 0.162182308193243,   0.450541598502498,   0.106652770180584,   0.431413827463545,   0.853031117721894,
   // 0.794284540683907,   0.083821377996933,   0.961898080855054,   0.910647594429523,   0.622055131485066,
   // 0.311215042044805,   0.228976968716819,   0.004634224134067,   0.181847028302852,   0.350952380892271,
   // 0.528533135506213,   0.913337361501670,   0.774910464711502,   0.263802916521990,   0.513249539867053,
   // 0.165648729499781,   0.152378018969223,   0.817303220653433,   0.145538980384717,   0.401808033751942,
   // 0.601981941401637,   0.825816977489547,   0.868694705363510,   0.136068558708664,   0.075966691690842,
   // 0.262971284540144,   0.538342435260057,   0.084435845510910,   0.869292207640089,   0.239916153553658,
   // 0.654079098476782,   0.996134716626885,   0.399782649098896,   0.579704587365570,   0.123318934835166,
   // 0.689214503140008,   0.078175528753184,   0.259870402850654,   0.549860201836332,   0.183907788282417};
   // double acheck[50]={0.337122644398882,   0.162182308193243,   0.794284540683907,   0.311215042044805,   0.528533135506213,   0.165648729499781,   0.601981941401637,   0.262971284540144,   0.654079098476782,   0.689214503140008,
   // 0.748151592823709,   0.450541598502498,   0.083821377996933,   0.228976968716819,   0.913337361501670,   0.152378018969223,   0.825816977489547,   0.538342435260057,   0.996134716626885,   0.078175528753184,
   // 0.442678269775446,   0.106652770180584,   0.961898080855054,   0.004634224134067,   0.774910464711502,   0.817303220653433,   0.868694705363510,   0.084435845510910,   0.399782649098896,   0.259870402850654,
   // 0.800068480224308,   0.431413827463545,   0.910647594429523,   0.181847028302852,   0.263802916521990,   0.145538980384717,   0.136068558708664,   0.869292207640089,   0.579704587365570,  0.549860201836332,
   // 0.144954798223727,   0.853031117721894,   0.622055131485066,   0.350952380892271,   0.513249539867053,   0.401808033751942,   0.075966691690842,   0.239916153553658,   0.123318934835166,   0.183907788282417};
   // double scheck[5], ucheck[100], vtcheck[25];
   // info = LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'A', 'A', 10, 5, acheck, 10,
   //                      scheck, ucheck, 10, vtcheck, 5, superb);
   // printf("singular values\n");
   // for (int i=0; i <5; i++){
   //  printf("%f\n",scheck[i]);
   // }
   // printf("U:\n");
   // for (int i=0; i <10; i++){
   //  for (int j=0; j <10; j++){
   //      printf("%f ",ucheck[i*10+j]);
   //  }
   //  printf("\n");
   // }
   // printf("V':\n");
   // for (int i=0; i <5; i++){
   //  for (int j=0; j <5; j++){
   //      printf("%f ",vtcheck[i*5+j]);
   //  }
   //  printf("\n");
   // }

   // double Utkcheck[50], Vkcheck[25], sikcheck[5];
   //  for (k = 0; k < 5; k++){
   //      sikcheck[k] = 1.0/scheck[k];
   //      for (i = 0; i < 5; i++){
   //          Vkcheck[k*5+i] = vtcheck[k+5*i]*sikcheck[k];
   //      }
   //  }
   //  for (i = 0; i < 10; i++){
   //      for (k = 0; k < 5; k++){
   //          Utkcheck[i*5+k] = ucheck[k*10+i];
   //      }
   //  }

   //  double *A_pinvcheck;

   //  A_pinvcheck = (double *) malloc(sizeof(double)* 50);

    
   //  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
   //              5, 10, 5, 1.0, Vkcheck, 5, Utkcheck, 5, 0.0, A_pinvcheck, 5);
   //  printf("pinv: \n");
   //  for (k = 0; k < 10; k++){
   //      for (i = 0; i < 5; i++){
   //          printf("%f ",A_pinvcheck[k*5+i]);
   //      }
   //      printf("\n");
   //  }
   //  exit(0);
    /////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
								//	SVD of design matrix
	///////////////////////////////////////////////////////////////////////////////////
	// info = LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'A', 'A', m, n, a_scaled1, lda,
 //                        s, u, ldu, vt, ldvt, superb );

	// /* Check for convergence */
 //        if( info > 0 ) {
 //                printf( "The algorithm computing SVD failed to converge.\n" );
 //                exit( 1 );
 //        }

 //    // Inverting the diagonal entries
 //    int rank=0;
 //    printf("svd values:\n");
 //    for (i =0; i<n; i++){
 //        printf("%10.16f\n",s[i]);
 //    	if (s[i]>1e-5){
 //    		rank++;
 //    	}
 //    }
 //    printf("rank of the design matrix: %d\n", rank);
 //    ///////////////////////////////////////////////////////////////////////////////////
 //    				// truncated SVD left and right singular vectors
 //    ///////////////////////////////////////////////////////////////////////////////////
 //    double Utk[rank*m], Vk[rank*n], sik[rank];
 //    for (k = 0; k < rank; k++){
 //    	sik[k] = 1.0/s[k];
 //    	for (i = 0; i < n; i++){
 //    		Vk[i+k*n] = vt[k+i*n]*sik[k];
 //    	}
 //    }
 //    for (i = 0; i < m; i++){
 //    	for (k = 0; k < rank; k++){
 //    		Utk[k+i*rank] = u[i+k*m];
 //    	}
 //    }
 //    ///////////////////////////////////////////////////////////////////////////////////
 //    					// V*Si*U' to calculate the pinv
 //    ///////////////////////////////////////////////////////////////////////////////////
 //    double *A_pinv;

 //    A_pinv = (double *) malloc(sizeof(double)* m*n);

 //    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
 //                n, m, rank, 1.0, Vk, n, Utk, rank, 0.0, A_pinv, n);

 //    FILE *Apinv;
 //    Apinv = fopen("A_pinv.txt","w");
 //    for (int i=0; i <n;i++){
 //        for (int j=0; j<m; j++){
 //            fprintf(Apinv,"%10.9f ",A_pinv[j*n+i]);
 //        }
 //        fprintf(Apinv,"\n");
 //    }
 //    fclose(Apinv);
 //    // weights by pinv(A) * b
 //    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
 //    			n, 1, m, 1.0, A_pinv, n, b_scaled, m, 0.0, mlff_str->weights, n);

    ///////////////////////////////////////////////////////////////////////////////////
                        // weights using SVD least square
    ///////////////////////////////////////////////////////////////////////////////////
    // FILE *fid;
    // fid = fopen("bvec.txt","w");
    // for (int i=0; i<m; i++){
    //     fprintf(fid,"%10.9f\n",b_scaled1[i]);
    // }
    // fclose(fid);

    int rank;
    LAPACKE_dgelsd(LAPACK_COL_MAJOR, m, n, 1, a_scaled1, m, b_scaled1, m, s, -1.0, &rank );
    printf("rank of the design matrix in least square: %d\n", rank);
    for (int i=0; i <n; i++){
        mlff_str->weights[i] = b_scaled1[i];
        // printf("%10.9f\n",mlff_str->weights[i]);
    }



    ///////////////////////////////////////////////////////////////////////////////////
    					// hyperparameter optimization
    ///////////////////////////////////////////////////////////////////////////////////

    double sigma_w0=mlff_str->sigma_w, sigma_v0=mlff_str->sigma_v, sigma_w, sigma_v;
    double error_w = 1, error_v=1, error_tol = 0.001;
    int count=0;

    //<calculate A'A>
    // <Q0, lambda_0>

    double *AtA, *AtA0;

    AtA = (double *) malloc(sizeof(double)* mlff_str->n_cols * mlff_str->n_cols);
    AtA0 = (double *) malloc(sizeof(double)* mlff_str->n_cols * mlff_str->n_cols);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
    			mlff_str->n_cols, mlff_str->n_cols, mlff_str->n_rows, 1.0, a_scaled, mlff_str->n_rows, a_scaled, mlff_str->n_rows, 0.0, AtA, mlff_str->n_cols);

    for (i = 0; i < mlff_str->n_cols * mlff_str->n_cols ; i++){
    	AtA0[i] = AtA[i];
        // printf("%10.9f %10.9f\n",AtA[i],AtA0[i]);
    }

    double *lambda_0, *lambda, gamma;

    lambda_0 = (double *) malloc(sizeof(double)* mlff_str->n_cols);
    lambda = (double *) malloc(sizeof(double)* mlff_str->n_cols);

    info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', mlff_str->n_cols, AtA, mlff_str->n_cols, lambda_0);


    double *e, *Ab;
    Ab = (double *) malloc(sizeof(double)*mlff_str->n_rows);
    e = (double *) malloc(sizeof(double)*mlff_str->n_rows);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
    			m, 1, n, 1.0, a_scaled, m, mlff_str->weights, n, 0.0, Ab, m);
    

    for (i = 0; i < mlff_str->n_rows; i++){
    	e[i] = Ab[i]-b_scaled[i];
    }


    double norm_er2 = dotProduct(e, e, mlff_str->n_cols);


    double norm_w = sqrt(dotProduct(mlff_str->weights, mlff_str->weights, mlff_str->n_cols));

    printf("sigma_v0: %f, sigma_w0: %f\n",sigma_v0,sigma_w0);
    while( (error_w > error_tol || error_v > error_tol) && count <100) {
    	gamma = 0;
    	for (k = 0; k < mlff_str->n_cols; k++){
    		lambda[k] = lambda_0[k]/(sigma_v0*sigma_v0);
    		if (lambda[k] > 0.0000000001){
    			gamma += lambda[k]/(lambda[k] + 1/(sigma_w0*sigma_w0));
    		}
    	}

    	sigma_w = sqrt(norm_w*norm_w/gamma);
    	sigma_v = sqrt(norm_er2 /(mlff_str->n_rows-gamma));
    	error_w = fabs(sigma_w-sigma_w0);
    	error_v = fabs(sigma_v-sigma_v0);
    	sigma_w0 = sigma_w;
    	sigma_v0 = sigma_v;
        printf("sigma_v: %f, sigma_w: %f\n",sigma_v,sigma_w);

      count++;
   }
    double *Q0_temp, *S;
    Q0_temp = (double *) malloc(sizeof(double)*mlff_str->n_cols*mlff_str->n_cols);
    S = (double *) malloc(sizeof(double)*mlff_str->n_cols*mlff_str->n_cols);


    for (i = 0; i < mlff_str->n_cols; i++){
    	for (j = 0; j < mlff_str->n_cols; j++){
    		Q0_temp[i*mlff_str->n_cols + j] = lambda_0[i]*(1/(1+(sigma_w0/sigma_v0)*(sigma_w0/sigma_v0)))*AtA[i*mlff_str->n_cols + j];
    	}
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 
    			mlff_str->n_cols, mlff_str->n_cols, mlff_str->n_cols, sigma_w*sigma_w, Q0_temp, mlff_str->n_cols, AtA, mlff_str->n_cols, 0.0, S, mlff_str->n_cols);
    
    free(mlff_str->cov_train);
    mlff_str->cov_train = (double *)malloc(mlff_str->n_cols*mlff_str->n_cols*sizeof(double));
    for (i = 0; i < mlff_str->n_cols; i++){
    	for (j = 0; j < mlff_str->n_cols; j++){
    		mlff_str->cov_train[i*mlff_str->n_cols+j] = AtA[i*mlff_str->n_cols+j];
    	}
    }
    if (error_w<0.001 &&  error_v <0.001){
	    printf("sigma_v and sigma_w converged in %d iteration\n",count);
    }
	else{
	    printf("sigma_v and sigma_w did not converge and error was %10.9f, %10.9f\n",error_v, error_w);
	}
	free(s);
	// free(u);
	// free(vt);
	free(a_scaled);
  free(a_scaled1);
  free(b_scaled);
  free(b_scaled1);
	// free(A_pinv);
	free(AtA);
  free(AtA0);
  free(lambda_0);
	free(lambda);
	free(Q0_temp);
	free(S);
  free(e);
  free(Ab);
}

void mlff_train_Bayesian(MLFF_Obj *mlff_str){

    int info, m = mlff_str->n_rows, n = mlff_str->n_cols;
    int i,j,k, count, count1;
    // FILE *f1, *f2;

    /////////////////////  CUR on whole dataset //////////////////////////////////////

    double **X2_cur, **X3_cur;
    int *org_idx;
    int ncols_to_remove;
    int *cols_to_remove;
    dyArray *highrank_ID_descriptors;
    highrank_ID_descriptors = (dyArray *) malloc(sizeof(dyArray)*mlff_str->nelem);

    for (int i=0; i < mlff_str->nelem; i++){
      init_dyarray(&highrank_ID_descriptors[i]);

      X2_cur = (double **) malloc(sizeof(double*)*mlff_str->natm_train_elemwise[i]);
      X3_cur = (double **) malloc(sizeof(double*)*mlff_str->natm_train_elemwise[i]);
      org_idx = (int *)malloc(sizeof(int)*mlff_str->natm_train_elemwise[i]);

      for (int j=0; j < mlff_str->natm_train_elemwise[i]; j++){
        X2_cur[j] = (double *) malloc(sizeof(double)*mlff_str->size_X2);
        X3_cur[j] = (double *) malloc(sizeof(double)*mlff_str->size_X3);
      }

      count=0;
      for (int j=0; j < mlff_str->natm_train_total; j++){
        if (mlff_str->natm_typ_train[j]==i){
          org_idx[count] = j;
          for (int k=0; k<mlff_str->size_X2; k++){
            // printf("i: %d, j: %d, count: %d, k: %d, mlff_str->X2_traindataset[j][k]: %f\n",i,j,count, k, mlff_str->X2_traindataset[j][k]);
            X2_cur[count][k] = mlff_str->X2_traindataset[j][k];

          }
          for (int k=0; k<mlff_str->size_X3; k++){
            // printf("i: %d, j: %d, count: %d, k: %d, mlff_str->X3_traindataset[j][k]: %f\n",i,j,count, k, mlff_str->X3_traindataset[j][k]);
            X3_cur[count][k] = mlff_str->X3_traindataset[j][k];
          }
          count++;
        }
      }

      if (count != mlff_str->natm_train_elemwise[i]){
        printf("something wrong in CUR in mlff_train_Bayesian\n");
        exit(1);
      }

      SOAP_CUR_sparsify(mlff_str->kernel_typ, X2_cur, X3_cur, mlff_str->natm_train_elemwise[i], mlff_str->size_X2, mlff_str->size_X3, mlff_str->beta_2, mlff_str->beta_3, mlff_str->xi_3, &highrank_ID_descriptors[i]);
      ncols_to_remove = mlff_str->natm_train_elemwise[i] - (&highrank_ID_descriptors[i])->len;
      printf("For element %d, CUR removed %d columns out of a total of %d columns\n",i,ncols_to_remove,mlff_str->natm_train_elemwise[i]);
      cols_to_remove = (int *) malloc(sizeof(int)*ncols_to_remove);

      count1=0;
      print_dyarray(&highrank_ID_descriptors[i]);
      for (int l=0; l < mlff_str->natm_train_elemwise[i]; l++){
        count=0;
        for (int k = 0; k <(&highrank_ID_descriptors[i])->len; k++){
          if (l==(&highrank_ID_descriptors[i])->array[k]){
            count++;
          }
        }
        if (count==0){
          cols_to_remove[count1] = org_idx[l];
          count1++;
        }

      }
      printf("\n");
      if (count1 != ncols_to_remove){
        printf("something wrong in CUR! \n");
        exit(1);
      }

      for (int k=0; k < ncols_to_remove; k++){
        remove_train_cols(mlff_str, cols_to_remove[k]);
      }
      
      delete_dyarray(&highrank_ID_descriptors[i]);

      for (int j=0; j < mlff_str->natm_train_elemwise[i]; j++){
        free(X2_cur[j]);
        free(X3_cur[j]);
      }

      free(X2_cur);
      free(X3_cur);
      free(org_idx);
      free(cols_to_remove);

    }

    free(highrank_ID_descriptors);
    /////////////////////  CUR on whole dataset //////////////////////////////////////

    double *a_scaled, *b_scaled;

    ///////////////////////////////////////////////////////////////////////////////////
                                //  Scaling applied to design matrices
    ///////////////////////////////////////////////////////////////////////////////////
    mlff_str->E_scale = 1.0;
    mlff_str->F_scale = mlff_str->std_E/mlff_str->std_F;
    mlff_str->stress_scale[0] = mlff_str->std_E/mlff_str->std_stress[0];
    mlff_str->stress_scale[1] = mlff_str->std_E/mlff_str->std_stress[1];
    mlff_str->stress_scale[2] = mlff_str->std_E/mlff_str->std_stress[2];
    mlff_str->stress_scale[3] = mlff_str->std_E/mlff_str->std_stress[3];
    mlff_str->stress_scale[4] = mlff_str->std_E/mlff_str->std_stress[4];
    mlff_str->stress_scale[5] = mlff_str->std_E/mlff_str->std_stress[5];


    a_scaled = (double *) malloc(m*n * sizeof(double));
    b_scaled = (double *) malloc(m * sizeof(double));

    double scale =0;
    for (i = 0; i < m; i++ ){
        int quot = i%(7+3*mlff_str->natom);
        if (quot==0){
            scale = mlff_str->E_scale ;
            b_scaled[i] = (1.0/mlff_str->std_E)*(mlff_str->b_no_norm[i] - mlff_str->mu_E);
        } else if (quot>0 && quot <= 3*mlff_str->natom){
            scale = mlff_str->F_scale* mlff_str->relative_scale_F;
            b_scaled[i] = (1.0/mlff_str->std_F)*(mlff_str->b_no_norm[i])* mlff_str->relative_scale_F;
        }else{
            scale = mlff_str->stress_scale[quot-3*mlff_str->natom-1]* mlff_str->relative_scale_stress[quot-3*mlff_str->natom-1];
            b_scaled[i] = (1.0/mlff_str->std_stress[quot-3*mlff_str->natom-1])*(mlff_str->b_no_norm[i])* mlff_str->relative_scale_stress[quot-3*mlff_str->natom-1];
        }

        for (j = 0; j < n; j++){
            a_scaled[j*m+i] = scale * mlff_str->K_train[i][j];
        }
    }

    double *AtA, *Atb, *AtA_h, *Atb_h;

    AtA = (double *) malloc(sizeof(double)* mlff_str->n_cols * mlff_str->n_cols);
    Atb = (double *) malloc(sizeof(double)* mlff_str->n_cols);

    AtA_h = (double *) malloc(sizeof(double)* mlff_str->n_cols * mlff_str->n_cols);
    Atb_h = (double *) malloc(sizeof(double)* mlff_str->n_cols);

    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
                mlff_str->n_cols, mlff_str->n_cols, mlff_str->n_rows, 1.0, a_scaled, mlff_str->n_rows, a_scaled, mlff_str->n_rows, 0.0, AtA, mlff_str->n_cols);

    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
                mlff_str->n_cols, 1, mlff_str->n_rows, 1.0, a_scaled, mlff_str->n_rows, b_scaled, mlff_str->n_rows, 0.0, Atb, mlff_str->n_cols);

   
    // FILE *f1, *f2;
    // f1 = fopen("a_scaled.txt","w");
    // f2 = fopen("b_scaled.txt","w");

    // for (i = 0; i < m; i++){
    //   fprintf(f2,"%10.9f\n",b_scaled[i]);
    //   for (j = 0; j < n; j++){
    //     fprintf(f1,"%10.9f ",a_scaled[j*m+i]);
    //   }
    //   fprintf(f1,"\n");
    // }
    // fclose(f1);
    // fclose(f2);


    for (int i=0; i < mlff_str->n_cols * mlff_str->n_cols; i++){
    	AtA_h[i] = AtA[i];
    }

    for (int i=0; i < mlff_str->n_cols; i++){
    	Atb_h[i] = Atb[i];
    }

    int dohyperparameter=1;
    double sigma_w, sigma_v;
    if (dohyperparameter){
    	hyperparameter_Bayesian(a_scaled, AtA_h, Atb_h, b_scaled, mlff_str);
    }

    for (int i =0; i <mlff_str->n_cols; i++){  
      AtA[i*mlff_str->n_cols+i] += (mlff_str->sigma_v*mlff_str->sigma_v)/(mlff_str->sigma_w*mlff_str->sigma_w);

    }



    int ipiv[mlff_str->n_cols];
    // info = LAPACKE_dgesv( LAPACK_COL_MAJOR, mlff_str->n_cols, 1, AtA, mlff_str->n_cols, &ipiv[0],
    //                      Atb, mlff_str->n_cols );

    // for (int i=0; i < mlff_str->n_cols; i++)
    // 	mlff_str->weights[i] = Atb[i];

    free(mlff_str->cov_train);
    mlff_str->cov_train = (double *)malloc(mlff_str->n_cols*mlff_str->n_cols*sizeof(double));

    // FILE *f1, *f2;
    // f1 = fopen("AtA.txt","w");
    // f2 = fopen("AtAinv.txt","w");
    // for (int i=0; i < mlff_str->n_cols; i++){
    //   for (int j=0; j <mlff_str->n_cols; j++){
    //     fprintf(f1,"%f ",AtA[j*mlff_str->n_cols + i]);
    //   }
    //   fprintf(f1,"\n");
    // }
    // fclose(f1);
    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, mlff_str->n_cols, mlff_str->n_cols, AtA, mlff_str->n_cols, &ipiv[0]);
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, mlff_str->n_cols, AtA, mlff_str->n_cols, &ipiv[0]);

    // for (int i=0; i < mlff_str->n_cols; i++){
    //   for (int j=0; j <mlff_str->n_cols; j++){
    //     fprintf(f2,"%f ",AtA[j*mlff_str->n_cols + i]);
    //   }
    //   fprintf(f2,"\n");
    // }
    // fclose(f2);

    // exit(1);

    for (int i = 0; i < mlff_str->n_cols*mlff_str->n_cols; i++){
    	mlff_str->cov_train[i] = AtA[i] * mlff_str->sigma_v*mlff_str->sigma_v;
    }


    free(a_scaled);
    free(b_scaled);
    free(AtA);
    free(Atb);
    free(AtA_h);
    free(Atb_h);
}

void hyperparameter_Bayesian(double *A, double *AtA, double *Atb, double *b, MLFF_Obj *mlff_str){


	double sigma_v0, sigma_w0;

	sigma_v0 = mlff_str->sigma_v;
	sigma_w0 = mlff_str->sigma_w;

	double *lambda_0, *lambda;

  lambda_0 = (double *) malloc(sizeof(double)* mlff_str->n_cols);
  lambda = (double *) malloc(sizeof(double)* mlff_str->n_cols);

	double *AtA_h, *Atb_h, *weights, *Ax, *e;

	AtA_h = (double *) malloc(mlff_str->n_cols*mlff_str->n_cols*sizeof(double));
	Atb_h = (double *) malloc(mlff_str->n_cols*sizeof(double));
	weights = (double *) malloc(mlff_str->n_cols*sizeof(double));
  Ax = (double *) malloc(sizeof(double)*mlff_str->n_rows);
  e = (double *) malloc(sizeof(double)*mlff_str->n_rows);

	for (int i =0; i <mlff_str->n_cols; i++){
		Atb_h[i] = Atb[i];
        for (int j=0; j <mlff_str->n_cols; j++){
            AtA_h[j*mlff_str->n_cols+i] = AtA[j*mlff_str->n_cols+i];//+ (sigma_v0*sigma_v0)/(sigma_w0*sigma_w0);
        }
    }

    int ipiv[mlff_str->n_cols], info;

    info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', mlff_str->n_cols, AtA_h, mlff_str->n_cols, lambda_0); // use LAPACKE_dsyevd
    double error_w = 1, error_v=1, error_tol_v = 0.001, error_tol_w = 1, gamma, sigma_v, sigma_w;
    int count=0;

    double regul_min = get_regularization_min(AtA, mlff_str->n_cols);
    printf("Regularization minimum: 1e%f\n", log10(regul_min));
    printf("Loop for Bayesian regression hyperparameters optimization: \n");
    while( (error_w > error_tol_w || error_v > error_tol_v) && count <100) {
    	gamma = 0;
    	for (int k = 0; k < mlff_str->n_cols; k++){
    		lambda[k] = lambda_0[k]/(sigma_v0*sigma_v0);
    		if (lambda[k] < 0.0000000001){
    			gamma += lambda[k]/(lambda[k] + 1.0/(sigma_w0*sigma_w0));
    		}
    	}

      for (int i = 0; i <mlff_str->n_cols*mlff_str->n_cols; i++ )
        AtA_h[i] = AtA[i];

      for (int i =0; i <mlff_str->n_cols; i++)
        AtA_h[i*mlff_str->n_cols+i] +=  (sigma_v0*sigma_v0)/(sigma_w0*sigma_w0);
        
        info = LAPACKE_dgesv( LAPACK_COL_MAJOR, mlff_str->n_cols, 1, AtA_h, mlff_str->n_cols, &ipiv[0],
                     Atb_h, mlff_str->n_cols );
        for (int i=0; i <mlff_str->n_cols; i++){
          weights[i] = Atb_h[i];
        }
            
        for (int i =0; i <mlff_str->n_cols; i++)
            Atb_h[i] = Atb[i];

        double norm_w2 = dotProduct(weights, weights, mlff_str->n_cols);

        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
                mlff_str->n_rows, 1, mlff_str->n_cols, 1.0, A, mlff_str->n_rows, weights, mlff_str->n_cols, 0.0, Ax, mlff_str->n_rows);
        
        for (int i = 0; i < mlff_str->n_rows; i++){
            e[i] = Ax[i]-b[i];
        }

        double norm_er2 = dotProduct(e, e, mlff_str->n_rows);


      // info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, mlff_str->n_cols, mlff_str->n_cols, AtA_h, mlff_str->n_cols, &ipiv[0]);
      // info = LAPACKE_dgecon(LAPACK_COL_MAJOR, 'I', mlff_str->n_cols, AtA_h, mlff_str->n_cols, &ipiv[0]);
      sigma_w = sqrt(norm_w2/gamma);
      sigma_v = sqrt(norm_er2/(mlff_str->n_rows - gamma));
      printf("sigma_v: %f and sigma_w: %f before check for cond\n",sigma_v,sigma_w);
      if ((sigma_v*sigma_v)/(sigma_w*sigma_w) < regul_min){
        sigma_w = sqrt(sigma_v*sigma_v/regul_min);
        printf("regularization: %f is too small as computed from the loop\n", (sigma_v*sigma_v/(sigma_w*sigma_w)));
      }

    	error_w = fabs(sigma_w-sigma_w0);
    	error_v = fabs(sigma_v-sigma_v0);
    	sigma_w0 = sigma_w;
    	sigma_v0 = sigma_v;
        printf("sigma_v: %f, sigma_w: %f\n",sigma_v,sigma_w);
      count++;
   }


   mlff_str->sigma_w = sigma_w;
   mlff_str->sigma_v = sigma_v;

   for (int i=0; i < mlff_str->n_cols; i++){
    mlff_str->weights[i] = weights[i];
   }

   if (error_w<0.001 &&  error_v <0.001){
        printf("sigma_v and sigma_w converged in %d iteration\n",count);
   }
    else{
        printf("sigma_v and sigma_w did not converge and error was %10.9f, %10.9f\n",error_v, error_w);
   }


   free(AtA_h);
   free(Atb_h);
   free(weights);
   free(Ax);
   free(e);
   free(lambda_0);
   free(lambda);
}

/*
mlff_predict predcits E, F, stress and error estimates for a new structure

[Input]
1. K_predict: Design matrix for the new structure
2. mlff_str: MLFF_Obj structure
7. natoms: number of atoms
[Output]
1. Etotal: predicted energy of the new struture
2. F: force of the predicted struture
3. stress: stress of the predicted struture
4. error_bayesian: Bayesian error in the predictions
*/

void mlff_predict(double *K_predict, MLFF_Obj *mlff_str, double *E,  double* F,
                 double* stress, double* error_bayesian, int natoms ){
	// K_predict is arranged column major

	int rows = 3*natoms + 1 + 6, cols = mlff_str->n_cols;
	double *b_predict, *cov_predict1, *error_predict;
	b_predict = (double *) malloc(rows *sizeof(double));
	cov_predict1 = (double *) malloc(rows *cols* sizeof(double));
	error_predict = (double *) malloc(rows *rows* sizeof(double));


	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
    			rows, 1, cols, 1.0, K_predict, rows, mlff_str->weights, cols, 0.0, b_predict, rows);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 
    			cols, rows, cols, 1.0, mlff_str->cov_train, cols, K_predict, rows, 0.0, cov_predict1, cols);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
    			rows, rows, cols, 1.0, K_predict, rows, cov_predict1, cols, 0.0, error_predict, rows);
 
	
	for (int i=0; i<rows; i++){
        int quot = i%(7+3*mlff_str->natom);
        if (quot==0){
            error_bayesian[i] = (mlff_str->std_E)*sqrt(error_predict[i*rows+i]);
        } else if (quot>0 && quot <= 3*mlff_str->natom){
            error_bayesian[i] = (1.0/mlff_str->relative_scale_F)*(mlff_str->std_F)*sqrt(error_predict[i*rows+i]);
        }else{
            error_bayesian[i] = (1.0/mlff_str->relative_scale_stress[quot-3*mlff_str->natom-1])*(mlff_str->std_stress[quot-3*mlff_str->natom-1])*sqrt(error_predict[i*rows+i]);
	      }
  }

	free(error_predict);
	free(cov_predict1);

	E[0] = b_predict[0] * mlff_str->std_E + mlff_str->mu_E;

	for (int i=0; i<natoms; i++){
		for (int j=0; j<3; j++){
			F[i*3+j] = b_predict[i*3+j+1] * mlff_str->std_F * (1.0/mlff_str->relative_scale_F);
		}
	}

	for (int i=0; i<6; i++){
		stress[i] = b_predict[1+3*natoms+i]*mlff_str->std_stress[i]*(1.0/mlff_str->relative_scale_stress[i]);
	}
  free(b_predict);
}


double get_regularization_min(double *A, int size){
  double *a_copy;
  int ipiv[size];
  double rcond;
  double reg_final;
  a_copy = (double*)malloc(size*size*sizeof(double));
  
  double reg_temp[9] = {1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7};

  for (int k=0; k < 9; k++){
    for (int i=0; i < size*size; i++){
      a_copy[i] =A[i];
    }
    for (int i=0; i < size; i++){
      a_copy[i+size*i] += reg_temp[k];
    }
    double norm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', size, size, a_copy, size);
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, size, size, a_copy, size, &ipiv[0]);
    LAPACKE_dgecon(LAPACK_COL_MAJOR, 'I', size, a_copy, size, norm, &rcond);
    printf("Regularization: 1e-%f, Condition number reciprocal: 1e%f\n",log10(reg_temp[k]), log10(rcond) );
    if (rcond > 1e-16){
      reg_final = reg_temp[k];
      break;
    }
  }
  free(a_copy);
  return reg_final;
}
