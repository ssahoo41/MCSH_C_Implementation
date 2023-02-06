#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "calculate_gmpordernorm.h"

const int NUM_IMPLEMENTED_TYPE = 73;
const int IMPLEMENTED_MCSH_TYPE[73][2] = {
    {0, 0},
    {1, 0},
    {2, 0},
    {3, 0},
    {4, 0},
    {5, 0},
    {6, 0},
    {7, 0},
    {8, 0},
    {9, 0},
    {0, 1},
    {1, 1},
    {2, 1},
    {3, 1},
    {4, 1},
    {5, 1},
    {6, 1},
    {7, 1},
    {8, 1},
    {9, 1},
    {0, 1},
    {1, 1},
    {2, 1},
    {2, 2},
    {3, 1},
    {3, 2},
    {3, 3},
    {4, 1},
    {4, 2},
    {4, 3},
    {4, 4},
    {5, 1},
    {5, 2},
    {5, 3},
    {5, 4},
    {5, 5},
    {6, 1},
    {6, 2},
    {6, 3},
    {6, 4},
    {6, 5},
    {6, 6},
    {6, 7},
    {7, 1},
    {7, 2},
    {7, 3},
    {7, 4},
    {7, 5},
    {7, 6},
    {7, 7},
    {7, 8},
    {8, 1},
    {8, 2},
    {8, 3},
    {8, 4},
    {8, 5},
    {8, 6},
    {8, 7},
    {8, 8},
    {8, 9},
    {8, 10},
    {9, 1},
    {9, 2},
    {9, 3},
    {9, 4},
    {9, 5},
    {9, 6},
    {9, 7},
    {9, 8},
    {9, 9},
    {9, 10},
    {9, 11},
    {9, 12}
};

int calculate_gmpordernorm_noderiv(double** cell, double** cart, double** scale, int* pbc_bools,
                                        int* atom_i, int natoms, int* cal_atoms, int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];


    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    // int **bin_i = (int **) malloc(natoms*sizeof(int*))
    int **bin_i = (int*) malloc(natoms*sizeof(int*));
    for (int i=0; i<natoms; i++) {
        bin_i[i] = (int*) malloc(sizeof(int)*4);
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    total_bins = 1;

    // calculate the inverse matrix of cell and the distance between cell plane
    cross[0][0] = cell[1][1]*cell[2][2] - cell[1][2]*cell[2][1];
    cross[0][1] = cell[1][2]*cell[2][0] - cell[1][0]*cell[2][2];
    cross[0][2] = cell[1][0]*cell[2][1] - cell[1][1]*cell[2][0];
    cross[1][0] = cell[2][1]*cell[0][2] - cell[2][2]*cell[0][1];
    cross[1][1] = cell[2][2]*cell[0][0] - cell[2][0]*cell[0][2];
    cross[1][2] = cell[2][0]*cell[0][1] - cell[2][1]*cell[0][0];
    cross[2][0] = cell[0][1]*cell[1][2] - cell[0][2]*cell[1][1];
    cross[2][1] = cell[0][2]*cell[1][0] - cell[0][0]*cell[1][2];
    cross[2][2] = cell[0][0]*cell[1][1] - cell[0][1]*cell[1][0];

    vol = cross[0][0]*cell[0][0] + cross[0][1]*cell[0][1] + cross[0][2]*cell[0][2];
    inv[0][0] = cross[0][0]/vol;
    inv[0][1] = cross[1][0]/vol;
    inv[0][2] = cross[2][0]/vol;
    inv[1][0] = cross[0][1]/vol;
    inv[1][1] = cross[1][1]/vol;
    inv[1][2] = cross[2][1]/vol;
    inv[2][0] = cross[0][2]/vol;
    inv[2][1] = cross[1][2]/vol;
    inv[2][2] = cross[2][2]/vol;

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = (int*) malloc(total_bins*sizeof(int));
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
        }
        bin_i[i][3] = bin_i[i][0] + nbins[0]*bin_i[i][1] + nbins[0]*nbins[1]*bin_i[i][2];
        atoms_bin[bin_i[i][3]]++;
    }

    max_atoms_bin = 0;
    for (int i=0; i < total_bins; ++i) {
        if (atoms_bin[i] > max_atoms_bin)
            max_atoms_bin = atoms_bin[i];
    }

    free(atoms_bin);

    // # of bins in each direction
    neigh_check_bins = 1;
    for (int i=0; i < 3; ++i) {
        bin_range[i] = ceil(cutoff * nbins[i] / plane_d[i]);
        neigh_check_bins *= 2*bin_range[i];
    }

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        int i=cal_atoms[ii];
        // calculate neighbor atoms
        double *nei_list_d = (double*) malloc(max_atoms_bin * 4 * neigh_check_bins*sizeof(double));
        int    *nei_list_i = (int*) malloc(max_atoms_bin * 2 * neigh_check_bins*sizeof(int));
        nneigh = 0;

        for (int j=0; j < 3; ++j) {
            max_bin[j] = bin_i[i][j] + bin_range[j];
            min_bin[j] = bin_i[i][j] - bin_range[j];
        }

        for (int dx=min_bin[0]; dx < max_bin[0]+1; ++dx) {
            for (int dy=min_bin[1]; dy < max_bin[1]+1; ++dy) {
                for (int dz=min_bin[2]; dz < max_bin[2]+1; ++dz) {
                    pbc_bin[0] = (dx%nbins[0] + nbins[0]) % nbins[0];
                    pbc_bin[1] = (dy%nbins[1] + nbins[1]) % nbins[1];
                    pbc_bin[2] = (dz%nbins[2] + nbins[2]) % nbins[2];
                    cell_shift[0] = (dx-pbc_bin[0]) / nbins[0];
                    cell_shift[1] = (dy-pbc_bin[1]) / nbins[1];
                    cell_shift[2] = (dz-pbc_bin[2]) / nbins[2];

                    bin_num = pbc_bin[0] + nbins[0]*pbc_bin[1] + nbins[0]*nbins[1]*pbc_bin[2];

                    for (int j=0; j < natoms; ++j) {
                        if (bin_i[j][3] != bin_num)
                            continue;

                        // same atom
                        // if (!(cell_shift[0] || cell_shift[1] || cell_shift[2]) && (i == j))
                        //     continue;

                        // take care of pbc
                        if (!pbc_bools[0] && cell_shift[0] != 0)
                            continue;

                        if (!pbc_bools[1] && cell_shift[1] != 0)
                            continue;

                        if (!pbc_bools[2] && cell_shift[2] != 0)
                            continue;

                        for (int a=0; a < 3; ++a) {
                            total_shift[a] = cell_shift[0]*cell[0][a] + cell_shift[1]*cell[1][a] + cell_shift[2]*cell[2][a]
                                             + cart[j][a] - cart[i][a];
                        }

                        // tmp = sqrt(total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2]);
                        tmp_r2 = total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2];
                        //printf("tmp_r2: %f\n",tmp_r2);
                        // if (tmp < cutoff) {
                        if (tmp_r2 < cutoff_sqr) {
                            for (int a=0; a < 3; ++a)
                                nei_list_d[nneigh*4 + a] = total_shift[a];
                            nei_list_d[nneigh*4 + 3] = tmp_r2;
                            nei_list_i[nneigh*2]    = atom_i[j];
                            nei_list_i[nneigh*2 + 1] = j;
                            nneigh++;
                        }
                    }
                }
            }
        }
        //printf("NNeighbors: %d\n",nneigh);
        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], inv_rs = params_d[m][5];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double sum_square = 0.0;
            // std::cout << "------------" << std::endl;
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                GMPFunctionNoderiv mcsh_function = get_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                // std::cout << "\t" << group_index  << "\t"<< mcsh_type << "\t" << group_coefficient<< std::endl;

                
                if (mcsh_type == 1){
                    double sum_desc = 0.0;
                    double m_desc[1];
                    //printf("atom %d\n",i);
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        //printf("Neigh: %d\n",neigh_atom_element_index);
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, m_desc);
                            sum_desc += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_desc * sum_desc;
                    printf("sum_miu: %f\n",sum_desc);
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0];
                            sum_miu2 += miu[1];
                            sum_miu3 += miu[2];
                            //printf("sum_miu1: %f sum_miu2: %f sum_miu3: %f\n",sum_miu1,sum_miu2,sum_miu3);
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                    //printf("sum_miu1: %f sum_miu2: %f sum_miu3: %f\n",sum_miu1,sum_miu2,sum_miu3);
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0];
                            sum_miu2 += miu[1];
                            sum_miu3 += miu[2];
                            sum_miu4 += miu[3];
                            sum_miu5 += miu[4];
                            sum_miu6 += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    printf("sum_miu1: %f sum_miu2: %f sum_miu3: %f sum_miu4: %f sum_miu5: %f sum_miu6: %f\n",sum_miu1,sum_miu2,sum_miu3,sum_miu4,sum_miu5,sum_miu6);
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

        }
        free(nei_list_d);
        free(nei_list_i);
    }
    /*
    for (int i=0; i<natoms; i++) {
        free(bin_i[i]);
    }
    */
    return 0;
}

