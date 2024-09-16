#ifndef MLFF_TYPES_H
#define MLFF_TYPES_H

typedef int value_type;
/*
* @brief  Data type for dynamic array.
*/


typedef struct NeighList NeighList;
typedef struct GMPObj GMPObj;
typedef struct MLFF_Obj MLFF_Obj;
typedef struct dyArray dyArray;

typedef struct dyArray{
    value_type *array;
    long int len;  // length of the array (used)
    long int capacity; // total capacity of memory available
} dyArray;

typedef struct NeighList
{
  // # of atoms 
  int natom;
  // # of element species
  int nelem;
  // cutoff distance (Bohr)
  double rcut;
  // length of fundamental unit cell - changed this to work with current format
  double cell[3];
  // # of atom for each eleemnt species
  int* natom_elem;
  // # of neighbours for each atoms
  int* Nneighbors;
  // # of neighbour in the fundamental cell for each atom (excludes neighbours in periodic images)
  int* unique_Nneighbors;
  // # of neighbours of a given element speciy for each atoms
  int** Nneighbors_elemWise;
  // # of neigbours of a given element speciy in the fundamental cell for each atom (excludes neighbours in periodic images) 
  int** unique_Nneighbors_elemWise;
  // list of neighbours of all atoms
  dyArray* neighborList;
  // list of neighbours of a given element speciy
  dyArray** neighborList_elemWise;
  // list of unique neighbours of all atoms
  dyArray* unique_neighborList;
  // list of unique neighbous of a given speciy
  dyArray** unique_neighborList_elemWise;
  // list of the image ID in x-direction of neighbours for all atoms
  dyArray* neighborList_imgX;
  // list of the image ID in y-direction of neighbours for all atoms
  dyArray* neighborList_imgY;
  // list of the image ID in z-direction of neighbours for all atoms
  dyArray* neighborList_imgZ;
  // list of the image ID in x-direction of neighbours of a given element speciy for all atoms
  dyArray** neighborList_elemWise_imgX;
  // list of the image ID in y-direction of neighbours of a given element speciy for all atoms
  dyArray** neighborList_elemWise_imgY;
  // list of the image ID in z-direction of neighbours of a given element speciy for all atoms
  dyArray** neighborList_elemWise_imgZ;
  // stores an ID of all the neighbours of the atoms. [IDs are integers ranging from 0 to N-1 whre N is number of unique neighbours of the concerning atom]
  dyArray* localID_neighbours;
  // stores an ID of all the unique neighbours of the atoms. [IDs are integers ranging from 0 to N-1 whre N is number of unique neighbours of the concerning atom]
  dyArray** localID_neighbours_elem;
  // list atom type of neighbours of a given atom
  dyArray* neighborAtmTyp;
} NeighList;


typedef struct GMPObj
{
  // # of atoms 
  int natom;
  // cutoff distance (Angstrom)
  double rcut;
  // dimensions of the fundamental unit cell (Angstrom)
  double cell[3];
  // size of GMP descriptor
  int size_X;
  // atoms for descriptor generation
  int *cal_atoms;
  // number of atoms for descriptor generation
  int cal_num;
  // GMP specific integer params
  int **params_i;
  // GMP specific double params
  double **params_d;
  // number of elements per descriptor
  int nmcsh;
  // parameters for primitive gaussians
  double** atom_gaussian;
  // number of gaussians per descriptor element
  int* ngaussians;
  // heleper tool
  int* element_index_to_order;
  // # of elements
  int nelem;
  // array of number of neighbours for all the atoms
  int* Nneighbors;
  // array of number of unqiue neighbours for all the atoms
  int* unique_Nneighbors;
  // array of number of atoms of different element species 
  int* natom_elem;
  // array of number of neighbours of different species for all the atoms
  int** unique_Nneighbors_elemWise;
  // list of neighbours of different species in the fundamental cell
  dyArray** unique_neighborList_elemWise;
  // list of neighbours of all the atoms
  dyArray* neighborList;
  // list of unique neighbours of all the atoms
  dyArray* unique_neighborList;
  // array of GMP descriptor, X[natom][size_X]
  double **X;
  // array of GMP descriptor derivatives for x-dir force, dX_dX[natom][1+neigh][size_X]
  double ***dX_dX;
  // array of GMP descriptor derivatives for y-dir force, dX_dY[natom][1+neigh][size_X]
  double ***dX_dY;
  // array of GMP descriptor derivatives for z-dir force, dX_dZ[natom][1+neigh][size_X]
  double ***dX_dZ;
} GMPObj;

typedef struct FeatureScaler
{
  // mean of each column corresponding to each combination of radial and angular probe
  double *mean_colWise;
  // stdev of each column corresponding to each combination of radial and angular probe
  double *stdev_colWise;
  // training data flag - indicates if mean and stdev need to be recalculated
  int train;
} FeatureScaler;

#endif 
