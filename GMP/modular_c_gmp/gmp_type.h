#ifndef GMP_TYPE_H
#define GMP_TYPE_H

typedef struct dyArray {
    value_type *array;
    size_t len;  // length of the array (used)
    size_t capacity; // total capacity of memory available
} dyArray;

//This basically allows you to refer to these types without `struct name` everytime
typedef struct NeighList NeighList;
typedef struct GMPObj GMPObj;

//start by passing in gmp struct instead of pointer to x
typedef struct NeighList
{
  // # of atoms 
  int natom;
  // # of element species
  int nelem;
  // cutoff distance (Bohr)
  double rcut;
  // length of fundamental unit cell
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
  // dimensions of the fundamental unit cell (Bohr)
  double cell[3];
  // # of atoms 
  int natom;
  // cutoff distance (Bohr)
  double rcut;
  
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
  // array of 2-body descriptor, X2[natom][size_X2]
  double **X2;
  // array of 3-body descriptor SOAP, X3[natom][size_X3]
  double ***dX2_dX;
  // array of 2-body descriptor derivatives for y-dir force, dX2_dY[natom][1+neigh][size_X2]
  double ***dX2_dY;
  // array of 2-body descriptor derivatives for z-dir force, dX2_dZ[natom][1+neigh][size_X2]
  double ***dX2_dZ;
} GMPObj;

#endif