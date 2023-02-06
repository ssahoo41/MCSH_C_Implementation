#ifndef CALCULATE_GMPORDERNORM_H
#define CALCULATE_GMPORDERNORM_H
#include <math.h>
#include "surface_harmonics.h"
#include "solid_harmonics.h"
#include <stdbool.h>

/*int calculate_gmpordernorm(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**, double**);
*/
int calculate_gmpordernorm_noderiv(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**);
/*
int calculate_solid_gmpordernorm(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**, double**);

int calculate_solid_gmpordernorm_noderiv(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**); 
*/
 // Change this when you implement new type!

#endif