#ifndef LOCAL_TYPE_TOOLS_H
#define LOCAL_TYPE_TOOLS_H
#define INIT_CAPACITY 4
#include "mlff_types.h"

void init_dyarray(dyArray *a);

/**
 * @brief Realloc the dynamic array to the given new capacity.
 *
 *        Note that if the array is extended, all the previous data
 *        are still preserved. If the array is shrinked, all the
 *        previous data up to the new capacity is preserved.
 */
void realloc_dyarray(dyArray *a, int new_capacity);

/**
 * @brief Double the capacity of the dynamic array.
 *
 *        Note that if the array is extended, all the previous data
 *        are still preserved. If the array is shrinked, all the
 *        previous data up to the new capacity is preserved.
 */
void dyarray_double_capacity(dyArray *a);

/**
 * @brief Half the capacity of the dynamic array.
 *
 *        Note that if the array is extended, all the previous data
 *        are still preserved. If the array is shrinked, all the
 *        previous data up to the new capacity is preserved.
 */
void dyarray_half_capacity(dyArray *a);

/**
 * @brief Append an element to the dynamic array.
 *
 */
void append_dyarray(dyArray *a, value_type element);

/**
 * @brief Pop the last element from the dynamic array.
 *
 */
value_type pop_dyarray(dyArray *a);

/**
 * @brief Clear the dynamic array.
 *
 *        This function does not destroy the array, it simply
 *        resets the lenght of the dynamic array to 0, and resets
 *        the capacity.
 */
void clear_dyarray(dyArray *a);

/**
 * @brief Delete the dynamic array.
 *
 */
void delete_dyarray(dyArray *a);

/**
 * @brief Print the dynamic array.
 *
 */
void print_dyarray(const dyArray *a);

// if array is too long, only show the first 5 and last 5
void show_dyarray(const dyArray *a); 

#endif