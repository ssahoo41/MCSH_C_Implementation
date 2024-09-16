#include "local_type_tools.h"
//#include "mlff_types.h"
#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))


/**
 * @brief Initialize the dynamic array, allocate initial
 *        memory and set size.
 *
 */
void init_dyarray(dyArray *a)
{
    assert(a != NULL);
    int initsize = INIT_CAPACITY;
    a->array = malloc(initsize * sizeof(*a->array));
    assert(a->array != NULL);
    a->len = 0;
    a->capacity = initsize;
}

/**
 * @brief Realloc the dynamic array to the given new capacity.
 *
 *        Note that if the array is extended, all the previous data
 *        are still preserved. If the array is shrinked, all the
 *        previous data up to the new capacity is preserved.
 */
void realloc_dyarray(dyArray *a, int new_capacity)
{
    assert(a != NULL);
    value_type *new_arr = realloc(a->array, new_capacity * sizeof(*a->array));
    assert(new_arr != NULL);
    a->array = new_arr;
    a->capacity = new_capacity;
}

/**
 * @brief Double the capacity of the dynamic array.
 *
 *        Note that if the array is extended, all the previous data
 *        are still preserved. If the array is shrinked, all the
 *        previous data up to the new capacity is preserved.
 */
void dyarray_double_capacity(dyArray *a) {
    assert(a != NULL);
    int new_capacity = a->capacity ? a->capacity << 1 : INIT_CAPACITY;
    new_capacity = max(new_capacity, INIT_CAPACITY);
    realloc_dyarray(a, new_capacity);
}

/**
 * @brief Half the capacity of the dynamic array.
 *
 *        Note that if the array is extended, all the previous data
 *        are still preserved. If the array is shrinked, all the
 *        previous data up to the new capacity is preserved.
 */
void dyarray_half_capacity(dyArray *a) {
    assert(a != NULL);
    int new_capacity = a->capacity >> 1;
    new_capacity = max(new_capacity, INIT_CAPACITY);
    realloc_dyarray(a, new_capacity);
}

/**
 * @brief Append an element to the dynamic array.
 *
 */
void append_dyarray(dyArray *a, value_type element)
{
    assert(a != NULL);

    // double the size of memory allocated if it's len up
    if (a->len == a->capacity) {
        dyarray_double_capacity(a);
    }

    // append the element to array
    a->array[a->len] = element;
    a->len++;
}

/**
 * @brief Pop the last element from the dynamic array.
 *
 */
value_type pop_dyarray(dyArray *a)
{
    assert(a != NULL);
    if (a->len < 1) {
        printf("Error: pop_dyarray target is empty!\n");
        exit(1);
    }

    a->len--; // reduce len by 1

    if (4 * a->len < a->capacity) {
        dyarray_half_capacity(a);
    }

    return a->array[a->len];
}

/**
 * @brief Clear the dynamic array.
 *
 *        This function does not destroy the array, it simply
 *        resets the lenght of the dynamic array to 0, and resets
 *        the capacity.
 */
void clear_dyarray(dyArray *a) {
    assert(a != NULL);
    int initsize = INIT_CAPACITY;
    realloc_dyarray(a, initsize);
    a->len = 0;
}

/**
 * @brief Delete the dynamic array.
 *
 */
void delete_dyarray(dyArray *a)
{
    assert(a != NULL);
    free(a->array);
    a->array = NULL;
    a->len = a->capacity = 0;
}


/**
 * @brief Print the dynamic array.
 *
 */
void print_dyarray(const dyArray *a) {
    printf("([");
    for (int i = 0; i < a->len; i++) {
        if (i > 0) printf(" ");
        printf("%d", a->array[i]);
    }
    printf("], len = %ld, capacity = %ld)\n",a->len,a->capacity);
}

// if array is too long, only show the first 5 and last 5
void show_dyarray(const dyArray *a) {
    if (a->len <= 10) {
        print_dyarray(a);
        return;
    }

    printf("([");
    for (int i = 0; i < 5; i++) {
        if (i > 0) printf(" ");
        printf("%d", a->array[i]);
    }
    printf(" ...");
    for (int i = a->len-5; i < a->len; i++) {
        if (i > 0) printf(" ");
        printf("%d", a->array[i]);
    }
    printf("], len = %ld, capacity = %ld)\n",a->len,a->capacity);
}