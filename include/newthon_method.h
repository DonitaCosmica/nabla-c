#pragma once
#include "matrix.h"

typedef double (*func_t)(double *);

size_t newton_method_nd(size_t iterations, size_t dimensions, matrix_t *table, func_t target_func);