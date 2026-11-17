#pragma once

#include <stdbool.h>
#include <stddef.h>

static const double eps = 1.0e-12;

typedef struct {
  size_t rows;
  size_t cols;
  double data[];
} matrix_t;

static inline double *matrix_at(matrix_t *matrix, size_t row, size_t col) {
  return &matrix->data[row * matrix->cols + col];
}

static inline double matrix_get(const matrix_t *matrix, size_t row, size_t col) {
  return matrix->data[row * matrix->cols + col];
}

double matrix_determinant(const matrix_t *matrix);
bool adjugate_matrix(matrix_t *matrix);