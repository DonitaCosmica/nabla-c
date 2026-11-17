#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "matrix.h"

double matrix_determinant(const matrix_t *matrix) {
  if (matrix == NULL || matrix->rows == 0 || matrix->rows != matrix->cols) {
    return 0.0;
  }

  size_t n = matrix->rows;

  if (n == 1) {
    return matrix_get(matrix, 0, 0);
  }

  if (n == 2) {
    return (matrix_get(matrix, 0, 0) * matrix_get(matrix, 1, 1)) -
           (matrix_get(matrix, 1, 0) * matrix_get(matrix, 0, 1));
  }

  if (n > (SIZE_MAX - sizeof(matrix_t)) / (n * sizeof(double))) {
    return 0.0;
  }

  size_t elems = n * n;
  matrix_t *work = malloc(sizeof(*work) + elems * sizeof(double));
  if (work == NULL) {
    return 0.0;
  }

  work->rows = n;
  work->cols = n;
  for (size_t i = 0; i < elems; ++i) {
    work->data[i] = matrix->data[i];
  }

  double det = 1.0;
  for (size_t k = 0; k < n; ++k) {
    size_t pivot_row = k;
    double max_val = fabs(*matrix_at(work, k, k));

    for (size_t row = k + 1; row < n; ++row) {
      double val = fabs(*matrix_at(work, row, k));
      if (val > max_val) {
        max_val = val;
        pivot_row = row;
      }
    }

    if (max_val <= eps) {
      free(work);
      return 0.0;
    }

    if (pivot_row != k) {
      for (size_t j = k; j < n; ++j) {
        double *a = matrix_at(work, k, j);
        double *b = matrix_at(work, pivot_row, j);
        double tmp = *a;
        *a = *b;
        *b = tmp;
      }

      det = -det;
    }

    double pivot = *matrix_at(work, k, k);
    det *= pivot;

    for (size_t i = k + 1; i < n; ++i) {
      double factor = *matrix_at(work, i, k) / pivot;
      *matrix_at(work, i, k) = 0.0;

      for (size_t j = k + 1; j < n; ++j) {
        *matrix_at(work, i, j) -= factor * *matrix_at(work, k, j);
      }
    }
  }

  free(work);
  return det;
}

bool adjugate_matrix(matrix_t *matrix) {
  if (matrix == NULL || matrix->rows == 0 || matrix->rows != matrix->cols) {
    return false;
  }

  size_t n = matrix->rows;

  if (n == 1) {
    *matrix_at(matrix, 0, 0) = 1.0;
    return true;
  }

  if (n == 2) {
    double a = *matrix_at(matrix, 0, 0);
    double b = *matrix_at(matrix, 0, 1);
    double c = *matrix_at(matrix, 1, 0);
    double d = *matrix_at(matrix, 1, 1);

    *matrix_at(matrix, 0, 0) = d;
    *matrix_at(matrix, 0, 1) = -b;
    *matrix_at(matrix, 1, 0) = -c;
    *matrix_at(matrix, 1, 1) = a;
    return true;
  }

  if (n > (SIZE_MAX - sizeof(matrix_t)) / (n * sizeof(double))) {
    return false;
  }

  size_t elems = n * n;
  size_t sub_n = n - 1;
  size_t sub_elems = sub_n * sub_n;

  matrix_t *adj = malloc(sizeof(*adj) + elems * sizeof(double));
  matrix_t *sub = malloc(sizeof(*sub) + sub_elems * sizeof(double));
  if (adj == NULL || sub == NULL) {
    free(adj);
    free(sub);
    return false;
  }

  adj->rows = n;
  adj->cols = n;
  sub->rows = sub_n;
  sub->cols = sub_n;

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      size_t sub_r = 0;
      for (size_t r = 0; r < n; ++r) {
        if (r == i) {
          continue;
        }

        size_t sub_c = 0;
        for (size_t c = 0; c < n; ++c) {
          if (c == j) {
            continue;
          }

          *matrix_at(sub, sub_r, sub_c) = *matrix_at(matrix, r, c);
          ++sub_c;
        }
        ++sub_r;
      }

      double minor_det = matrix_determinant(sub);
      double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;

      *matrix_at(adj, j, i) = sign * minor_det;
    }
  }

  for (size_t k = 0; k < elems; ++k) {
    matrix->data[k] = adj->data[k];
  }

  free(sub);
  free(adj);
  return true;
}