#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <newthon_method.h>

static const double h1 = 1.0e-6;
static const double h2 = 1.0e-4;

static double first_derive(double *p, int var, func_t func) {
  double orig = p[var];
  p[var] = orig + h1;

  double f_plus = func(p);
  p[var] = orig - h1;

  double f_minus = func(p);
  p[var] = orig;

  return (f_plus - f_minus) / (2.0 * h1);
}

static double second_derive(double *p, int var1, int var2, func_t func) {
  double orig1 = p[var1];
  double orig2 = p[var2];

  // fxx, fyy, fzz, ... fórmula directa para derivada pura
  if (var1 == var2) {
    double f_center = func(p);
    p[var1] = orig1 + h2;

    double f_plus = func(p);
    p[var1] = orig1 - h2;

    double f_minus = func(p);
    p[var1] = orig1;

    return (f_plus - 2.0 * f_center + f_minus) / (h2 * h2);
  }

  p[var1] = orig1 + h2; p[var2] = orig2 + h2;

  double f_pp = func(p);
  p[var1] = orig1 + h2; p[var2] = orig2 - h2;

  double f_pm = func(p);
  p[var1] = orig1 - h2; p[var2] = orig2 + h2;

  double f_mp = func(p);
  p[var1] = orig1 - h2; p[var2] = orig2 - h2;

  double f_mm = func(p);
  p[var1] = orig1;
  p[var2] = orig2;

  return (f_pp - f_pm - f_mp + f_mm) / (4.0 * h2 * h2);
}

size_t newton_method_nd(size_t iterations, size_t dimensions, matrix_t *table, func_t target_func) {
  if (dimensions == 0 || table == NULL || target_func == NULL) {
    return 0;
  }

  if (dimensions > SIZE_MAX / dimensions) {
    fprintf(stderr, "Desbordamiento al calcular dimensiones de la matriz\n");
    return 0;
  }

  size_t square_matrix = dimensions * dimensions;
  if (dimensions > (SIZE_MAX - sizeof(matrix_t)) / sizeof(double) ||
      square_matrix > (SIZE_MAX - sizeof(matrix_t)) / sizeof(double)) {
    fprintf(stderr, "Desbordamiento de tamanio en malloc para estructuras n-D\n");
    return 0;
  }

  matrix_t *w_t = malloc(sizeof(*w_t) + dimensions * sizeof(double));
  matrix_t *H = malloc(sizeof(*H) + square_matrix * sizeof(double));
  matrix_t *gradient = malloc(sizeof(*gradient) + dimensions * sizeof(double));

  if (w_t == NULL || H == NULL || gradient == NULL) {
    fprintf(stderr, "No se pudo asignar memoria en newton_method_nd\n");
    free(w_t);
    free(H);
    free(gradient);
    return 0;
  }

  w_t->rows = dimensions;
  w_t->cols = 1;

  H->rows = dimensions;
  H->cols = dimensions;

  gradient->rows = dimensions;
  gradient->cols = 1;

  puts("Ingrese los valores iniciales:");
  for (size_t i = 0; i < dimensions; ++i) {
    printf("x%zu: ", i + 1);
    if (scanf("%lf", matrix_at(w_t, i, 0)) != 1) {
      fprintf(stderr, "Entrada invalida en x%zu\n", i + 1);
      free(w_t);
      free(H);
      free(gradient);
      return 0;
    }
  }

  printf("f(");
  for(size_t d = 0; d < dimensions; ++d) {
    if(d != dimensions - 1) printf("x%zu, ", d + 1);
    else printf("x%zu", d + 1);
  }

  printf("): %.6f\n", target_func(w_t->data));

  size_t completed = 0;

  for (size_t iter = 0; iter < iterations; ++iter) {
    for (size_t i = 0; i < dimensions; ++i) {
      *matrix_at(gradient, i, 0) = first_derive(w_t->data, (int)i, target_func);
    }

    for (size_t r = 0; r < dimensions; ++r) {
      for (size_t c = 0; c < dimensions; ++c) {
        *matrix_at(H, r, c) = second_derive(w_t->data, (int)r, (int)c, target_func);
      }
    }

    double det_H = matrix_determinant(H);
    if (fabs(det_H) <= eps) {
      fprintf(stderr, "Hessiano singular o mal condicionado en iteracion %zu (det = %.3e)\n",
              iter + 1, det_H);
      break;
    }

    if(!adjugate_matrix(H)) {
      fprintf(stderr, "Error al calcular la matriz adjunta en iteracion %zu\n", iter + 1);
      break;
    }

    for (size_t r = 0; r < dimensions; ++r) {
      double acc = 0.0;
      for (size_t c = 0; c < dimensions; ++c) {
        acc += (*matrix_at(H, r, c)) * (*matrix_at(gradient, c, 0));
      }

      double dx = -(acc / det_H);
      *matrix_at(w_t, r, 0) += dx;
    }

    *matrix_at(table, completed, 0) = (double)(iter + 1);
    for (size_t d = 0; d < dimensions; ++d) {
      *matrix_at(table, completed, d + 1) = *matrix_at(w_t, d, 0);
    }

    *matrix_at(table, completed, dimensions + 1) = target_func(w_t->data);

    ++completed;
  }

  free(w_t);
  free(H);
  free(gradient);

  return completed;
}