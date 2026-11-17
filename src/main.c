#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "newthon_method.h"

// x^4 - 3x^3 + 2
double f(double *points) {
  return points[0] * points[0] * points[0] * points[0] -
         3.0 * points[0] * points[0] * points[0] + 2.0;
}

// x^3 + x^2y^3 - 2y^2
double g(double *points) {
  return points[0] * points[0] * points[0] +
         points[0] * points[0] * points[1] * points[1] * points[1] -
         2.0 * points[1] * points[1];
}

// x^2 + y^2 + z^2 + xy + yz + z^4
double h(double *points) {
  return points[0] * points[0] +
         2.0 * points[1] * points[1] +
         3.0 * points[2] * points[2] +
         points[0] * points[1] +
         points[1] * points[2] +
         points[2] * points[2] * points[2] * points[2];
}

void export_malla_nd(size_t dimensions, size_t completed, const matrix_t *table, func_t selected_func);

int main(void) {
  func_t funcs[] = { f, g, h };
  size_t dimensions = 0;
  size_t iterations = 0;
  size_t max_funcs = sizeof(funcs) / sizeof(funcs[0]);

  printf("Ingrese el numero de dimensiones de su funcion: ");
  if (scanf("%zu", &dimensions) != 1 || dimensions == 0 || dimensions > max_funcs) {
    fprintf(stderr, "Entrada invalida o dimension no soportada (1-%zu)\n", max_funcs);
    return EXIT_FAILURE;
  }

  func_t selected_func = funcs[dimensions - 1];
  if (selected_func == NULL) {
    fprintf(stderr, "No hay una funcion definida para %zu dimensiones\n", dimensions);
    return EXIT_FAILURE;
  }

  printf("Ingrese el numero de iteraciones: ");
  if (scanf("%zu", &iterations) != 1 || iterations == 0) {
    fprintf(stderr, "Entrada invalida\n");
    return EXIT_FAILURE;
  }

  if (dimensions > SIZE_MAX - 1) {
    fprintf(stderr, "Desbordamiento en numero de columnas\n");
    return EXIT_FAILURE;
  }

  size_t cols = dimensions + 2;
  if (iterations > SIZE_MAX / cols) {
    fprintf(stderr, "Desbordamiento en total de elementos de table\n");
    return EXIT_FAILURE;
  }

  size_t table_elems = iterations * cols;
  if (table_elems > (SIZE_MAX - sizeof(matrix_t)) / sizeof(double)) {
    fprintf(stderr, "Desbordamiento de tamanio en malloc para table\n");
    return EXIT_FAILURE;
  }

  matrix_t *table = malloc(sizeof(*table) + table_elems * sizeof(double));
  if (table == NULL) {
    fprintf(stderr, "No se pudo asignar memoria para table\n");
    return EXIT_FAILURE;
  }

  table->rows = iterations;
  table->cols = cols;

  size_t completed = newton_method_nd(iterations, dimensions, table, selected_func);
  if(completed == 0) {
    free(table);
    return EXIT_FAILURE;
  }

  FILE *data = fopen("docs/data.csv", "w");
  if(data == NULL) {
    fprintf(stderr, "Error al abrir archivo\n");
    free(table);
    return EXIT_FAILURE;
  }

  printf("\nIter");
  fprintf(data, "iter,");
  for (size_t d = 0; d < dimensions; ++d) {
    printf(" |       x%-3zu", d + 1);
    fprintf(data, "%zu,", d + 1);
  }

  printf(" | f(");
  for(size_t d = 0; d < dimensions; ++d) {
    printf("x%zu%s", d + 1, (d + 1 < dimensions) ? ", " : "");
  }

  printf(")\n");
  fprintf(data, "f\n");

  for(size_t i = 0; i < completed; ++i) {
    double iter_val = matrix_get(table, i, 0);
    double fn_value = matrix_get(table, i, dimensions + 1);

    printf("%4.0f", iter_val);
    fprintf(data, "%.0f,", iter_val);
    for(size_t d = 0; d < dimensions; ++d) {
      double value = matrix_get(table, i, d + 1);
      printf(" | %10.6f", value);
      fprintf(data, "%.6f,", value);
    }

    printf(" | %10.6f\n", fn_value);
    fprintf(data, "%.6f\n", fn_value);
  }

  if(fclose(data) != 0) {
    fprintf(stderr, "Error al cerrar data.csv");
    free(table);
    return EXIT_FAILURE;
  }

  export_malla_nd(dimensions, completed, table, selected_func);

  free(table);
  return EXIT_SUCCESS;
}

void export_malla_nd(size_t dimensions, size_t completed, const matrix_t *table, func_t selected_func) {
  if (dimensions > 4) {
    printf("Nota: No se genera malla.csv para %zuD por tamano exponencial.\n", dimensions);
    return;
  }

  size_t steps = 100;
  if (dimensions == 1) steps = 250;
  else if (dimensions == 2) steps = 100; // 10,000 puntos
  else if (dimensions == 3) steps = 25;  // 15,625 puntos
  else if (dimensions == 4) steps = 12;  // 20,736 puntos

  double *min_bounds = malloc(dimensions * sizeof(double));
  double *max_bounds = malloc(dimensions * sizeof(double));
  double *point = malloc(dimensions * sizeof(double));
  size_t *indices = calloc(dimensions, sizeof(size_t));

  if (!min_bounds || !max_bounds || !point || !indices) {
    fprintf(stderr, "Error de memoria al generar malla.csv\n");
    free(min_bounds); free(max_bounds); free(point); free(indices);
    return;
  }

  for (size_t d = 0; d < dimensions; ++d) {
    double min_v = matrix_get(table, 0, d + 1);
    double max_v = min_v;

    for (size_t i = 1; i < completed; ++i) {
      double val = matrix_get(table, i, d + 1);
      if (val < min_v) min_v = val;
      if (val > max_v) max_v = val;
    }

    double margin = (max_v - min_v) * 0.35;
    if (margin < 1.0) margin = 1.5;

    min_bounds[d] = min_v - margin;
    max_bounds[d] = max_v + margin;
  }

  FILE *fm = fopen("docs/malla.csv", "w");
  if (!fm) {
    fprintf(stderr, "Error al abrir docs/malla.csv\n");
    free(min_bounds); free(max_bounds); free(point); free(indices);
    return;
  }

  for (size_t d = 0; d < dimensions; ++d) {
    fprintf(fm, "%zu,", d + 1);
  }
  
  fprintf(fm, "f\n");

  while (1) {
    for (size_t d = 0; d < dimensions; ++d) {
      double t = (double)indices[d] / (double)(steps - 1);
      point[d] = min_bounds[d] + t * (max_bounds[d] - min_bounds[d]);
      fprintf(fm, "%.6f,", point[d]);
    }

    double f_val = selected_func(point);
    fprintf(fm, "%.6f\n", f_val);

    size_t dim = dimensions;
    while (dim > 0) {
      dim--;
      indices[dim]++;
      if (indices[dim] < steps) {
        break;
      }
      indices[dim] = 0;
      if (dim == 0) {
        goto done;
      }
    }
  }

done:
  fclose(fm);
  free(min_bounds);
  free(max_bounds);
  free(point);
  free(indices);
}