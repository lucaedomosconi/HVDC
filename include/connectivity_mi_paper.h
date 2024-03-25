#ifndef MI_PAPER_CONNECTIVITY_H
#define MI_PAPER_CONNECTIVITY_H
#include <p8est.h>

constexpr p4est_topidx_t simple_conn_num_vertices = 8;
constexpr p4est_topidx_t simple_conn_num_trees = 1;

const double simple_conn_p[simple_conn_num_vertices*3] = 
  {-8.0e-4, -5.e-4, -1.5e-4,
  8.0e-4, -5.e-4, -1.5e-4,
  -8.0e-4, 5.e-4, -1.5e-4,
  8.0e-4, 5.e-4, -1.5e-4,
  -8.0e-4, -5.e-4, 1.5e-4,
  8.0e-4, -5.e-4, 1.5e-4,
  -8.0e-4, 5.e-4, 1.5e-4,
  8.0e-4, 5.e-4, 1.5e-4};

const p4est_topidx_t simple_conn_t[simple_conn_num_trees*9] = 
  {1, 2, 3, 4, 5, 6, 7, 8, 1};
#endif