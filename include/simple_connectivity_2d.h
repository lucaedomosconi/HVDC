#ifndef SIMPLE_CONNECTIVITY_H
#define SIMPLE_CONNECTIVITY_H

#include <p4est.h>

constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] = 
  {0., 0., 1.e-3, 0., 1.e-3, 1.e-3, 0., 1.e-3};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] = 
  {1, 2, 3, 4, 1};
#endif