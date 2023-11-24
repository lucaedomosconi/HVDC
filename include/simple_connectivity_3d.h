#ifndef SIMPLE_CONNECTIVITY_H
#define SIMPLE_CONNECTIVITY_H
#include <p8est.h>

constexpr p4est_topidx_t simple_conn_num_vertices = 8;
constexpr p4est_topidx_t simple_conn_num_trees = 1;

const double simple_conn_p[simple_conn_num_vertices*3] = 
  {0., 0., 0.,
  1.e-3, 0., 0.,
  0., 1.e-3, 0.,
  1.e-3, 1.e-3, 0.,
  0., 0., 1.e-3,
  1.e-3, 0., 1.e-3,
  0., 1.e-3, 1.e-3,
  1.e-3, 1.e-3, 1.e-3};
/*
  const double simple_conn_p[simple_conn_num_vertices*3] = 
  {0., 0., 0.,
  5.e-3, 0., 0.,
  0., 5.e-3, 0.,
  5.e-3, 5.e-3, 0.,
  0., 0., 1.e-3,
  5.e-3, 0., 1.e-3,
  0., 5.e-3, 1.e-3,
  5.e-3, 5.e-3, 1.e-3};
*/
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*9] = 
  {1, 2, 3, 4, 5, 6, 7, 8, 1};
#endif