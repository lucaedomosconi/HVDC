#ifndef CURRENTS_COMP_H
#define CURRENTS_COMP_H

#include <array>
#include <vector>
#include <set>
#include <mpi.h>

#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <quad_operators_3d.h>
#include <connectivity_paper_layer.h>

using q1_vector = q1_vec<distributed_vector>;


class conduction_current {
private:
// Matrices
  distributed_sparse_matrix phi_stiffness;
  distributed_sparse_matrix rho_mass;
// Boundary nodes
  std::array<std::set<size_t>,6> border_nodes;
  //std::array<std::set<size_t>,6> border_nodes_adjacent;
  std::array<bool,6> border_built;
// Problem data
  int ln_nodes, N_vars, ord_phi, ord_rho;
  tmesh_3d & tmsh;
// support vectors
    q1_vector drho_dt, phi;
public:
  conduction_current (tmesh_3d & tmsh_,
                      const std::vector<double> & sigma,
                      const std::vector<double> & unitary_vec,
                      const q1_vector & null_q1_vec,
                      const q1_vector & unitary_q1_vec,
                      int ln_nodes_,
                      int ord_rho_ = 0, int ord_phi_ = 1, int N_vars_ = 5) :
                          ln_nodes{ln_nodes_}, N_vars{N_vars_}, tmsh(tmsh_),
                          ord_phi(ord_phi_), ord_rho(ord_rho_), drho_dt(ln_nodes_), phi(ln_nodes_) {

    phi_stiffness.set_ranges(ln_nodes);
    rho_mass.set_ranges(ln_nodes);

    bim3a_advection_diffusion (tmsh, sigma, null_q1_vec, phi_stiffness, true);
    bim3a_reaction (tmsh, unitary_vec, unitary_q1_vec, rho_mass);

    phi_stiffness.assemble();
    rho_mass.assemble();

    border_built.fill(false);
  }
  void sigma_update (const std::vector<double> & sigma,
                     const q1_vector & null_q1_vec) {
    phi_stiffness.reset ();
    bim3a_advection_diffusion (tmsh, sigma, null_q1_vec, phi_stiffness, true);
    phi_stiffness.assemble();
  }

  double operator() (const q1_vector & sol1,
                     const q1_vector & sol0,
                     double dt, int side, bool only_dphi = false) {
        // Preparing vectors
    drho_dt.get_owned_data().assign(ln_nodes, 0.0);
    phi.get_owned_data().assign(ln_nodes, 0.0);
    for (int ii = 0; ii < ln_nodes; ++ii) {
      drho_dt[drho_dt.get_range_start()+ii] = only_dphi ? 0 : (sol1.get_owned_data()[ii*N_vars+ord_rho]-sol0.get_owned_data()[ii*N_vars+ord_rho]) / dt;
      phi[phi.get_range_start()+ii] = sol1.get_owned_data()[ii*N_vars+ord_phi];
    }
    // Matrix * vector
    q1_vector I_term1 = phi_stiffness * phi;
    q1_vector I_term2 = rho_mass * drho_dt;

    // Storing nodes on boundary side if not already stored
    if (!border_built[side]) {
      border_built[side] = true;
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
                      quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
        for (int ii = 0; ii < 8; ii++)
//          if (!quadrant->is_hanging(ii)) {
            if (quadrant->e(ii) == side)
              border_nodes[side].insert(quadrant->gt (ii));
//            if (side == 5 && quadrant->p(2,ii) > 1e-3 *12/15)
//              border_nodes_adjacent[side].insert(quadrant->gt(ii));
//          }          
    }

    // Summing terms of I_term1 and I_term2 on nodes of boundary "side"
    double ret_value = 0.0;
    for (auto it = border_nodes[side].cbegin(); it != border_nodes[side].cend(); ++it)
      ret_value += I_term1[*it] + I_term2[*it];
//    for (auto it = border_nodes_adjacent[side].cbegin(); it != border_nodes_adjacent[side].cend(); ++it)
//      ret_value += I_term2[*it];
    
    // Summing with border parts owned by other ranks
    MPI_Allreduce(MPI_IN_PLACE, &ret_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return ret_value;
  }
};

template <int... orders>
class diff_pol_total_charge {
private:
  int ln_nodes;
  static const int N_vars = 2+sizeof...(orders);
  tmesh_3d & tmsh;
  std::vector<q1_vector> p_charge, p_charge_abs;
  const std::vector<double> & unitary_vec;
  double old_p = 0;

  template<int t, int... ts>
  void build_rhs (q1_vector & vec_mass) {
    bim3a_rhs (tmsh, unitary_vec, p_charge[t-2], vec_mass, dof_ordering<N_vars,t>);
    if constexpr (sizeof...(ts))
      build_rhs<ts...>(vec_mass);
  }
  template<int t, int... ts>
  void build_rhs_a (q1_vector & vec_mass) {
    bim3a_rhs (tmsh, unitary_vec, p_charge_abs[t-2], vec_mass, dof_ordering<N_vars,t>);
    if constexpr (sizeof...(ts))
      build_rhs_a<ts...>(vec_mass);
  }

public:
  diff_pol_total_charge (tmesh_3d & tmsh_,
                         const std::vector<double> & unitary_vec_,
                         int ln_nodes_) :
                            ln_nodes{ln_nodes_}, unitary_vec(unitary_vec_),
                            tmsh(tmsh_) {
    for (auto i = 0; i < sizeof...(orders); ++i) {
      p_charge.emplace_back(ln_nodes);
      p_charge_abs.emplace_back(ln_nodes);
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
              quadrant != tmsh.end_quadrant_sweep (); ++quadrant) {
        for (auto ii = 0; ii < 8; ++ii)
          if (!quadrant->is_hanging(ii)) {
            p_charge.back()[quadrant->gt(ii)] = 0.0;
            p_charge_abs.back()[quadrant->gt(ii)] = 0.0;
          }
          else
            for (auto jj = 0; jj < quadrant->num_parents (ii); ++jj) {
              p_charge.back()[quadrant->gparent(jj, ii)] += 0.0;
              p_charge_abs.back()[quadrant->gparent(jj, ii)] += 0.0;
            }
      }
      p_charge.back().assemble(replace_op);
      p_charge_abs.back().assemble(replace_op);
    }
  }

  void operator () (const q1_vector & sol,
                    double dt, double & polcharge, 
                    double & abs_polcharge, double & charge_variation) {
    for (auto i = 0; i < sizeof...(orders); ++i) {
      for (auto j = 0; j < ln_nodes; ++j) {
        p_charge[i].get_owned_data()[j] = sol.get_owned_data()[i+2+j*N_vars];
        p_charge_abs[i].get_owned_data()[j] = std::fabs(sol.get_owned_data()[i+2+j*N_vars]);
      }
      p_charge[i].assemble(replace_op);
      p_charge_abs[i].assemble(replace_op);
    }

    q1_vector p_charge_mass(ln_nodes*N_vars), p_charge_abs_mass(ln_nodes*N_vars);
    p_charge_mass.get_owned_data().assign(ln_nodes*N_vars, 0.0);
    p_charge_abs_mass.get_owned_data().assign(ln_nodes*N_vars, 0.0);
    build_rhs<orders...>(p_charge_mass);
    build_rhs_a<orders...>(p_charge_abs_mass);
    p_charge_mass.assemble();
    p_charge_abs_mass.assemble();
    polcharge = 0; abs_polcharge = 0;
    for (auto i = 0; i < ln_nodes*N_vars; ++i) {
      polcharge += p_charge_mass.get_owned_data()[i];
      abs_polcharge += p_charge_abs_mass.get_owned_data()[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &polcharge, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &abs_polcharge, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    charge_variation = polcharge - old_p;
    old_p = polcharge;
  }
};
#endif