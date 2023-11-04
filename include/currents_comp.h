#ifndef CURRENTS_COMP_H
#define CURRENTS_COMP_H

#include <array>
#include <vector>
#include <set>
#include <mpi.h>

#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <quad_operators_3d.h>
#include <simple_connectivity_3d.h>

using q1_vector = q1_vec<distributed_vector>;


class conduction_current {
private:
// Matrices
  distributed_sparse_matrix phi_stiffness;
  distributed_sparse_matrix rho_mass;
// Boundary nodes
  std::array<std::set<size_t>,6> border_nodes;
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
  
  double operator() (const q1_vector & sol1,
                     const q1_vector & sol0,
                     double dt, int side) {
        // Preparing vectors
    drho_dt.get_owned_data().assign(ln_nodes, 0.0);
    phi.get_owned_data().assign(ln_nodes, 0.0);
    for (int ii = 0; ii < ln_nodes; ++ii) {
      drho_dt[drho_dt.get_range_start()+ii] = (sol1.get_owned_data()[ii*N_vars+ord_rho]-sol0.get_owned_data()[ii*N_vars+ord_rho]) / dt;
      phi[phi.get_range_start()+ii] = sol1.get_owned_data()[ii*N_vars+ord_phi];
    }
    // Matrix * vector
    q1_vector I_term1 = phi_stiffness * phi;
    q1_vector I_term2 = rho_mass * drho_dt;

    // Storing nodes on boundary "side" if not already stored
    if (!border_built[side]) {
        border_built[side] = true;
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
                        quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
          for (int ii = 0; ii < 8; ii++)
            if (quadrant->e(ii) == side)
              border_nodes[side].insert(quadrant->gt (ii));
    }

    // Summing terms of I_term1 and I_term2 on nodes of boundary "side"
    double ret_value = 0.0;
    for (auto it = border_nodes[side].cbegin(); it != border_nodes[side].cend(); ++it)
      ret_value += I_term1[*it] + I_term2[*it];
    
    // Summing with border parts owned by other ranks
    MPI_Allreduce(MPI_IN_PLACE, &ret_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return ret_value;
  }
};

#endif