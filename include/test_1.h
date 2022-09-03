/*!
\file test_1.h
Test case with uniform sigma and epsilon in the domain
*/
#include <tmesh_3d.h>
#include <simple_connectivity_3d_thin.h>

constexpr int NUM_REFINEMENTS = 4;

constexpr double DELTAT = 0.5;
constexpr double T = 5;
constexpr double tau = 1.0; 
constexpr bool save_sol = true;

// Problem parameters
constexpr double epsilon_0 = 8.8542e-12;
constexpr double epsilon_r = 2.0;         // permittivity
constexpr double sigma_ = 3.21e-14;           // conducivity coeff

constexpr size_t N_rhos = 1;
std::vector<size_t> rho_idx;
std::vector<std::vector<double>> points{{0.0005,0.0005,0.0}};
std::vector<std::vector<double>> tols{{1e-4,1e-4,1e-9}};

bool extra_refinement = false;

static int
uniform_refinement (tmesh_3d::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
refinement (tmesh_3d::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
coarsening (tmesh_3d::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

double epsilon_fun(const double & x, const double & y, const double & z)
{return epsilon_0 * epsilon_r;}

double sigma_fun(const double & x, const double & y, const double & z)
{return sigma_ * DELTAT;}

std::vector<size_t> find_idx(tmesh_3d &tmsh,std::vector<std::vector<double>> &points,std::vector<std::vector<double>> &tols, const size_t &N_rhos)
{
  std::vector<size_t> id(N_rhos,0);

  for (size_t i = 0; i < N_rhos; i++) {
    bool found = false;

    for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant){

        for (int ii = 0; ii < 8; ++ii) {

          if (fabs(quadrant->p(0,ii)-points[i][0])<tols[i][0] && fabs(quadrant->p(1,ii)-points[i][1])<tols[i][1] && fabs(quadrant->p(2,ii)-points[i][2])<tols[i][2]) {
            id[i] = quadrant->t(ii);
            found = true;
            break;
          }
        }

        if (found)
          break;
    }
  }
  return id;
}