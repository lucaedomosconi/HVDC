#include <tmesh.h>
#include <simple_connectivity_2d.h>
#include <nlohmann/json.hpp>
constexpr int NUM_REFINEMENTS = 4;

// constexpr double DELTAT = 0.25;
double T;
double tau;
double tau_p1, tau_p2, tau_p3;
bool save_sol;

// Problem parameters
double epsilon_0;
double epsilon_inf;       // permittivity at infinite frequency
double csi1, csi2, csi3;
double sigma_;            // conducivity coeff

constexpr size_t N_rhos = 1;
std::vector<size_t> rho_idx;
std::vector<std::vector<double>> points{{0.0005,0.0005}};
std::vector<std::vector<double>> tols{{1e-4,1e-4}};

bool extra_refinement = false;


static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
coarsening (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

double epsilon_fun(const double & x, const double & y)
{return epsilon_0 * epsilon_inf;}

double epsilon_inf_fun(const double & x, const double & y)
{return epsilon_inf;}

double csi_1_fun(const double & x, const double & y)
{return csi1;}

double csi_2_fun(const double & x, const double & y)
{return csi2;}

double csi_3_fun(const double & x, const double & y)
{return csi3;}

double sigma_fun(const double & x, const double & y, const double & DT)
{return sigma_ * DT;}

std::vector<size_t> find_idx(tmesh &tmsh,std::vector<std::vector<double>> &points,std::vector<std::vector<double>> &tols, const size_t &N_rhos)
{
  std::vector<size_t> id(N_rhos,0);

  for (size_t i = 0; i < N_rhos; i++) {
    bool found = false;

    for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant){

        for (int ii = 0; ii < 4; ++ii) {

          if (fabs(quadrant->p(0,ii)-points[i][0])<tols[i][0] && fabs(quadrant->p(1,ii)-points[i][1])<tols[i][1]) {
            id[i] = quadrant->t(ii);
            std::cout << quadrant->p(0,ii) << " " << quadrant->p(1,ii) << " " << quadrant->p(2,ii) << std::endl;
            found = true;
            std::cout << "Point " << i+1 << ": x= " << quadrant->p(0,ii) << ", y= " << quadrant->p(1,ii) << std::endl;
            break;
          }
        }

        if (found)
          break;
    }
    if(!found)
        std::cout << "Node " << i << " not found in current rank" << std::endl; 
  }
  return id;
}