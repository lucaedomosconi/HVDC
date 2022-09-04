#include <tmesh_3d.h>
#include <simple_connectivity_3d_thin.h>

constexpr double b_x{0.0005}, b_y{0.0}, b_z{0.0005}, b_r{0.0001};

constexpr int NUM_REFINEMENTS = 4;
constexpr int maxlevel_ball = 8;

constexpr double DELTAT = 2.0;
constexpr double T = 200;
constexpr double tau = 5.0; 
constexpr bool save_sol = true;
constexpr double tol = 2e-6;

// Problem parameters
constexpr double epsilon_0 = 8.8542e-12;
constexpr double epsilon_r_1 = 2.0;         // permittivity
constexpr double sigma_1 = 3.21e-14;           // conducivity coeff
constexpr double sigma_2 = 1e-12;

constexpr size_t N_rhos = 4;
std::vector<size_t> rho_idx;
std::vector<std::vector<double>> points{{5e-4,5e-4,0.0},{5e-4,0.0,4e-4},{5e-4,0.0,6e-4},{5e-4,5e-4,1e-3}};
std::vector<std::vector<double>> tols{{1e-4,1e-4,tol},{tol,tol,tol},{tol,tol,tol},{1e-4,1e-4,tol}};
bool extra_refinement = true;

static int
uniform_refinement (tmesh_3d::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
refinement (tmesh_3d::quadrant_iterator quadrant)
{
  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double xcoord, ycoord, zcoord, dist = .0;
  int retval = 0;
  for (int ii = 0; ii < 8; ++ii)
    {
      xcoord = quadrant->p(0, ii);
      ycoord = quadrant->p(1, ii);
      zcoord = quadrant->p(2, ii);
      
      dist = std::sqrt (std::pow (xcoord - b_x, 2) +
                        std::pow (ycoord - b_y, 2) +
                        std::pow (zcoord - b_z, 2));


      if (dist < 1.1 * b_r)
        {
          retval = maxlevel_ball - currentlevel;
          break;
        }
    }

  if (currentlevel >= maxlevel_ball)
    retval = 0;
      
  return (retval);
}

static int
coarsening (tmesh_3d::quadrant_iterator quadrant)
{
  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double xcoord, ycoord, zcoord, dist = .0;
  int retval = currentlevel - NUM_REFINEMENTS;
  for (int ii = 0; ii < 8; ++ii)
    {
      xcoord = quadrant->p(0, ii);
      ycoord = quadrant->p(1, ii);
      zcoord = quadrant->p(2, ii);
      
      dist = std::sqrt (std::pow (xcoord - b_x, 2) +
                        std::pow (ycoord - b_y, 2) +
                        std::pow (zcoord - b_z, 2));


      if (dist < 1.1 * b_r && dist > 0.9 * b_r)
        {
          retval = 0;
          break;
        }
    }

  if (currentlevel <= NUM_REFINEMENTS)
    retval = 0;
      
  return (retval);
}

double epsilon_fun(const double & x, const double & y, const double & z)
{
  double dist = std::sqrt (std::pow (x - b_x, 2) +
                           std::pow (y - b_y, 2) +
                           std::pow (z - b_z, 2));


  return dist < b_r ? epsilon_0 : epsilon_0 * epsilon_r_1;
}

double sigma_fun(const double & x, const double & y, const double & z)
{
  double dist = std::sqrt (std::pow (x - b_x, 2) +
                           std::pow (y - b_y, 2) +
                           std::pow (z - b_z, 2));


  return dist < b_r ? sigma_2 * DELTAT : sigma_1 * DELTAT;
}


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
                        std::cout << "Point " << i+1 << ": x= " << quadrant->p(0,ii) << ", y= " << quadrant->p(1,ii) << ", z= " << quadrant->p(2,ii) << std::endl;
            break;
          }
        }

        if (found)
          break;
    }
    if(!found)
        std::cout << "Node " << i+1 << " not found in current rank" << std::endl; 
  }
  return id;
}