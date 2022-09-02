#include <tmesh_3d.h>
#include <simple_connectivity_3d_thin.h>

constexpr int NUM_REFINEMENTS = 4;
constexpr int maxlevel = 6;

constexpr double DELTAT = 1.0;
constexpr double T = 100;
constexpr double tau = 1.0; 
constexpr bool save_sol = true;

// Problem parameters
constexpr double epsilon_0 = 8.8542e-12;
constexpr double epsilon_r_1 = 2.0;         // permittivity bottom half
constexpr double epsilon_r_2 = 4.0;         // permittivity upper half
constexpr double sigma_ = 3.21e-14;           // conducivity coeff

bool extra_refinement = true;

static int
uniform_refinement (tmesh_3d::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
refinement (tmesh_3d::quadrant_iterator quadrant)
{
  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double zcoord;
  int retval = 0;
  for (int ii = 0; ii < 8; ++ii)
    {
      zcoord = quadrant->p(2, ii);
      
      if (fabs(zcoord - 0.0005) < 0.000001)
        {
          retval = maxlevel - currentlevel;
          break;
        }
    }

  if (currentlevel >= maxlevel)
    retval = 0;
      
  return retval;
}

static int
coarsening (tmesh_3d::quadrant_iterator quadrant)
{
  return 0;
}

double epsilon_fun(const double & x, const double & y, const double & z)
{return z < 0.0005 ? epsilon_0 * epsilon_r_1 : epsilon_0 * epsilon_r_2;}

double sigma_fun(const double & x, const double & y, const double & z)
{return sigma_ * DELTAT;}
