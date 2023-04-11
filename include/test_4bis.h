#include <tmesh.h>
#include <simple_connectivity_2d.h>

constexpr int NUM_REFINEMENTS = 5;
constexpr int maxlevel = 6;
constexpr int minlevel = 4;

constexpr double DELTAT = 50.0;
constexpr double T = 5000;
constexpr double tau = 50.0; 
constexpr bool save_sol = true;

// Problem parameters
constexpr double epsilon_0 = 8.8542e-12;
constexpr double epsilon_r_1 = 2.0;         // permittivity oil
constexpr double epsilon_r_2 = 4.0;         // permittivity paper
constexpr double sigma_1 = 3.21e-14;        // conducivity coeff
constexpr double sigma_2 = 1e-12;

constexpr double y_oil = 5e-5;              //thickness oil layer
constexpr double y_paper = 3e-4;            // thickness paper layer and butt gaps side length
constexpr double tol = 1e-5;

constexpr size_t N_rhos = 6;
std::vector<size_t> rho_idx;
std::vector<std::vector<double>> points{{5e-4,1e-3}, {5e-4,y_paper}, {0.0,y_oil+y_paper}, {y_paper/2,5e-4,y_oil+y_paper},{y_paper,y_oil+y_paper},{y_paper,y_oil+y_paper}};
std::vector<std::vector<double>> tols{{1e-4,tol},{1e-4,tol},{tol,tol},{1e-4,tol},{tol,tol},{tol,tol}};
bool extra_refinement = true;

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
refinement (tmesh::quadrant_iterator quadrant)
{
  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double ycoord;
  int retval = 0;
  for (int ii = 0; ii < 4; ++ii)
    {
      ycoord = quadrant->p(1, ii);
      
      if (ycoord > y_paper - tol || ycoord < 2*y_paper+2*y_oil+tol)
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
coarsening (tmesh::quadrant_iterator quadrant)
{
  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double xcoord,ycoord;
  int retval = currentlevel - minlevel;
  for (int ii = 0; ii < 4; ++ii)
    { 
      xcoord = quadrant->p(0, ii);
      ycoord = quadrant->p(1, ii);    

      if (fabs(ycoord - y_paper) < tol || fabs(ycoord - y_paper-y_oil) < tol || fabs(ycoord - 2*y_paper-y_oil) < tol || fabs(ycoord - 2*y_paper-2*y_oil) < tol || 
         (xcoord<y_paper && fabs(ycoord-5e-4+y_paper/2)<tol && ycoord>5e-4-y_paper/2 && ycoord<5e-4+y_paper/2) ||
         (xcoord<y_paper && fabs(ycoord-5e-4-y_paper/2)<tol && ycoord>5e-4-y_paper/2 && ycoord<5e-4+y_paper/2) ||
         (fabs(xcoord-y_paper)<tol && ycoord>5e-4-y_paper/2 && ycoord<5e-4+y_paper/2 && ycoord>5e-4-y_paper/2 && ycoord<5e-4+y_paper/2))
        {
          retval = 0;
          break;
        }
    }

  if (currentlevel <= minlevel)
    retval = 0;
      
  return (retval);
}

double epsilon_fun(const double & x, const double & y)
{   
/*
    if((y > y_paper && y<y_paper+y_oil) ||
    (y>2*y_paper+y_oil && y<2*(y_paper+y_oil)))
        return  epsilon_0 * epsilon_r_1;
    
    if((y > y_paper+y_oil && y<2*y_paper+y_oil) &&
    x<y_paper)
        return  epsilon_0;
*/   
    return epsilon_0 * epsilon_r_2;
}

double sigma_fun(const double & x, const double & y)
{
/*
  if((y > y_paper+y_oil && y<2*y_paper+y_oil) && x<y_paper)
        return  sigma_2 * DELTAT;
*/
  return sigma_1 * DELTAT;
}

std::vector<size_t> find_idx(tmesh &tmsh,
                            std::vector<std::vector<double>> &points,
                            std::vector<std::vector<double>> &tols,
                            const size_t &N_rhos)
{
  std::vector<size_t> id(N_rhos,0);

  for (size_t i = 0; i < N_rhos; i++) {
    bool found = false;

    for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant){

        for (int ii = 0; ii < 4; ++ii) {

          if (fabs(quadrant->p(0,ii)-points[i][0])<tols[i][0] &&
          fabs(quadrant->p(1,ii)-points[i][1])<tols[i][1]) {
            id[i] = quadrant->t(ii);
            found = true;
            std::cout << "Point " << i+1 << ": x= " << quadrant->p(0,ii) << ", y= " << quadrant->p(1,ii) << std::endl;
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