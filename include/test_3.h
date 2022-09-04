/*!
\file test_3.h
\brief Test case with oil and paper layers and oil filled cubic butt gaps.
*/
#include <tmesh_3d.h>
#include <simple_connectivity_3d_thin.h> // Loading library with domain geometry

constexpr int NUM_REFINEMENTS = 5;      /*!<\brief Level for global refinement */
constexpr int maxlevel = 6;             /*!<\brief Level for local refinement */
constexpr int minlevel = 4;             /*!<\brief Level for global coarsening */

constexpr double DELTAT = 50.0;         /*!<\brief Temporal time step */
constexpr double T = 5000;              /*!<\brief Final time of simulation */
constexpr double tau = 50.0;            /*!<\brief Time constant for boundary conditions */
constexpr bool save_sol = true;         /*!<\brief If set to 'true' saves data for Paraview visualization */

// Problem parameters
constexpr double epsilon_0 = 8.8542e-12;    /*!<\brief Permittivity vacuum */
constexpr double epsilon_r_1 = 2.0;         /*!<\brief Permittivity oil */
constexpr double epsilon_r_2 = 4.0;         /*!<\brief Permittivity paper */
constexpr double sigma_ = 3.21e-14;         /*!<\brief Conducivity */

constexpr double z_oil = 5e-5;              /*!<\brief Thickness oil layer */
constexpr double z_paper = 3e-4;            /*!<\brief Thickness paper layer and butt gaps side length*/
constexpr double tol = 1e-5;                /*!<\brief Tolerance value for refinement and point selection */

constexpr size_t N_rhos = 6;                /*!<\brief Number of poiunts to select for output */
std::vector<size_t> rho_idx;
std::vector<std::vector<double>> points{{5e-4,5e-4,0.0},{5e-4,5e-4,z_paper},{5e-4,5e-4,z_oil+z_paper},{5e-4,5e-4,z_oil+2*z_paper},{5e-4,5e-4,2*z_oil+2*z_paper},{5e-4,5e-4,1e-3}}; /*!<\brief Coordinates of the selected points */
std::vector<std::vector<double>> tols{{1e-4,1e-4,tol},{1e-4,1e-4,tol},{1e-4,1e-4,tol},{1e-4,1e-4,tol},{1e-4,1e-4,tol},{1e-4,1e-4,tol}}; /*!<\brief Tolerance around the selected points */
bool extra_refinement = true;       /*!<\brief true for all cases except Test 1 */

static int
uniform_refinement (tmesh_3d::quadrant_iterator q)
{ return NUM_REFINEMENTS; }     /*!<\brief Uniform refinement function */

static int
refinement (tmesh_3d::quadrant_iterator quadrant)
{
  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double zcoord;
  int retval = 0;
  for (int ii = 0; ii < 8; ++ii)
    {
      zcoord = quadrant->p(2, ii);
      
      if (zcoord > z_paper - tol || zcoord < 2*z_paper+2*z_oil+tol)
        {
          retval = maxlevel - currentlevel;
          break;
        }
    }

  if (currentlevel >= maxlevel)
    retval = 0;
      
  return retval;
}       /*!<\brief Performs local refinement at the interfaces */

static int
coarsening (tmesh_3d::quadrant_iterator quadrant)
{
  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double xcoord,ycoord,zcoord;
  int retval = currentlevel - minlevel;
  for (int ii = 0; ii < 8; ++ii)
    { 
      xcoord = quadrant->p(0, ii);
      ycoord = quadrant->p(1, ii);    
      zcoord = quadrant->p(2, ii);

      if (fabs(zcoord - z_paper) < tol || fabs(zcoord - z_paper-z_oil) < tol || fabs(zcoord - 2*z_paper-z_oil) < tol || fabs(zcoord - 2*z_paper-2*z_oil) < tol || 
         (xcoord<z_paper && fabs(ycoord-5e-4+z_paper/2)<tol && zcoord>5e-4-z_paper/2 && zcoord<5e-4+z_paper/2) ||
         (xcoord<z_paper && fabs(ycoord-5e-4-z_paper/2)<tol && zcoord>5e-4-z_paper/2 && zcoord<5e-4+z_paper/2) ||
         (fabs(xcoord-z_paper)<tol && ycoord>5e-4-z_paper/2 && ycoord<5e-4+z_paper/2 && zcoord>5e-4-z_paper/2 && zcoord<5e-4+z_paper/2))
        {
          retval = 0;
          break;
        }
    }

  if (currentlevel <= minlevel)
    retval = 0;
      
  return (retval);
}       /*!<\brief Performs local coarsening, leaves interfaces refined */

double epsilon_fun(const double & x, const double & y, const double & z)
{   
    if((z > z_paper && z<z_paper+z_oil) || (z>2*z_paper+z_oil && z<2*(z_paper+z_oil)))
        return  epsilon_0 * epsilon_r_1;
    
    if((z > z_paper+z_oil && z<2*z_paper+z_oil) && x<z_paper && (y>5e-4-z_paper/2 && y<5e-4+z_paper/2))
        return  epsilon_0 * epsilon_r_1;
    
    return epsilon_0 * epsilon_r_2;
}       /*!<\brief Return the value of epsilon at the point (x,y,z) */

double sigma_fun(const double & x, const double & y, const double & z)
{return sigma_ * DELTAT;}       /*!<\brief Return the value of sigma at the point (x,y,z) */

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
}       /*!<\brief Find the global index of the points given by the vector 'points' */
