#ifndef MI_PAPER_TEST
#define MI_PAPER_TEST

#include "generic_test.h"

using json = nlohmann::json;

namespace {
  int num_refinements;
  int maxlevel;
  int minlevel;
}

namespace tests
{
  class MI_paper : public generic_test {
    private:
      double epsilon_r_1, epsilon_r_2, epsilon_r_gap;
      double chi1, chi2, chi3, chi4, chi5;
      double tau_p1, tau_p2, tau_p3, tau_p4, tau_p5;
      double sigma_1, sigma_2, sigma_gap;
      double tol = 1e-5;
      double z_oil = 5e-5;
      double z_paper = 3e-4;
    public:

      MI_paper() {extra_refinement = true;}
      
      void import_params (json & data) {
        epsilon_r_1 = data["physics"]["plugin_params"]["epsilon_r_1"];
        epsilon_r_2 = data["physics"]["plugin_params"]["epsilon_r_2"];
        epsilon_r_gap = data["physics"]["plugin_params"]["epsilon_r_gap"];
        chi1 = data["physics"]["plugin_params"]["chi1"];
        chi2 = data["physics"]["plugin_params"]["chi2"];
        chi3 = data["physics"]["plugin_params"]["chi3"];
        chi4 = data["physics"]["plugin_params"]["chi4"];
        chi5 = data["physics"]["plugin_params"]["chi5"];
        tau_p1 = data["physics"]["plugin_params"]["tau_p1"];
        tau_p2 = data["physics"]["plugin_params"]["tau_p2"];
        tau_p3 = data["physics"]["plugin_params"]["tau_p3"];
        tau_p4 = data["physics"]["plugin_params"]["tau_p4"];
        tau_p5 = data["physics"]["plugin_params"]["tau_p5"];
        sigma_1 = data["physics"]["plugin_params"]["sigma_1"];
        sigma_2 = data["physics"]["plugin_params"]["sigma_2"];
        sigma_gap = data["physics"]["plugin_params"]["sigma_gap"];
        num_refinements = data["algorithm"]["num_refinements"];
        maxlevel = data["algorithm"]["maxlevel"];
        minlevel = data["algorithm"]["minlevel"];
        return;
      }

      int	uniform_refinement (tmesh_3d::quadrant_iterator q) const
        { return num_refinements; }

      int refinement (tmesh_3d::quadrant_iterator q) const {
        int currentlevel = static_cast<int> (q->the_quadrant->level);
        double zcoord;
        int retval = 0;
        for (int ii = 0; ii < 8; ++ii) {
          zcoord = q->p(2, ii);
            
          if (zcoord > z_paper - tol || zcoord < 2*z_paper+2*z_oil+tol) {
            retval = maxlevel - currentlevel;
            break;
          }
        }
        
        if (currentlevel >= maxlevel)
          retval = 0;
              
        return retval;
      }

      int coarsening (tmesh_3d::quadrant_iterator q) const {
        int currentlevel = static_cast<int> (q->the_quadrant->level);
        double xcoord,ycoord,zcoord;
        int retval = currentlevel - minlevel;
        for (int ii = 0; ii < 8; ++ii) {
          xcoord = q->p(0, ii);
          ycoord = q->p(1, ii);    
          zcoord = q->p(2, ii);

          if (fabs(zcoord - z_paper) < tol || fabs(zcoord - z_paper-z_oil) < tol || fabs(zcoord - 2*z_paper-z_oil) < tol || fabs(zcoord - 2*z_paper-2*z_oil) < tol || 
             (xcoord<z_paper && fabs(ycoord-5e-4+z_paper/2)<tol && zcoord>5e-4-z_paper/2 && zcoord<5e-4+z_paper/2) ||
             (xcoord<z_paper && fabs(ycoord-5e-4-z_paper/2)<tol && zcoord>5e-4-z_paper/2 && zcoord<5e-4+z_paper/2) ||
             (fabs(xcoord-z_paper)<tol && ycoord>5e-4-z_paper/2 && ycoord<5e-4+z_paper/2 && zcoord>5e-4-z_paper/2 && zcoord<5e-4+z_paper/2)) {
            retval = 0;
            break;
          }
        }

        if (currentlevel <= minlevel)
          retval = 0;
      
        return (retval);
      }
        
      double epsilon_fun (double x, double y, double z) const {   
        if((z > z_paper && z<z_paper+z_oil) || (z>2*z_paper+z_oil && z<2*(z_paper+z_oil)))
          return  epsilon_0 * epsilon_r_2;
    
        if((z > z_paper+z_oil && z<2*z_paper+z_oil) && x<z_paper && (y>5e-4-z_paper/2 && y<5e-4+z_paper/2))
          return  epsilon_0 * epsilon_r_gap;
    
        return epsilon_0 * epsilon_r_1;
      }

      double tau_p1_fun (double x, double y, double z) const {
        return tau_p1;
      }

      double tau_p2_fun (double x, double y, double z) const {
        return tau_p2;
      }

      double tau_p3_fun (double x, double y, double z) const {
        return tau_p3;
      }

      double tau_p4_fun (double x, double y, double z) const {
        return tau_p4;
      }

      double tau_p5_fun (double x, double y, double z) const {
        return tau_p5;
      }

      double chi_1_fun (double x, double y, double z) const {
        if((z > z_paper && z<z_paper+z_oil) || (z>2*z_paper+z_oil && z<2*(z_paper+z_oil)))
          return  0;
    
        if((z > z_paper+z_oil && z<2*z_paper+z_oil) && x<z_paper && (y>5e-4-z_paper/2 && y<5e-4+z_paper/2))
          return  0;
    
        return chi1;
      }

      double chi_2_fun (double x, double y, double z) const {
        if((z > z_paper && z<z_paper+z_oil) || (z>2*z_paper+z_oil && z<2*(z_paper+z_oil)))
          return  0;
    
        if((z > z_paper+z_oil && z<2*z_paper+z_oil) && x<z_paper && (y>5e-4-z_paper/2 && y<5e-4+z_paper/2))
          return  0;
    
        return chi2;
      }

      double chi_3_fun (double x, double y, double z) const {
        if((z > z_paper && z<z_paper+z_oil) || (z>2*z_paper+z_oil && z<2*(z_paper+z_oil)))
          return  0;
    
        if((z > z_paper+z_oil && z<2*z_paper+z_oil) && x<z_paper && (y>5e-4-z_paper/2 && y<5e-4+z_paper/2))
          return  0;
    
        return chi3;
      }

      
      double chi_4_fun (double x, double y, double z) const {
        if((z > z_paper && z<z_paper+z_oil) || (z>2*z_paper+z_oil && z<2*(z_paper+z_oil)))
          return  chi4;
    
        if((z > z_paper+z_oil && z<2*z_paper+z_oil) && x<z_paper && (y>5e-4-z_paper/2 && y<5e-4+z_paper/2))
          return  chi4;
    
        return 0;
      }

      double chi_5_fun (double x, double y, double z) const {
        if((z > z_paper && z<z_paper+z_oil) || (z>2*z_paper+z_oil && z<2*(z_paper+z_oil)))
          return  chi5;
    
        if((z > z_paper+z_oil && z<2*z_paper+z_oil) && x<z_paper && (y>5e-4-z_paper/2 && y<5e-4+z_paper/2))
          return  chi5;
    
        return 0;
      }

      double sigma_fun (double x, double y, double z, double DT) const {
         if((z > z_paper+z_oil && z<2*z_paper+z_oil) && x<z_paper && (y>5e-4-z_paper/2 && y<5e-4+z_paper/2))
          return  sigma_2 * DT;
        
        if((z > z_paper && z<z_paper+z_oil) || (z>2*z_paper+z_oil && z<2*(z_paper+z_oil)))
          return  sigma_gap * DT;
    
        return sigma_1 * DT;
      }
  };

} // namespace tests

#endif