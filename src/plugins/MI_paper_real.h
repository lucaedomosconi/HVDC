#ifndef MI_PAPER_REAL_TEST
#define MI_PAPER_REAL_TEST

#include "generic_test.h"

using json = nlohmann::json;

namespace {
  int num_refinements;
  int maxlevel;
  int minlevel;
}


#define MI_HOLE 2

namespace tests
{
  class MI_paper_real : public generic_test {
    private:
      double epsilon_r_paper, epsilon_r_oil, epsilon_r_bubble;
      double chi1, chi2, chi3, chi4, chi5, chi6;
      double tau_p1, tau_p2, tau_p3, tau_p4, tau_p5, tau_p6;
      double sigma_paper, sigma_oil, sigma_bubble;
      double eta = 2.52e-8;
      double tol = 1.5e-5;
      double z_oil = 5e-5;
      double half_gap_height = 5e-5;
      double half_gap_length = 5e-4;
      double radius = 4e-5;
      double displ_bubbles = 1.5*radius;
      double z_paper = 0;

      int refinement_1 (tmesh_3d::quadrant_iterator q) const {
        int currentlevel = static_cast<int> (q->the_quadrant->level);
        double xcoord, ycoord, zcoord;
        int retval = 0;
        for (int ii = 0; ii < 8; ++ii) {
          xcoord = q->p(0, ii);
          ycoord = q->p(1, ii);
          zcoord = q->p(2, ii);
            
          if ((std::fabs(zcoord) > half_gap_height - tol && std::fabs(zcoord) < half_gap_height + tol && std::fabs(xcoord) < half_gap_length + tol) ||
              (std::fabs(zcoord) < half_gap_height + tol && std::fabs(xcoord) > half_gap_length - tol && std::fabs(xcoord) < half_gap_length + tol)
#if MI_HOLE == 1
                || (xcoord*xcoord + ycoord*ycoord + zcoord*zcoord < radius*radius)
#elif MI_HOLE == 2
                || ((xcoord-displ_bubbles)*(xcoord-displ_bubbles) + ycoord*ycoord + zcoord*zcoord < radius*radius)
                || ((xcoord+displ_bubbles)*(xcoord+displ_bubbles) + ycoord*ycoord + zcoord*zcoord < radius*radius)
#endif
              ) {
            return 1;
          }
        }
        return retval;
      }

      int refinement_2 (tmesh_3d::quadrant_iterator q) const {
        int currentlevel = static_cast<int> (q->the_quadrant->level);
        double xcoord, ycoord, zcoord;
        int retval = 0;
        for (int ii = 0; ii < 8; ++ii) {
          xcoord = q->p(0, ii);
          ycoord = q->p(1, ii);
          zcoord = q->p(2, ii);
#if MI_HOLE == 1
          if ((xcoord*xcoord + ycoord*ycoord + zcoord*zcoord < 1.2 * radius*radius) &&
              (xcoord*xcoord + ycoord*ycoord + zcoord*zcoord > 0.5 * radius*radius)) {
            return 1;
          }
#elif MI_HOLE == 2
          if ((((xcoord-displ_bubbles)*(xcoord-displ_bubbles) + ycoord*ycoord + zcoord*zcoord < 1.2 * radius*radius) &&
               ((xcoord-displ_bubbles)*(xcoord-displ_bubbles) + ycoord*ycoord + zcoord*zcoord > 0.5 * radius*radius)) ||
              (((xcoord+displ_bubbles)*(xcoord+displ_bubbles) + ycoord*ycoord + zcoord*zcoord < 1.2 * radius*radius) &&
               ((xcoord+displ_bubbles)*(xcoord+displ_bubbles) + ycoord*ycoord + zcoord*zcoord > 0.5 * radius*radius))) {
                return 1;
              }
#endif
        }
        return retval;
      }

    public:

      MI_paper_real() {extra_refinement = true;}
      
      void import_params (json & data) {
        epsilon_r_paper = data["physics"]["plugin_params"]["epsilon_r_paper"];
        epsilon_r_oil = data["physics"]["plugin_params"]["epsilon_r_oil"];
        epsilon_r_bubble = data["physics"]["plugin_params"]["epsilon_r_bubble"];
        chi1 = data["physics"]["plugin_params"]["chi1"];
        chi2 = data["physics"]["plugin_params"]["chi2"];
        chi3 = data["physics"]["plugin_params"]["chi3"];
        chi4 = data["physics"]["plugin_params"]["chi4"];
        chi5 = data["physics"]["plugin_params"]["chi5"];
        chi6 = data["physics"]["plugin_params"]["chi6"];
        tau_p1 = data["physics"]["plugin_params"]["tau_p1"];
        tau_p2 = data["physics"]["plugin_params"]["tau_p2"];
        tau_p3 = data["physics"]["plugin_params"]["tau_p3"];
        tau_p4 = data["physics"]["plugin_params"]["tau_p4"];
        tau_p5 = data["physics"]["plugin_params"]["tau_p5"];
        tau_p6 = data["physics"]["plugin_params"]["tau_p6"];
        sigma_paper = data["physics"]["plugin_params"]["sigma_paper"];
        sigma_oil = data["physics"]["plugin_params"]["sigma_oil"];
        sigma_bubble = data["physics"]["plugin_params"]["sigma_bubble"];
        num_refinements = data["algorithm"]["num_refinements"];
        maxlevel = data["algorithm"]["maxlevel"];
        minlevel = data["algorithm"]["minlevel"];
        return;
      }

      int uniform_refinement (tmesh_3d::quadrant_iterator q) const
        { return num_refinements; }

      int refinement (tmesh_3d &tmsh) const {
        for (auto ii = 0; ii < 2; ++ii) {
          tmsh.set_refine_marker([this](tmesh_3d::quadrant_iterator q) {return this->refinement_1(q);});
          tmsh.refine (1);
        }
#ifdef MI_HOLE
        for (auto ii = 0; ii < 1; ++ii) {
          tmsh.set_refine_marker([this](tmesh_3d::quadrant_iterator q) {return this->refinement_2(q);});
          tmsh.refine (1);
        }
#endif
        return 0;
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
             (fabs(xcoord-5e-4+z_paper/2)<tol && zcoord>5e-4-z_paper/2 && zcoord<5e-4+z_paper/2) ||
             (fabs(xcoord-5e-4-z_paper/2)<tol && zcoord>5e-4-z_paper/2 && zcoord<5e-4+z_paper/2) ||
             (xcoord>5e-4-z_paper/2 && xcoord<5e-4+z_paper/2 && zcoord>5e-4-z_paper/2 && zcoord<5e-4+z_paper/2)) {
            retval = 0;
            break;
          }
        }

        if (currentlevel <= minlevel)
          retval = 0;
        retval = 0;
        return (retval);
      }
        
      double epsilon_fun (double x, double y, double z) const {
#if MI_HOLE == 1
        if(x*x+y*y+z*z < radius*radius)
          return  epsilon_0 * epsilon_r_bubble;
#elif MI_HOLE == 2
        if((x-displ_bubbles)*(x-displ_bubbles)+y*y+z*z < radius*radius ||
           (x+displ_bubbles)*(x+displ_bubbles)+y*y+z*z < radius*radius)
          return  epsilon_0 * epsilon_r_bubble;
#endif
        if(std::fabs(z) < half_gap_height && std::fabs(x) < half_gap_length)
          return  epsilon_0 * epsilon_r_oil;
    
        return epsilon_0 * epsilon_r_paper;
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

      double tau_p6_fun (double x, double y, double z) const {
        return tau_p6;
      }

      double chi_1_fun (double x, double y, double z) const {
#if MI_HOLE == 1
        if(x*x+y*y+z*z < radius*radius)
          return 0;
#elif MI_HOLE == 2
        if((x-displ_bubbles)*(x-displ_bubbles)+y*y+z*z < radius*radius ||
           (x+displ_bubbles)*(x+displ_bubbles)+y*y+z*z < radius*radius)
          return 0;
#endif
    
        if(x*x < half_gap_length*half_gap_length && z*z < half_gap_height*half_gap_height)
          return  0;
    
        return chi1;
      }

      double chi_2_fun (double x, double y, double z) const {
#if MI_HOLE == 1
        if(x*x+y*y+z*z < radius*radius)
          return 0;
#elif MI_HOLE == 2
        if((x-displ_bubbles)*(x-displ_bubbles)+y*y+z*z < radius*radius ||
           (x+displ_bubbles)*(x+displ_bubbles)+y*y+z*z < radius*radius)
          return 0;
#endif
    
        if(x*x < half_gap_length*half_gap_length && z*z < half_gap_height*half_gap_height)
          return 0;
    
        return chi2;
      }

      double chi_3_fun (double x, double y, double z) const {
#if MI_HOLE == 1
        if(x*x+y*y+z*z < radius*radius)
          return 0;
#elif MI_HOLE == 2
        if((x-displ_bubbles)*(x-displ_bubbles)+y*y+z*z < radius*radius ||
           (x+displ_bubbles)*(x+displ_bubbles)+y*y+z*z < radius*radius)
          return 0;
#endif
    
        if(x*x < half_gap_length*half_gap_length && z*z < half_gap_height*half_gap_height)
          return 0;
    
        return chi3;
      }

      
      double chi_4_fun (double x, double y, double z) const {
#if MI_HOLE == 1
        if(x*x+y*y+z*z < radius*radius)
          return 0;
#elif MI_HOLE == 2
        if((x-displ_bubbles)*(x-displ_bubbles)+y*y+z*z < radius*radius ||
           (x+displ_bubbles)*(x+displ_bubbles)+y*y+z*z < radius*radius)
          return 0;
#endif
    
        if(x*x < half_gap_length*half_gap_length && z*z < half_gap_height*half_gap_height)
          return chi4;
    
        return 0;
      }

      double chi_5_fun (double x, double y, double z) const {
#if MI_HOLE == 1
        if(x*x+y*y+z*z < radius*radius)
          return 0;
#elif MI_HOLE == 2
        if((x-displ_bubbles)*(x-displ_bubbles)+y*y+z*z < radius*radius ||
           (x+displ_bubbles)*(x+displ_bubbles)+y*y+z*z < radius*radius)
          return 0;
#endif
    
        if(x*x < half_gap_length*half_gap_length && z*z < half_gap_height*half_gap_height)
          return chi5;
    
        return 0;
      }

      double chi_6_fun (double x, double y, double z) const {
#if MI_HOLE == 1
        if(x*x+y*y+z*z < radius*radius)
          return 0;
#elif MI_HOLE == 2
        if((x-displ_bubbles)*(x-displ_bubbles)+y*y+z*z < radius*radius ||
           (x+displ_bubbles)*(x+displ_bubbles)+y*y+z*z < radius*radius)
          return 0;
#endif
    
        if(x*x < half_gap_length*half_gap_length && z*z < half_gap_height*half_gap_height)
          return chi6;
    
        return 0;
      }

      double sigma_fun (double x, double y, double z, double DT, double E = 0) const {
#if MI_HOLE == 1
        if(x*x+y*y+z*z < radius*radius)
          return sigma_bubble * DT;
#elif MI_HOLE == 2
        if((x-displ_bubbles)*(x-displ_bubbles)+y*y+z*z < radius*radius ||
           (x+displ_bubbles)*(x+displ_bubbles)+y*y+z*z < radius*radius)
          return sigma_bubble * DT;
#endif
        if(x*x < half_gap_length*half_gap_length && z*z < half_gap_height*half_gap_height)
          return sigma_oil * DT;
    
        return DT * sigma_paper * std::exp(eta * E);
      }
  };

} // namespace tests

#endif