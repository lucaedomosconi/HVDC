#ifndef PARTICULAR_TEST_HPP
#define PARTICULAR_TEST_HPP
#include "generic_test.h"
using json = nlohmann::json;

namespace {
  int num_refinements = 4;
  int maxlevel = 6;
}

namespace tests {

  class homogeneous : public generic_test {
    private:
      double epsilon_r1;		                // permittivity at infinite frequency
      double chi1, chi2, chi3;
      double tau_p1, tau_p2, tau_p3;
      double sigma_;            					    // conductivity coeff
    public:

      homogeneous() {extra_refinement = false;}
      
      void import_params (json & data) {
        try{epsilon_r1 = data["physics"]["plugin_params"]["epsilon_r1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][epsilon_r1]");}
        try{chi1 = data["physics"]["plugin_params"]["chi1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi1]");}
        try{chi2 = data["physics"]["plugin_params"]["chi2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi2]");}
        try{chi3 = data["physics"]["plugin_params"]["chi3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi3]");}
        try{tau_p1 = data["physics"]["plugin_params"]["tau_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_p1]");}
        try{tau_p2 = data["physics"]["plugin_params"]["tau_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_p2]");}
        try{tau_p3 = data["physics"]["plugin_params"]["tau_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_p3]");}
        try{sigma_ = data["physics"]["plugin_params"]["sigma"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][sigma]");}
        try{num_refinements = data["algorithm"]["num_refinements"];}
        catch(...){throw std::runtime_error("[algorithm][num_refinements]");}
        try{maxlevel = data["algorithm"]["maxlevel"];}
        catch(...){throw std::runtime_error("[algorithm][maxlevel]");}
        
        return;
      }

      int	uniform_refinement (tmesh_3d::quadrant_iterator q) const
        { return num_refinements; }

      int refinement (tmesh_3d::quadrant_iterator q) const
        { return num_refinements; }

      int coarsening (tmesh_3d::quadrant_iterator q) const
        { return num_refinements; }
        
      double epsilon_fun (double x, double y, double z) const
        {return epsilon_0 * epsilon_r1;}

      double chi_1_fun (double x, double y, double z) const
        {return chi1;}

      double chi_2_fun (double x, double y, double z) const
        {return chi2;}

      double chi_3_fun (double x, double y, double z) const
        {return chi3;}

      double chi_4_fun (double x, double y, double z) const
        {return 0;}

      double chi_5_fun (double x, double y, double z) const
        {return 0;}

      double tau_p1_fun (double x, double y, double z) const
        {return tau_p1;}

      double tau_p2_fun (double x, double y, double z) const
        {return tau_p2;}

      double tau_p3_fun (double x, double y, double z) const
        {return tau_p3;}

      double tau_p4_fun (double x, double y, double z) const
        {return 1;}

      double tau_p5_fun (double x, double y, double z) const
        {return 1;}

      double sigma_fun (double x, double y, double z, double DT, double E = 0) const
        {return sigma_ * DT;}
  };

  class two_phase_serial : public generic_test {
    private:
      double epsilon_r1, epsilon_r2;		// permittivity at infinite frequency
      double chi_m1_p1, chi_m2_p1, chi_m1_p2, chi_m2_p2, chi_m1_p3, chi_m2_p3, chi_m1_p4, chi_m2_p4, chi_m1_p5, chi_m2_p5; 
      double tau_m1_p1, tau_m2_p1, tau_m1_p2, tau_m2_p2, tau_m1_p3, tau_m2_p3, tau_m1_p4, tau_m2_p4, tau_m1_p5, tau_m2_p5;
      double sigma_1, sigma_2;            					    // conductivity coeff
    public:

      two_phase_serial() {extra_refinement = true;}
      
      void import_params (json & data) {
        try{epsilon_r1 = data["physics"]["plugin_params"]["epsilon_r1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][epsilon_r1]");}
        try{epsilon_r2 = data["physics"]["plugin_params"]["epsilon_r2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][epsilon_r1]");}
        try{chi_m1_p1 = data["physics"]["plugin_params"]["chi_m1_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p1]");}
        try{chi_m2_p1 = data["physics"]["plugin_params"]["chi_m2_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p1]");}
        try{chi_m1_p2 = data["physics"]["plugin_params"]["chi_m1_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p2]");}
        try{chi_m2_p2 = data["physics"]["plugin_params"]["chi_m2_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p2]");}
        try{chi_m1_p3 = data["physics"]["plugin_params"]["chi_m1_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p3]");}
        try{chi_m2_p3 = data["physics"]["plugin_params"]["chi_m2_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p3]");}
        try{chi_m1_p4 = data["physics"]["plugin_params"]["chi_m1_p4"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p4]");}
        try{chi_m2_p4 = data["physics"]["plugin_params"]["chi_m2_p4"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p4]");}
        try{chi_m1_p5 = data["physics"]["plugin_params"]["chi_m1_p5"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p5]");}
        try{chi_m2_p5 = data["physics"]["plugin_params"]["chi_m2_p5"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p5]");}
        try{tau_m1_p1 = data["physics"]["plugin_params"]["tau_m1_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p1]");}
        try{tau_m2_p1 = data["physics"]["plugin_params"]["tau_m2_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p1]");}
        try{tau_m1_p2 = data["physics"]["plugin_params"]["tau_m1_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p2]");}
        try{tau_m2_p2 = data["physics"]["plugin_params"]["tau_m2_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p2]");}
        try{tau_m1_p3 = data["physics"]["plugin_params"]["tau_m1_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p3]");}
        try{tau_m2_p3 = data["physics"]["plugin_params"]["tau_m2_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p3]");}
        try{tau_m1_p4 = data["physics"]["plugin_params"]["tau_m1_p4"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p4]");}
        try{tau_m2_p4 = data["physics"]["plugin_params"]["tau_m2_p4"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p4]");}
        try{tau_m1_p5 = data["physics"]["plugin_params"]["tau_m1_p5"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p5]");}
        try{tau_m2_p5 = data["physics"]["plugin_params"]["tau_m2_p5"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p5]");}
        try{sigma_1 = data["physics"]["plugin_params"]["sigma1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][sigma1]");}
        try{sigma_2 = data["physics"]["plugin_params"]["sigma2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][sigma2]");}
        try{num_refinements = data["algorithm"]["num_refinements"];}
        catch(...){throw std::runtime_error("[algorithm][num_refinements]");}
        try{maxlevel = data["algorithm"]["maxlevel"];}
        catch(...){throw std::runtime_error("[algorithm][maxlevel]");}
        
        return;
      }

      int	uniform_refinement (tmesh_3d::quadrant_iterator q) const
        { return num_refinements; }

      int refinement (tmesh_3d::quadrant_iterator q) const
        {
          int currentlevel = static_cast<int> (q->the_quadrant->level);
            double zcoord;
            int retval = 0;
            for (int ii = 0; ii < 8; ++ii)
              {
                zcoord = q->p(2, ii);

                if (fabs(zcoord - 0.0005) < 1e-9)
                  {
                    retval = maxlevel - currentlevel;
                    break;
                  }
              }

            if (currentlevel >= maxlevel)
              retval = 0;

            return retval;
        }

      int coarsening (tmesh_3d::quadrant_iterator q) const
        {
          int currentlevel = static_cast<int> (q->the_quadrant->level);
            double zcoord;
            int retval = currentlevel - num_refinements;
            for (int ii = 0; ii < 8; ++ii)
              {     
                  zcoord = q->p(2, ii);

                  if (fabs(zcoord - 0.0005) < 1e-9)
                  {
                      retval = 0;
                      break;
                  }
              }

          if (currentlevel <= num_refinements)
              retval = 0;
      
            return (retval);
        }
        
      double epsilon_fun (double x, double y, double z) const
        {return z < 0.0005 ? epsilon_0 * epsilon_r1 : epsilon_0 * epsilon_r2;}

      double chi_1_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? chi_m1_p1 : chi_m2_p1;
        }

      double chi_2_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? chi_m1_p2 : chi_m2_p2;
        }

      double chi_3_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? chi_m1_p3 : chi_m2_p3;
        }

      double chi_4_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? chi_m1_p4 : chi_m2_p4;
        }
      
      double chi_5_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? chi_m1_p5 : chi_m2_p5;
        }

      double tau_p1_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? tau_m1_p1 : tau_m2_p1;
        }

      double tau_p2_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? tau_m1_p2 : tau_m2_p2;
        }

      double tau_p3_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? tau_m1_p3 : tau_m2_p3;
        }

      double tau_p4_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? tau_m1_p4 : tau_m2_p4;
        }

      double tau_p5_fun (double x, double y, double z) const
        {
          return z < 0.0005 ? tau_m1_p5 : tau_m2_p5;
        }

      double sigma_fun (double x, double y, double z, double DT, double E = 0) const
        {return z < 0.0005 ? sigma_1 * DT : sigma_2 * DT;}
  };
  
  class two_phase_parallel : public generic_test {
    private:
      double epsilon_r1, epsilon_r2;		// permittivity at infinite frequency
      double chi1, chi2, chi3;
      double tau1_p1, tau2_p1, tau_p2, tau_p3;
      double sigma_1, sigma_2;							// conductivity coeff
    public:

      two_phase_parallel() {extra_refinement = true;}
      
      void import_params (json & data) {
        try{epsilon_r1 = data["physics"]["plugin_params"]["epsilon_r1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][epsilon_r1]");}
        try{epsilon_r2 = data["physics"]["plugin_params"]["epsilon_r2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][epsilon_r1]");}
        try{chi1 = data["physics"]["plugin_params"]["chi1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi1]");}
        try{chi2 = data["physics"]["plugin_params"]["chi2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi2]");}
        try{chi3 = data["physics"]["plugin_params"]["chi3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi3]");}
        try{tau1_p1 = data["physics"]["plugin_params"]["tau1_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau1_p1]");}
        try{tau2_p1 = data["physics"]["plugin_params"]["tau2_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau2_p1]");}
        try{tau_p2 = data["physics"]["plugin_params"]["tau_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_p2]");}
        try{tau_p3 = data["physics"]["plugin_params"]["tau_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_p3]");}
        try{sigma_1 = data["physics"]["plugin_params"]["sigma1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][sigma1]");}
        try{sigma_2 = data["physics"]["plugin_params"]["sigma2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][sigma2]");}
        try{num_refinements = data["algorithm"]["num_refinements"];}
        catch(...){throw std::runtime_error("[algorithm][num_refinements]");}
        try{maxlevel = data["algorithm"]["maxlevel"];}
        catch(...){throw std::runtime_error("[algorithm][maxlevel]");}
        
        return;
      }

      int	uniform_refinement (tmesh_3d::quadrant_iterator q) const
        { return num_refinements; }

      int refinement (tmesh_3d::quadrant_iterator q) const
        {
          int currentlevel = static_cast<int> (q->the_quadrant->level);
            double xcoord, zcoord;
            int retval = 0;
            for (int ii = 0; ii < 8; ++ii)
              {
                xcoord = q->p(0, ii);
                zcoord = q->p(2, ii);

                if (fabs(xcoord - 0.0005) < 1e-9 /*|| fabs(zcoord) < 1e-9 || fabs(zcoord - 1e-3) < 1e-9*/)
                  {
                    retval = maxlevel - currentlevel;
                    break;
                  }
              }

            if (currentlevel >= maxlevel)
              retval = 0;

            return retval;
        }

      int coarsening (tmesh_3d::quadrant_iterator q) const
        {
          int currentlevel = static_cast<int> (q->the_quadrant->level);
            double xcoord, zcoord;
            int retval = currentlevel - num_refinements;
            for (int ii = 0; ii < 8; ++ii)
              {
                xcoord = q->p(0, ii);
                zcoord = q->p(2, ii);

                  if (fabs(xcoord - 0.0005) < 1e-9 /*|| fabs(zcoord) < 1e-9 || fabs(zcoord - 1e-3) < 1e-9*/)
                    {
                      retval = 0;
                      break;
                    }
              }

          if (currentlevel <= num_refinements)
              retval = 0;
      
            return (retval);
        }
        
      double epsilon_fun (double x, double y, double z) const
        {return x < 0.0005 ? epsilon_0 * epsilon_r1 : epsilon_0 * epsilon_r2;}

      double chi_1_fun (double x, double y, double z) const
        {return chi1;}

      double chi_2_fun (double x, double y, double z) const
        {return chi2;}

      double chi_3_fun (double x, double y, double z) const
        {return chi3;}

      double chi_4_fun (double x, double y, double z) const
        {return 0;}

      double chi_5_fun (double x, double y, double z) const
        {return 0;}

      double tau_p1_fun (double x, double y, double z) const
        {return x < 0.0005 ? tau1_p1 : tau2_p1;}

      double tau_p2_fun (double x, double y, double z) const
        {return tau_p2;}

      double tau_p3_fun (double x, double y, double z) const
        {return tau_p3;}

      double tau_p4_fun (double x, double y, double z) const
        {return 1;}

      double tau_p5_fun (double x, double y, double z) const
        {return 1;}

      double sigma_fun (double x, double y, double z, double DT, double E = 0) const
        {return x < 0.0005 ? sigma_1 * DT : sigma_2 * DT;}
  };



  class hole : public generic_test {
    private:
      double epsilon_r1, epsilon_r2;		// permittivity at infinite frequency
      double chi_m1_p1, chi_m2_p1, chi_m1_p2, chi_m2_p2, chi_m1_p3, chi_m2_p3, chi_m1_p4, chi_m2_p4, chi_m1_p5, chi_m2_p5; 
      double tau_m1_p1, tau_m2_p1, tau_m1_p2, tau_m2_p2, tau_m1_p3, tau_m2_p3, tau_m1_p4, tau_m2_p4, tau_m1_p5, tau_m2_p5;
      static constexpr double hole_rad_2 = 1e-8;
      double sigma_1, sigma_2;            					    // conductivity coeff
    public:

      hole() {extra_refinement = true;}
      
      void import_params (json & data) {
        try{epsilon_r1 = data["physics"]["plugin_params"]["epsilon_r1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][epsilon_r1]");}
        try{epsilon_r2 = data["physics"]["plugin_params"]["epsilon_r2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][epsilon_r1]");}
        try{chi_m1_p1 = data["physics"]["plugin_params"]["chi_m1_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p1]");}
        try{chi_m2_p1 = data["physics"]["plugin_params"]["chi_m2_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p1]");}
        try{chi_m1_p2 = data["physics"]["plugin_params"]["chi_m1_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p2]");}
        try{chi_m2_p2 = data["physics"]["plugin_params"]["chi_m2_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p2]");}
        try{chi_m1_p3 = data["physics"]["plugin_params"]["chi_m1_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p3]");}
        try{chi_m2_p3 = data["physics"]["plugin_params"]["chi_m2_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p3]");}
        try{chi_m1_p4 = data["physics"]["plugin_params"]["chi_m1_p4"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p4]");}
        try{chi_m2_p4 = data["physics"]["plugin_params"]["chi_m2_p4"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p4]");}
        try{chi_m1_p5 = data["physics"]["plugin_params"]["chi_m1_p5"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m1_p5]");}
        try{chi_m2_p5 = data["physics"]["plugin_params"]["chi_m2_p5"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][chi_m2_p5]");}
        try{tau_m1_p1 = data["physics"]["plugin_params"]["tau_m1_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p1]");}
        try{tau_m2_p1 = data["physics"]["plugin_params"]["tau_m2_p1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p1]");}
        try{tau_m1_p2 = data["physics"]["plugin_params"]["tau_m1_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p2]");}
        try{tau_m2_p2 = data["physics"]["plugin_params"]["tau_m2_p2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p2]");}
        try{tau_m1_p3 = data["physics"]["plugin_params"]["tau_m1_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p3]");}
        try{tau_m2_p3 = data["physics"]["plugin_params"]["tau_m2_p3"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p3]");}
        try{tau_m1_p4 = data["physics"]["plugin_params"]["tau_m1_p4"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p4]");}
        try{tau_m2_p4 = data["physics"]["plugin_params"]["tau_m2_p4"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p4]");}
        try{tau_m1_p5 = data["physics"]["plugin_params"]["tau_m1_p5"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m1_p5]");}
        try{tau_m2_p5 = data["physics"]["plugin_params"]["tau_m2_p5"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][tau_m2_p5]");}
        try{sigma_1 = data["physics"]["plugin_params"]["sigma1"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][sigma1]");}
        try{sigma_2 = data["physics"]["plugin_params"]["sigma2"];}
        catch(...){throw std::runtime_error("[physics][plugin_params][sigma2]");}
        try{num_refinements = data["algorithm"]["num_refinements"];}
        catch(...){throw std::runtime_error("[algorithm][num_refinements]");}
        try{maxlevel = data["algorithm"]["maxlevel"];}
        catch(...){throw std::runtime_error("[algorithm][maxlevel]");}
        
        return;
      }

      int	uniform_refinement (tmesh_3d::quadrant_iterator q) const
        { return num_refinements; }

      int refinement (tmesh_3d::quadrant_iterator q) const
        {
          int currentlevel = static_cast<int> (q->the_quadrant->level);
            double xcoord, ycoord, zcoord;
            int retval = 0;
            for (int ii = 0; ii < 8; ++ii)
              {
                xcoord = q->p(0, ii);
                ycoord = q->p(1, ii);
                zcoord = q->p(2, ii);

                if ((xcoord-2.5e-3)*(xcoord-2.5e-3)+(ycoord-2.5e-3)*(ycoord-2.5e-3)+(zcoord-0.5e-3)*(zcoord-0.5e-3) < 1.1*hole_rad_2)
                  {
                    retval = maxlevel - currentlevel;
                    break;
                  }
              }

            if (currentlevel >= maxlevel)
              retval = 0;

            return retval;
        }

      int coarsening (tmesh_3d::quadrant_iterator q) const
        {
          int currentlevel = static_cast<int> (q->the_quadrant->level);
            double xcoord, ycoord, zcoord;
            int retval = currentlevel - num_refinements;
            for (int ii = 0; ii < 8; ++ii)
              {     
                xcoord = q->p(0, ii);
                ycoord = q->p(1, ii);
                zcoord = q->p(2, ii);


                  if ((xcoord-2.5e-3)*(xcoord-2.5e-3)+(ycoord-2.5e-3)*(ycoord-2.5e-3)+(zcoord-0.5e-3)*(zcoord-0.5e-3) < 0.8*hole_rad_2)
                  {
                      retval = 0;
                      break;
                  }
              }

          if (currentlevel <= num_refinements)
              retval = 0;
      
            return (retval);
        }
        
      double epsilon_fun (double x, double y, double z) const
        {return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? epsilon_0 * epsilon_r1 : epsilon_0 * epsilon_r2;}

      double chi_1_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? chi_m1_p1 : chi_m2_p1;
        }

      double chi_2_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? chi_m1_p2 : chi_m2_p2;
        }

      double chi_3_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? chi_m1_p3 : chi_m2_p3;
        }

      double chi_4_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? chi_m1_p4 : chi_m2_p4;
        }
      
      double chi_5_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? chi_m1_p5 : chi_m2_p5;
        }

      double tau_p1_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? tau_m1_p1 : tau_m2_p1;
        }

      double tau_p2_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? tau_m1_p2 : tau_m2_p2;
        }

      double tau_p3_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? tau_m1_p3 : tau_m2_p3;
        }

      double tau_p4_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? tau_m1_p4 : tau_m2_p4;
        }

      double tau_p5_fun (double x, double y, double z) const
        {
          return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? tau_m1_p5 : tau_m2_p5;
        }

      double sigma_fun (double x, double y, double z, double DT, double E = 0) const
        {return (x-2.5e-3)*(x-2.5e-3)+(y-2.5e-3)*(y-2.5e-3)+(z-0.5e-3)*(z-0.5e-3) < hole_rad_2 ? sigma_1 * DT : sigma_2 * DT;}
  };
  
}

#endif