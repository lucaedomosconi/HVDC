#ifndef PARTICULAR_TEST_HPP
#define PARTICULAR_TEST_HPP
#include "generic_test.h"
using json = nlohmann::json;

namespace {
  int num_refinements = 4;
  int maxlevel = 6;
}

namespace tests {

  class test1 : public generic_test {
    private:
      double epsilon_r1;		                // permittivity at infinite frequency
      double chi1, chi2, chi3;
      double tau_p1, tau_p2, tau_p3;
      double sigma_;            					    // conductivity coeff
    public:

      test1() {extra_refinement = false;}
      
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

      double tau_p1_fun (double x, double y, double z) const
        {return tau_p1;}

      double tau_p2_fun (double x, double y, double z) const
        {return tau_p2;}

      double tau_p3_fun (double x, double y, double z) const
        {return tau_p3;}

      double sigma_fun (double x, double y, double z, double DT) const
        {return sigma_ * DT;}
  };

  class test2 : public generic_test {
    private:
      double epsilon_r1, epsilon_r2;		// permittivity at infinite frequency
      double chi1, chi2, chi3;
      double tau_p1, tau_p2, tau_p3;
      double sigma_;            					    // conductivity coeff
    public:

      test2() {extra_refinement = true;}
      
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
        {return chi1;}

      double chi_2_fun (double x, double y, double z) const
        {return chi2;}

      double chi_3_fun (double x, double y, double z) const
        {return chi3;}

      double tau_p1_fun (double x, double y, double z) const
        {return tau_p1;}

      double tau_p2_fun (double x, double y, double z) const
        {return tau_p2;}

      double tau_p3_fun (double x, double y, double z) const
        {return tau_p3;}

      double sigma_fun (double x, double y, double z, double DT) const
        {return sigma_ * DT;}
  };
  
  class test3 : public generic_test {
    private:
      double epsilon_r1, epsilon_r2;		// permittivity at infinite frequency
      double chi1, chi2, chi3;
      double tau_p1, tau_p2, tau_p3;
      double sigma_;							// conductivity coeff
    public:

      test3() {extra_refinement = true;}
      
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
        {
          int currentlevel = static_cast<int> (q->the_quadrant->level);
            double xcoord, zcoord;
            int retval = 0;
            for (int ii = 0; ii < 8; ++ii)
              {
                xcoord = q->p(0, ii);
                zcoord = q->p(2, ii);

                if (fabs(xcoord - 0.0005) < 1e-9 || fabs(zcoord) < 1e-9 || fabs(zcoord - 1e-3) < 1e-9)
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

                  if (fabs(xcoord - 0.0005) < 1e-9 || fabs(zcoord) < 1e-9 || fabs(zcoord - 1e-3) < 1e-9)
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

      double tau_p1_fun (double x, double y, double z) const
        {return tau_p1;}

      double tau_p2_fun (double x, double y, double z) const
        {return tau_p2;}

      double tau_p3_fun (double x, double y, double z) const
        {return tau_p3;}

      double sigma_fun (double x, double y, double z, double DT) const
        {return sigma_ * DT;}
  };
}

#endif