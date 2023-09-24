#ifndef PARTICULAR_TEST_HPP
#define PARTICULAR_TEST_HPP
#include "test_factory.h"
#include "generic_test.h"
using json = nlohmann::json;

int NUM_REFINEMENTS = 4;
int maxlevel = 6;




namespace tests {

  class test1 : public generic_test {
    private:
      double epsilon_inf_1;		// permittivity at infinite frequency
      double csi1, csi2, csi3;
      double tau_p1, tau_p2, tau_p3;
      double sigma_;            					    // conducivity coeff
    public:

			test1() {extra_refinement = false;}
			
			void import_params (const json &data) {
				epsilon_inf_1 = data["test1"]["physics_grid"]["plugin_params"]["epsilon_inf_1"];
				csi1 = data["test1"]["physics_grid"]["plugin_params"]["csi1"];
				csi2 = data["test1"]["physics_grid"]["plugin_params"]["csi2"];
				csi3 = data["test1"]["physics_grid"]["plugin_params"]["csi3"];
        tau_p1 = data["test1"]["physics_grid"]["plugin_params"]["tau_p1"];
        tau_p2 = data["test1"]["physics_grid"]["plugin_params"]["tau_p2"];
        tau_p3 = data["test1"]["physics_grid"]["plugin_params"]["tau_p3"];
				sigma_ = data["test1"]["physics_grid"]["plugin_params"]["sigma"];
        NUM_REFINEMENTS = data["test1"]["algorithm"]["NUM_REFINEMENTS"];
        maxlevel = data["test1"]["algorithm"]["maxlevel"];
				return;
			}

			int	uniform_refinement (tmesh_3d::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }

			int refinement (tmesh_3d::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }

			int coarsening (tmesh_3d::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }
				
			double epsilon_fun (const double & x, const double & y, const double & z) const
				{return epsilon_0 * epsilon_inf_1;}

			double csi_1_fun (const double & x, const double & y, const double & z) const
				{return csi1;}

			double csi_2_fun (const double & x, const double & y, const double & z) const
				{return csi2;}

			double csi_3_fun (const double & x, const double & y, const double & z) const
				{return csi3;}

      double tau_p1_fun (const double & x, const double & y, const double & z) const
				{return tau_p1;}

			double tau_p2_fun (const double & x, const double & y, const double & z) const
				{return tau_p2;}

			double tau_p3_fun (const double & x, const double & y, const double & z) const
				{return tau_p3;}

			double sigma_fun (const double & x, const double & y, const double & z, const double & DT) const
				{return sigma_ * DT;}
  };

	class test2 : public generic_test {
    private:
      double epsilon_inf_1, epsilon_inf_2;		// permittivity at infinite frequency
      double csi1, csi2, csi3;
      double tau_p1, tau_p2, tau_p3;
      double sigma_;            					    // conducivity coeff
    public:

			test2() {extra_refinement = true;}
			
			void import_params (const json &data) {
				epsilon_inf_1 = data["test2"]["physics_grid"]["plugin_params"]["epsilon_inf_1"];
				epsilon_inf_2 = data["test2"]["physics_grid"]["plugin_params"]["epsilon_inf_2"];
				csi1 = data["test2"]["physics_grid"]["plugin_params"]["csi1"];
				csi2 = data["test2"]["physics_grid"]["plugin_params"]["csi2"];
				csi3 = data["test2"]["physics_grid"]["plugin_params"]["csi3"];
        tau_p1 = data["test2"]["physics_grid"]["plugin_params"]["tau_p1"];
        tau_p2 = data["test2"]["physics_grid"]["plugin_params"]["tau_p2"];
        tau_p3 = data["test2"]["physics_grid"]["plugin_params"]["tau_p3"];
				sigma_ = data["test2"]["physics_grid"]["plugin_params"]["sigma"];
        NUM_REFINEMENTS = data["test2"]["algorithm"]["NUM_REFINEMENTS"];
        maxlevel = data["test2"]["algorithm"]["maxlevel"];
				return;
			}

			int	uniform_refinement (tmesh_3d::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }

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
  					int retval = currentlevel - NUM_REFINEMENTS;
  					for (int ii = 0; ii < 8; ++ii)
    					{     
	      					zcoord = q->p(2, ii);

      						if (fabs(zcoord - 0.0005) < 1e-9)
        					{
	          					retval = 0;
          						break;
        					}
    					}

					if (currentlevel <= NUM_REFINEMENTS)
		    			retval = 0;
      
  					return (retval);
				}
				
			double epsilon_fun (const double & x, const double & y, const double & z) const
				{return z < 0.0005 ? epsilon_0 * epsilon_inf_1 : epsilon_0 * epsilon_inf_2;}

			double csi_1_fun (const double & x, const double & y, const double & z) const
				{return csi1;}

			double csi_2_fun (const double & x, const double & y, const double & z) const
				{return csi2;}

			double csi_3_fun (const double & x, const double & y, const double & z) const
				{return csi3;}

      double tau_p1_fun (const double & x, const double & y, const double & z) const
				{return tau_p1;}

			double tau_p2_fun (const double & x, const double & y, const double & z) const
				{return tau_p2;}

			double tau_p3_fun (const double & x, const double & y, const double & z) const
				{return tau_p3;}

			double sigma_fun (const double & x, const double & y, const double & z, const double & DT) const
				{return sigma_ * DT;}
  };
}

#endif