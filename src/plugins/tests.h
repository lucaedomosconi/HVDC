#ifndef PARTICULAR_TEST_HPP
#define PARTICULAR_TEST_HPP
#include "factory.h"
#include "generic_test.h"
using json = nlohmann::json;

const int NUM_REFINEMENTS = 4;
double epsilon_inf_1, epsilon_inf_2;		// permittivity at infinite frequency
double csi1, csi2, csi3;
double sigma_;            					// conducivity coeff
constexpr int maxlevel = 5;

namespace tests{

  class test1 : public generic_test {
    public:

			test1() {extra_refinement = false;}
			
			void import_params(const json &data) const {
				epsilon_inf_1 = data["test1"]["epsilon_inf_1"];
				csi1 = data["test1"]["csi1"];
				csi2 = data["test1"]["csi2"];
				csi3 = data["test1"]["csi3"];
        tau_p1 = data["test1"]["tau_p1"];
        tau_p2 = data["test1"]["tau_p2"];
        tau_p3 = data["test1"]["tau_p3"];
				sigma_ = data["test1"]["sigma_"];
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
    public:

			test2() {extra_refinement = true;}
			
			void import_params (const json &data) const {
				epsilon_inf_1 = data["test2"]["epsilon_inf_1"];
				epsilon_inf_2 = data["test2"]["epsilon_inf_2"];
				csi1 = data["test2"]["csi1"];
				csi2 = data["test2"]["csi2"];
				csi3 = data["test2"]["csi3"];
        tau_p1 = data["test2"]["tau_p1"];
        tau_p2 = data["test2"]["tau_p2"];
        tau_p3 = data["test2"]["tau_p3"];
				sigma_ = data["test2"]["sigma_"];
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