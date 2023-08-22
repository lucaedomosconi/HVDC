#ifndef PARTICULAR_TEST_HPP
#define PARTICULAR_TEST_HPP
#include "factory.h"
#include "generic_test.h"
using json = nlohmann::json;


double epsilon_inf_1;       // permittivity at infinite frequency
double csi1, csi2, csi3;
double sigma_;            // conducivity coeff


namespace tests{

  class test1 : public generic_test {
    public:

			test1() {extra_refinement = false;}
			
			void import_params(const json &data) const {
				epsilon_inf_1 = data["test1"]["epsilon_inf_1"];
				csi1 = data["test1"]["csi1"];
				csi2 = data["test1"]["csi2"];
				csi3 = data["test1"]["csi3"];
				sigma_ = data["test1"]["sigma_"];
				return;
			}

      bool works() const {return true;} 

			int	uniform_refinement (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }

			int refinement (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }

			int coarsening (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }
				
			double epsilon_fun(const double & x, const double & y) const
				{return epsilon_0 * epsilon_inf_1;}

			double epsilon_inf_1_fun(const double & x, const double & y) const
				{return epsilon_inf_1;}

			double csi_1_fun(const double & x, const double & y) const
				{return csi1;}

			double csi_2_fun(const double & x, const double & y) const
				{return csi2;}

			double csi_3_fun(const double & x, const double & y) const
				{return csi3;}

			double sigma_fun(const double & x, const double & y, const double & DT) const
				{return sigma_ * DT;}
  };
}

#endif