#ifndef PARTICULAR_TEST_HPP
#define PARTICULAR_TEST_HPP
#include "factory.h"
#include "generic_test.h"
using json = nlohmann::json;

namespace tests{

  class test1 : public generic_test {
    public:

			test1() {extra_refinement = false;}
			
			void import_params(const json &data) const {
				epsilon_inf = data["test1"]["epsilon_inf"];
				csi1 = data["test1"]["csi1"];
				csi2 = data["test1"]["csi2"];
				csi3 = data["test1"]["csi3"];
				sigma_ = data["test1"]["sigma_"];
				return;
			}

      bool works() override {return true;} 

			int	uniform_refinement (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS+1; }

			int refinement (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS+1; }

			int coarsening (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS+1; }
				
			double epsilon_fun(const double & x, const double & y) const
				{return epsilon_0 * epsilon_inf;}

			double epsilon_inf_fun(const double & x, const double & y) const
				{return epsilon_inf;}

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