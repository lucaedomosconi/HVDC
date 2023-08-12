#ifndef PARTICULAR_TEST_HPP
#define PARTICULAR_TEST_HPP
#include "factory.h"
#include "generic_test.h"


namespace tests{

  class test1 : public generic_test {
    public:

      bool works() {return true;} const

			bool extra_refinement = false;

			int	uniform_refinement (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }

			int refinement (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }

			int coarsening (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }
				
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