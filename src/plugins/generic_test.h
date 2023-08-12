#ifndef GENERIC_TEST_HPP
#define GENERIC_TEST_HPP
#include <tmesh.h>
#include <simple_connectivity_2d.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
using json = nlohmann::json;


const int NUM_REFINEMENTS = 4;
double T;
double tau;
double tau_p1, tau_p2, tau_p3;
bool save_sol;
// Problem parameters
double epsilon_0;
double epsilon_inf;       // permittivity at infinite frequency
double csi1, csi2, csi3;
double sigma_;            // conducivity coeff

namespace tests{
	class generic_test {
    public:
			virtual void import_params(const json &data) const {return;}
			
      virtual bool works() = 0;

			bool extra_refinement;

			virtual int	uniform_refinement (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }

			virtual int refinement (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }
			
			virtual int coarsening (tmesh::quadrant_iterator q) const
				{ return NUM_REFINEMENTS; }

			virtual double epsilon_fun(const double & x, const double & y) const
				{return epsilon_0 * epsilon_inf;}

			virtual double epsilon_inf_fun(const double & x, const double & y) const
				{return epsilon_inf;}

			virtual double csi_1_fun(const double & x, const double & y) const
				{return csi1;}

			virtual double csi_2_fun(const double & x, const double & y) const
				{return csi2;}

			virtual double csi_3_fun(const double & x, const double & y) const
				{return csi3;}

			virtual double sigma_fun(const double & x, const double & y, const double & DT) const
				{return sigma_ * DT;}

  };
}

#endif