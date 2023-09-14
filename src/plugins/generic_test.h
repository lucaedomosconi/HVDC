#ifndef GENERIC_TEST_HPP
#define GENERIC_TEST_HPP
#include <tmesh_3d.h>
#include <simple_connectivity_3d.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
using json = nlohmann::json;


extern const int NUM_REFINEMENTS;
double T;
double T_discharge;
double tau;
double tau_p1, tau_p2, tau_p3;
bool save_sol;
// Problem parameters
double epsilon_0;


namespace tests{
	class generic_test {
    public:
			virtual void import_params(const json &data) const = 0;

			bool extra_refinement;

			virtual int	uniform_refinement (tmesh_3d::quadrant_iterator q) const = 0;

			virtual int refinement (tmesh_3d::quadrant_iterator q) const = 0;
			
			virtual int coarsening (tmesh_3d::quadrant_iterator q) const = 0;

			virtual double epsilon_fun(const double & x, const double & y, const double & z) const = 0;

			virtual double csi_1_fun(const double & x, const double & y, const double & z) const = 0;

			virtual double csi_2_fun(const double & x, const double & y, const double & z) const = 0;

			virtual double csi_3_fun(const double & x, const double & y, const double & z) const = 0;

			virtual double tau_p1_fun(const double & x, const double & y, const double & z) const = 0;

			virtual double tau_p2_fun(const double & x, const double & y, const double & z) const = 0;

			virtual double tau_p3_fun(const double & x, const double & y, const double & z) const = 0;

			virtual double sigma_fun(const double & x, const double & y, const double & z, const double & DT) const = 0;

  };
}

#endif