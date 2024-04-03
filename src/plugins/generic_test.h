#ifndef GENERIC_TEST_HPP
#define GENERIC_TEST_HPP
#include <tmesh_3d.h>
#include <connectivity_mi_paper.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
#include <cmath>
using json = nlohmann::json;


// Problem parameters
double epsilon_0;


namespace tests{
  class generic_test {
    public:
      virtual void import_params(json & data) = 0;

      bool extra_refinement;

      virtual int	uniform_refinement (tmesh_3d::quadrant_iterator q) const = 0;

      virtual int refinement (tmesh_3d &tmsh) const = 0;
      
      virtual int coarsening (tmesh_3d::quadrant_iterator q) const = 0;

      virtual double epsilon_fun(double x, double y, double z) const = 0;

      virtual double chi_1_fun(double x, double y, double z) const = 0;

      virtual double chi_2_fun(double x, double y, double z) const = 0;

      virtual double chi_3_fun(double x, double y, double z) const = 0;

      virtual double chi_4_fun(double x, double y, double z) const = 0;

      virtual double chi_5_fun(double x, double y, double z) const = 0;

      virtual double chi_6_fun(double x, double y, double z) const = 0;

      virtual double tau_p1_fun(double x, double y, double z) const = 0;

      virtual double tau_p2_fun(double x, double y, double z) const = 0;

      virtual double tau_p3_fun(double x, double y, double z) const = 0;

      virtual double tau_p4_fun(double x, double y, double z) const = 0;

      virtual double tau_p5_fun(double x, double y, double z) const = 0;

      virtual double tau_p6_fun(double x, double y, double z) const = 0;

      virtual double sigma_fun(double x, double y, double z, double DT, double E = 0) const = 0;

  };
}

#endif