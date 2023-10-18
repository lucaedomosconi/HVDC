#ifndef GENERIC_TEST_HPP
#define GENERIC_TEST_HPP
#include <tmesh_3d.h>
#include <simple_connectivity_3d.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
using json = nlohmann::json;


// Problem parameters
double epsilon_0;


namespace tests{
  class generic_test {
    public:
      virtual void import_params(json & data) = 0; // Necessary to pass by non const reference to let exceptions work fine

      bool extra_refinement;

      virtual int	uniform_refinement (tmesh_3d::quadrant_iterator q) const = 0;

      virtual int refinement (tmesh_3d::quadrant_iterator q) const = 0;
      
      virtual int coarsening (tmesh_3d::quadrant_iterator q) const = 0;

      virtual double epsilon_fun(double x, double y, double z) const = 0;

      virtual double csi_1_fun(double x, double y, double z) const = 0;

      virtual double csi_2_fun(double x, double y, double z) const = 0;

      virtual double csi_3_fun(double x, double y, double z) const = 0;

      virtual double tau_p1_fun(double x, double y, double z) const = 0;

      virtual double tau_p2_fun(double x, double y, double z) const = 0;

      virtual double tau_p3_fun(double x, double y, double z) const = 0;

      virtual double sigma_fun(double x, double y, double z, double DT) const = 0;

  };
}

#endif