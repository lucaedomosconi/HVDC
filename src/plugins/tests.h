#ifndef PARTICULAR_TEST_HPP
#define PARTICULAR_TEST_HPP
#include "factory.h"
#include "generic_test.h"
using json = nlohmann::json;


double epsilon_inf_1;       // permittivity at infinite frequency
double csi1, csi2, csi3;
double sigma_;            // conducivity coeff
double sigma_0;
double sigma_1;
double a_bar;


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

			double sigma_fun(const double & x, const double & y, const double & a_local, const double & grad_phi_module, const double & DT) const
				{return sigma_ * DT;}
  };

  class test2 : public generic_test {
    public:

			test2() {extra_refinement = false;}
			
			void import_params(const json &data) const {
				epsilon_inf_1 = data["test2"]["epsilon_inf_1"];
				csi1 = data["test2"]["csi1"];
				csi2 = data["test2"]["csi2"];
				csi3 = data["test2"]["csi3"];
				sigma_0 = data["test2"]["sigma_0"];
				sigma_1 = data["test2"]["sigma_1"];
				E_tr = data["test2"]["E_tr"];
				a_bar = data["test2"]["a_bar"];
				tau_a = data["test2"]["tau_a"];
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

			double sigma_fun(const double & x, const double & y, const double & a_local, const double & grad_phi_module, const double & DT) const
				{
					return ((1-a_local)*sigma_1 + a_local*sigma_0) * DT;
				}
  };
}

#endif