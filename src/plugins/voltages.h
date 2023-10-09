#ifndef VOLTAGES_HPP
#define VOLTAGES_HPP

#include "generic_voltage.h"

using json = nlohmann::json;


namespace voltages{
	class voltage1 : public generic_voltage {
    private:
      double T_discharge;
      double tau;
    public:
      void import_params (const std::string &test_name, json &data) {
        try{T_discharge = data[test_name]["algorithm"]["voltage_plugin_params"]["T_discharge"];}
        catch(...){std::cerr << "Error: Impossible to read object [" << test_name << "][algorithm][voltage_plugin_params][T_discharge]" << std::endl; throw;}
        try{tau = data[test_name]["algorithm"]["voltage_plugin_params"]["tau"];}
        catch(...){std::cerr << "Error: Impossible to read object [" << test_name << "][algorithm][voltage_plugin_params][tau]" << std::endl; throw;}
      }
      double V_in_time (int contact, double time, double x, double y, double z) {
        if (contact == 5)
          return time < T_discharge+1e-10 ? 1.5e4 * (1 - exp(-time/tau)) : 0.0;
        return 0.0;
      }
  };
}

#endif