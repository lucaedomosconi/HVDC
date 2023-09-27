#ifndef VOLTAGES_HPP
#define VOLTAGES_HPP
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
#include "voltage_factory.h"
#include "generic_voltage.h"

using json = nlohmann::json;


namespace voltages{
	class voltage1 : public generic_voltage {
    private:
      double T_discharge;
      double tau;
    public:
      void import_params (const std::string &test_name, const json &data) {
        T_discharge = data[test_name]["algorithm"]["voltage_plugins_params"]["T_discharge"];
        tau = data[test_name]["algorithm"]["voltage_plugins_params"]["tau"];
      }
      double V_in_time (int contact, double time, double x, double y, double z) {
        if (contact == 5)
          return time < T_discharge ? 1.5e4 * (1 - exp(-time/tau)) : 0.0;
        return 0.0;
      }
  };
}

#endif