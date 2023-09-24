#ifndef GENERIC_VOLTAGE_HPP
#define GENERIC_VOLTAGE_HPP
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>
using json = nlohmann::json;



namespace voltages{
	class generic_voltage {
    public:
      virtual void import_params (const std::string &test_name, const json &data) = 0;
      virtual double V_in_time (int contact, double time, double x, double y, double z) = 0;
  };
}

#endif