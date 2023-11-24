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
      void import_params (json &data) {
        try{T_discharge = data["T_discharge"];}
        catch(...){throw std::runtime_error("[T_discharge]");}
        try{tau = data["tau"];}
        catch(...){throw std::runtime_error("[tau]");}
      }
      double V_in_time (int contact, double time, double x, double y, double z) {
        if (contact == 5)
          return time < T_discharge+1e-10 ? 1.5e4 * (1 - exp(-time/tau)) : 0.0;
        return 0.0;
      }
  };

  
	class voltage2 : public generic_voltage {
    private:
      double T_discharge;
      double tau;
    public:
      void import_params (json &data) {
        try{T_discharge = data["T_discharge"];}
        catch(...){throw std::runtime_error("[T_discharge]");}
        try{tau = data["tau"];}
        catch(...){throw std::runtime_error("[tau]");}
      }
      double V_in_time (int contact, double time, double x, double y, double z) {
        if (contact == 5)
          return time < T_discharge+1e-10 ? 1.5e4 * (1 - exp(-time/tau)) : 1.5e4 * (1 - exp(-T_discharge/tau)) * exp(-time/0.5);
        return 0.0;
      }
  };
}

#endif