#ifndef VOLTAGE_FACTORY_HPP
#define VOLTAGE_FACTORY_HPP

#include <string>
#include <functional>
#include <map>
#include <memory>
#include <dlfcn.h>
#include "generic_voltage.h"

namespace voltages{
  using voltage_builder = std::function<std::unique_ptr<generic_voltage>()>;
  using voltage_id = std::string;
  using voltage_factory = std::map<voltage_id, voltage_builder>;
  extern voltage_factory V_factory;
}


#endif