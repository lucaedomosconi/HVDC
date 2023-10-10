#include "voltages.h"
#include "generic_factory.h"
#include <iostream>

namespace voltages {

namespace {
  struct Loadvoltage {
    Loadvoltage() {
      using voltagefactory = Factory<generic_voltage, std::function<std::unique_ptr<generic_voltage>()>>;
      voltagefactory & V_factory = voltagefactory::Instance();
      V_factory.add("voltage1", []() {return std::make_unique<voltage1>();});
    }
  };
  const Loadvoltage loadvoltage;
}

}
