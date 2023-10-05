#include "voltages.h"
#include "voltage_factory.h"
#include <iostream>

namespace voltages {

namespace {
  struct Loadvoltage {
    Loadvoltage() {
      voltages::V_factory["voltage1"] = []() {return std::make_unique<voltage1>();};
    }
  };
  const Loadvoltage loadvoltage;
}

}
