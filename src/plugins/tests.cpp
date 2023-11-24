#include "tests.h"
#include "generic_factory.h"
#include <iostream>
#include <functional>
namespace tests {

namespace {
  struct LoadTest {
    LoadTest() {
      using testfactory = Factory<generic_test, std::function<std::unique_ptr<generic_test>()>>;
      testfactory & T_factory = testfactory::Instance();
      T_factory.add("homogeneous", []() {return std::make_unique<homogeneous>();});
      T_factory.add("two_phase_serial", []() {return std::make_unique<two_phase_serial>();});
      T_factory.add("two_phase_parallel", []() {return std::make_unique<two_phase_parallel>();});
      T_factory.add("hole", []() {return std::make_unique<hole>();});
      }
    };
  const LoadTest loadtest;
}

}
