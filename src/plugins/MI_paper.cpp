#include "MI_paper.h"
#include "generic_factory.h"
#include <iostream>
#include <functional>
namespace tests {

namespace {
  struct LoadTest {
    LoadTest() {
      using testfactory = Factory<generic_test, std::function<std::unique_ptr<generic_test>()>>;
      testfactory & T_factory = testfactory::Instance();
      T_factory.add("MI_paper", []() {return std::make_unique<MI_paper>();});
      }
    };
  const LoadTest loadtest;
}

}
