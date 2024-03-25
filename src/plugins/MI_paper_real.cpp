#include "MI_paper_real.h"
#include "generic_factory.h"
#include <iostream>
#include <functional>
namespace tests {

namespace {
  struct LoadTest {
    LoadTest() {
      using testfactory = Factory<generic_test, std::function<std::unique_ptr<generic_test>()>>;
      testfactory & T_factory = testfactory::Instance();
      T_factory.add("MI_paper_real", []() {return std::make_unique<MI_paper_real>();});
      }
    };
  const LoadTest loadtest;
}

}
