#include "paper_layer.h"
#include "generic_factory.h"
#include <iostream>
#include <functional>
namespace tests {

namespace {
  struct LoadTest {
    LoadTest() {
      using testfactory = Factory<generic_test, std::function<std::unique_ptr<generic_test>()>>;
      testfactory & T_factory = testfactory::Instance();
      T_factory.add("paper_layer", []() {return std::make_unique<paper_layer>();});
      }
    };
  const LoadTest loadtest;
}

}
