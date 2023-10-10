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
      T_factory.add("test1", []() {return std::make_unique<test1>();});
      T_factory.add("test2", []() {return std::make_unique<test2>();});
      T_factory.add("test3", []() {return std::make_unique<test3>();});
      }
    };
  const LoadTest loadtest;
}

}
