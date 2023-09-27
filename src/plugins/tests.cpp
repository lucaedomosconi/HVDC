#include "tests.h"
#include "test_factory.h"
#include <iostream>

namespace tests {

namespace {
  struct LoadTest {
    LoadTest() {
      tests::T_factory["test1"] = []() {return std::make_unique<test1>();};
      tests::T_factory["test2"] = []() {return std::make_unique<test2>();};
      tests::T_factory["test3"] = []() {return std::make_unique<test3>();};
      }
    };
  const LoadTest loadtest;
}

}
