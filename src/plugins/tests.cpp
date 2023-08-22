#include "tests.h"
#include "factory.h"
#include <iostream>

namespace tests {

namespace {
    struct LoadTest {
        LoadTest() {
            tests::factory["test1"] = []() {return std::make_unique<test1>();};
            std::cout << "factory loaded" << std::endl;
        }
    };
    const LoadTest loadtest;
}

}
