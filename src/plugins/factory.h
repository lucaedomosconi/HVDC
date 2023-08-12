#ifndef TEST_FACTORY_HPP
#define TEST_FACTORY_HPP

#include <string>
#include <functional>
#include <map>
#include <memory>
#include <dlfcn.h>
#include "generic_test.h"

namespace tests{
    using test_builder = std::function<std::unique_ptr<generic_test>()>;
    using test_id = std::string;
    using test_factory = std::map<test_id, test_builder>;
    extern test_factory factory;
}


#endif