#ifndef GENERIC_FACTORY_H
#define GENERIC_FACTORY_H


#include <dlfcn.h>
#include <string>
#include <map>
#include <memory>
#include "generic_test.h"
#include "generic_voltage.h"


template
<class object, class builder>
class Factory {
private:
  Factory() = default;
  Factory(Factory const &) = delete;
  Factory & operator = (Factory const &) = delete;
  std::map<std::string, builder> Container;

public:
  void add (std::string const & identifier, builder const & Builder_fun) {
    Container[identifier] = Builder_fun;
    return;
  }
  std::unique_ptr<object> create (std::string const & identifier) const {
    auto which_test = Container.find(identifier);
    return which_test->second();
  }
  static Factory & Instance();
  
};

#endif