#include "generic_factory.h"

template <class object, class builder>
Factory<object, builder>& Factory<object, builder>::Instance(){
  static Factory<object, builder> TheFactory;
  std::cout << "instance called" << std::endl;
  return TheFactory;
}

using testfactory = Factory<tests::generic_test, std::function<std::unique_ptr<tests::generic_test>()>>;
using voltagefactory = Factory<voltages::generic_voltage, std::function<std::unique_ptr<voltages::generic_voltage>()>>;

template
testfactory & testfactory::Instance();

template
voltagefactory & voltagefactory::Instance();
