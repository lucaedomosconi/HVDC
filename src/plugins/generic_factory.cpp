#include "generic_factory.h"

using testfactory = Factory<tests::generic_test, std::function<std::unique_ptr<tests::generic_test>()>>;
using voltagefactory = Factory<voltages::generic_voltage, std::function<std::unique_ptr<voltages::generic_voltage>()>>;

template
testfactory & testfactory::Instance();

template
voltagefactory & voltagefactory::Instance();

