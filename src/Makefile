CXX=mpic++

CXXFLAGS= -std=c++17 -ggdb -O0

CPPFLAGS=-I$(mkP4estInc) -I$(mkMumpsInc) -I$(mkLisInc) -I/u/software/octave_file_io/1.0.91/include/ \
-I/u/software/octave/6.4.0/include/octave-6.4.0 \
-I../include -I../json/include -I../addons \
-I/u/software/bimpp/dev/include -DHAVE_OCTAVE_44 -DOMPI_SKIP_MPICXX -DBIM_TIMING -DUSE_MPI

LDFLAGS=-L/u/software/octave_file_io/1.0.91/lib -Wl,-rpath, -L.\
-L/u/software/bimpp/dev/lib -L$(mkLisLib) -L$(mkMumpsLib) -L$(mkScotchLib) \
-Wl,-rpath,. -ldl

LIBS= -lbim -lbimmumps -lbimlis -lbimp4est -lbimlinalg -llis -ldmumps -lmumps_common \
-lscotcherr -lbz2 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgenericfactory

DOXYFILE= ../doc/Doxyfile

include local_settings.mk

plugs= libgenericfactory.so libVoltages.so libPaperLayer.so libMIPaperReal.so

all : $(plugs) HVDC_main

include HVDC_main.d

libPhysics.so : plugins/tests.cpp plugins/tests.h plugins/generic_test.h
	g++ $(CXXFLAGS) $(CPPFLAGS) -fPIC -c plugins/tests.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -shared tests.o -o libPhysics.so

libMIPaper.so : plugins/MI_paper.cpp plugins/MI_paper.h plugins/generic_test.h
	g++ $(CXXFLAGS) $(CPPFLAGS) -fPIC -c plugins/MI_paper.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -shared MI_paper.o -o libMIPaper.so

libMIPaperReal.so : plugins/MI_paper_real.cpp plugins/MI_paper_real.h plugins/generic_test.h
	g++ $(CXXFLAGS) $(CPPFLAGS) -fPIC -c plugins/MI_paper_real.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -shared MI_paper_real.o -o libMIPaperReal.so

libPaperLayer.so : plugins/paper_layer.cpp plugins/paper_layer.h plugins/generic_test.h
	g++ $(CXXFLAGS) $(CPPFLAGS) -fPIC -c plugins/paper_layer.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -shared paper_layer.o -o libPaperLayer.so

libVoltages.so : plugins/voltages.cpp plugins/voltages.h plugins/generic_voltage.h
	g++ $(CXXFLAGS) $(CPPFLAGS) -fPIC -c plugins/voltages.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -shared voltages.o -o libVoltages.so

libgenericfactory.so : plugins/generic_factory.cpp plugins/generic_factory.h plugins/generic_voltage.h plugins/generic_test.h
	g++ $(CXXFLAGS) $(CPPFLAGS) -fPIC -c plugins/generic_factory.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -shared generic_factory.o -o libgenericfactory.so

%.d : %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $< -MF $@

HVDC_main : HVDC_main.o $(patsubst %.cpp, %.o, $(wildcard *.cpp))
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $? -o $@ $(LIBS)

.PHONY : clean distclean

clean :
	$(RM) *.o
	$(RM) *.d
	$(RM) *.so

distclean : clean
	$(RM) HVDC_main

doc :
	doxygen $(DOXYFILE)
