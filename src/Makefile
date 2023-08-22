CXX=mpic++

CXXFLAGS= -std=c++17 -ggdb -O0

CPPFLAGS=-I$(mkP4estInc) -I$(mkMumpsInc) -I$(mkLisInc) -I/u/software/octave_file_io/1.0.91/include/ \
-I/u/software/octave/6.4.0/include/octave-6.4.0 \
-I../include -I../json/include -I../addons \
-I/u/software/bimpp/dev/include -DHAVE_OCTAVE_44 -DOMPI_SKIP_MPICXX -DBIM_TIMING -DUSE_MPI

LDFLAGS=-L/u/software/octave_file_io/1.0.91/lib -Wl,-rpath, -L.\
-L/u/software/bimpp/dev/lib -L$(mkLisLib) -L$(mkMumpsLib) -L$(mkScotchLib)

LIBS=-lNanoShaper -lbim -lbimmumps -lbimlis -lbimp4est -lbimlinalg -llis -ldmumps -lmumps_common \
-lscotcherr -lbz2 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lfactory

DOXYFILE= ../doc/Doxyfile

include local_settings.mk

plugs=libfactory.so libTests.so

all : $(plugs) HVDC_main

include HVDC_main.d

libTests.so : libfactory.so plugins/tests.cpp plugins/tests.h plugins/generic_test.h
	g++ $(CXXFLAGS) $(CPPFLAGS) -fPIC -c plugins/tests.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -shared tests.o -o libTests.so

libfactory.so : plugins/factory.cpp plugins/factory.h plugins/generic_test.h
	g++ $(CXXFLAGS) $(CPPFLAGS) -fPIC -c plugins/factory.cpp
	g++ $(CXXFLAGS) $(LDFLAGS) -shared factory.o -o libfactory.so

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
