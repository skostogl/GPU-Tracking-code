# location of the Python header files

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

# location of the Boost Python include files and library

BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

# compile mesh classes
NAFF.so: interface_NAFF.o
	g++ -O3 -shared -Wl,--export-dynamic interface_NAFF.o -L$(BOOST_LIB) -lboost_python -lfftw3 -lm -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION) -o NAFF.so

interface_NAFF.o: interface_NAFF.cpp NAFF.h
	g++ -O3 -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c interface_NAFF.cpp
