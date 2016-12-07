include ../Makefile.common

NAFF.so: interface_NAFF.o
	echo "Linking NAFF..."
	$(CXX) $(OPTS) -shared -Wl,--export-dynamic interface_NAFF.o $(BOOST_LDFLAGS) -lboost_python -lfftw3 -lm $(PYTHON_LDFLAGS) -lpython$(PYTHON_VERSION) -o NAFF.so

interface_NAFF.o: interface_NAFF.cpp NAFF.h
	echo "Compiling NAFF..."
	$(CXX) $(OPTS) $(PYTHON_CPPFLAGS) $(BOOST_CPPFLAGS) -fPIC -c interface_NAFF.cpp

clean:
	rm -f interface_NAFF.o NAFF.so
