include ../Makefile.common

# files required to compile the cuda kernels
CUDA_HEADERS = elements.h complex.h float_type.h particle.h
CUDA_NAMES = $(basename $(CUDA_HEADERS))

# compile mesh classes
tracker.so: interface.o
	echo "Linking tracker..."
	$(CXX) $(OPTS) -shared -Wl,--export-dynamic interface.o $(BOOST_LDFLAGS) -lboost_python $(PYTHON_LDFLAGS) -lpython$(PYTHON_VERSION) $(CUDA_LDFLAGS) -lnvrtc -lcuda -lcudart -Wl,-rpath,$(CUDA_LDFLAGS) -o tracker.so

interface.o: interface.cpp $(wildcard *.h)
	echo "Compiling tracker..."
	$(foreach h,$(CUDA_NAMES), xxd -i < $(h).h > $(h).xxd && echo ', 0' >> $(h).xxd;)
	$(CXX) $(OPTS) $(PYTHON_CPPFLAGS) $(BOOST_CPPFLAGS) $(CUDA_CPPFLAGS) -fPIC -c interface.cpp
	$(foreach h,$(CUDA_NAMES), rm $(h).xxd;)

clean:
	rm -f interface.o tracker.so
