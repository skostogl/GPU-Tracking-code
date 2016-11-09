# location of the Python header files

PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

# location of the Boost Python include files and library

BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

# location of CUDA files

CUDA_PATH=/opt/cuda
CUDA_HEADERS=elements complex float_type

# compile mesh classes
tracker.so: interface.o
	g++ -shared -Wl,--export-dynamic interface.o -L$(BOOST_LIB) -lboost_python -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION) -L $(CUDA_PATH)/lib64 -lnvrtc -lcuda -lcudart -Wl,-rpath,$(CUDA_PATH)/lib64 -o tracker.so

interface.o: interface.cpp $(wildcard *.h)
	$(foreach h,$(CUDA_HEADERS), xxd -i < $(h).h > $(h).xxd && echo ', 0' >> $(h).xxd;)
	g++ -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -I$(CUDA_PATH)/include -fPIC -c interface.cpp
	$(foreach h,$(CUDA_HEADERS), rm $(h).xxd;)
