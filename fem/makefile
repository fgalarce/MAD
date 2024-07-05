include ../makefile_header

compile_lib:
ifeq ($(UNAME_S), Darwin)
	$(MAKE) -dynamiclib -o $(FEM)/libfem.dylib *.cpp $(LINK_EXTERNAL)\
	 -L$(SIMPLEX) -lsimplex\
	 -L$(IO) -lio\
	 -L$(TOOLS) -ltools
endif

ifeq ($(UNAME_S), Linux)
	$(MAKE) -Wall -fPIC -c *.cpp
	$(MAKE) -shared -o libfem.so *.o
	rm *.o
endif

clean:
ifeq ($(UNAME_S), Darwin)
	rm *.o *.dylib
endif

ifeq ($(UNAME_S), Linux)
	rm *.o *.so
endif
