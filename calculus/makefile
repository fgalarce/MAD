include ../makefile_header

compile_lib:
ifeq ($(UNAME_S), Darwin)
	$(MAKE) -dynamiclib -o $(CALCULUS)/libcalculus.dylib *.cpp $(LINK_EXTERNAL)\
	 -L$(SIMPLEX) -lsimplex\
	 -L$(FEM) -lfem\
	 -L$(IO) -lio\
	 -L$(TOOLS) -ltools
endif

ifeq ($(UNAME_S), Linux)
	$(MAKE) -Wall -fPIC -c *.cpp
	$(MAKE) -shared -o libcalculus.so *.o
	rm *.o
endif

clean:
ifeq ($(UNAME_S), Darwin)
	rm *.o *.dylib
endif

ifeq ($(UNAME_S), Linux)
	rm *.o *.so
endif
