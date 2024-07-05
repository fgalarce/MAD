include ../makefile_header

compile_lib:
ifeq ($(UNAME_S), Darwin)
	$(MAKE) -dynamiclib -o $(ASSIMILATION)/libassimilation.dylib *.cpp $(LINK_EXTERNAL)\
	 -L$(TOOLS) -ltools\
	 -L$(IO) -lio\
	 -L$(SIMPLEX) -lsimplex\
	 -L$(FEM) -lfem\
	 -L$(LINEAR_ALGEBRA) -llinear_algebra\
	 -L$(CALCULUS) -lcalculus
endif

ifeq ($(UNAME_S), Linux)
	$(MAKE) -Wall -fPIC -c *.cpp
	$(MAKE) -shared -o libassimilation.so *.o
	rm *.o
endif

clean:
ifeq ($(UNAME_S), Darwin)
	rm *.o *.dylib
endif

ifeq ($(UNAME_S), Linux)
	rm *.o *.so
endif
