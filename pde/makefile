include ../makefile_header

compile_lib:
ifeq ($(UNAME_S), Darwin)
	$(MAKE) -dynamiclib -o $(PDE)/libpde.dylib *.cpp $(LINK_EXTERNAL)\
	 -L$(FEM) -lfem\
	 -L$(SIMPLEX) -lsimplex\
	 -L$(TOOLS) -ltools\
	 -L$(CALCULUS) -lcalculus
endif

ifeq ($(UNAME_S), Linux)
	$(MAKE) -Wall -fPIC -c *.cpp
	$(MAKE) -shared -o libpde.so *.o
	rm *.o
endif

clean:
ifeq ($(UNAME_S), Darwin)
	rm *.o *.dylib
endif

ifeq ($(UNAME_S), Linux)
	rm *.o *.so
endif
