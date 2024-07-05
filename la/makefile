include ../makefile_header

compile_lib:
ifeq ($(UNAME_S), Darwin)
	$(MAKE) -dynamiclib -o $(LINEAR_ALGEBRA)/liblinear_algebra.dylib *.cpp $(LINK_EXTERNAL)\
	 -L$(TOOLS) -ltools\
	 -L$(FEM) -lfem\
	 -L$(SIMPLEX) -lsimplex
endif

ifeq ($(UNAME_S), Linux)
	$(MAKE) -Wall -fPIC -c *.cpp
	$(MAKE) -shared -o liblinear_algebra.so *.o
	rm *.o
endif

clean:
ifeq ($(UNAME_S), Darwin)
	rm *.o *.dylib
endif

ifeq ($(UNAME_S), Linux)
	rm *.o *.so
endif
