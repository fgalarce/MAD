include ../makefile_header

compile_lib:
ifeq ($(UNAME_S), Darwin)
	$(MAKE) -dynamiclib -o $(IO)/libio.dylib *.cpp\
	 -L$(TOOLS) -ltools $(LINK_EXTERNAL) 
endif

ifeq ($(UNAME_S), Linux)
	$(MAKE) -Wall -fPIC -c *.cpp
	$(MAKE) -shared -o libio.so *.o
	rm *.o
endif

clean:
ifeq ($(UNAME_S), Darwin)
	rm *.o *.dylib
endif

ifeq ($(UNAME_S), Linux)
	rm *.o *.so
endif
