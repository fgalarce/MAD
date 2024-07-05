include ../makefile_header

create_lib: 
ifeq ($(UNAME_S), Darwin)
	$(MAKE) -dynamiclib -o $(TOOLS)/libtools.dylib *.cpp $(LINK_EXTERNAL)
endif
ifeq ($(UNAME_S), Linux)
	$(MAKE) -Wall -fPIC -c *.cpp
	$(MAKE) -shared -o libtools.so *.o 
	rm *.o
endif
