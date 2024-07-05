include $(MAD_ROOT)/makefile_header

project_name = normals_split_2D

all: $(project_name).cpp
ifeq ($(UNAME_S), Linux)
	$(MAKE) -o $(project_name).exe $(project_name).cpp $(LINK)
endif
ifeq ($(UNAME_S), Darwin)
	$(MAKE) -o $(project_name).exe $(project_name).cpp $(LINK_EXTERNAL) $(LINK)
endif

clean:
	rm *.o *.exe
