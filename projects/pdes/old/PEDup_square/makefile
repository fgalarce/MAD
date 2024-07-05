include ../../../makefile_header

project_name = poro

all: $(project_name).cpp
ifeq ($(UNAME_S), Linux)
	$(MAKE) -o $(project_name).exe $(project_name).cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -L$(CFD) -lcfd -L$(ASSIMILATION) -lassimilation -L$(FEM) -lfem -L$(SIMPLEX) -lsimplex -L$(IO) -lio -L$(TOOLS) -ltools -L$(LINEAR_ALGEBRA) -llinear_algebra -L$(PDE) -lpde
endif
ifeq ($(UNAME_S), Darwin)
	$(MAKE) -o $(project_name).exe $(project_name).cpp -L$(PETSC_LIB) -lpetsc -lslepc -L$(MPI_LIB) -lmpi -headerpad_max_install_names -L$(CFD) -lcfd -L$(ASSIMILATION) -lassimilation -L$(PDE) -lpde -L$(FEM) -lfem -L$(SIMPLEX) -lsimplex -L$(IO) -lio -L$(TOOLS) -ltools -L$(LINEAR_ALGEBRA) -llinear_algebra
	install_name_tool -change libcfd.dylib $(MAD_ROOT)/cfd/libcfd.dylib ./$(project_name).exe
	install_name_tool -change libassimilation.dylib $(MAD_ROOT)/assimilation/libassimilation.dylib ./$(project_name).exe
	install_name_tool -change libpde.dylib $(MAD_ROOT)/pde/libpde.dylib ./$(project_name).exe
	install_name_tool -change libfem.dylib $(MAD_ROOT)/fem/libfem.dylib ./$(project_name).exe
	install_name_tool -change libsimplex.dylib $(MAD_ROOT)/simplex/libsimplex.dylib ./$(project_name).exe
	install_name_tool -change libio.dylib $(MAD_ROOT)/io/libio.dylib ./$(project_name).exe
	install_name_tool -change libtools.dylib $(MAD_ROOT)/tools/libtools.dylib ./$(project_name).exe
	install_name_tool -change liblinear_algebra.dylib $(MAD_ROOT)/la/liblinear_algebra.dylib ./$(project_name).exe
endif

clean:
	rm *.o *.exe
